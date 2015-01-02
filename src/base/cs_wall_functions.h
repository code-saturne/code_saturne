#ifndef __CS_WALL_FUNCTIONS_H__
#define __CS_WALL_FUNCTIONS_H__

/*============================================================================
 * Wall functions
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/* Wall functions descriptor */
/*--------------------------------------------*/

typedef struct {

  int           ideuch;       /* wall functions
                                 - 0: one scale of friction velocities
                                 - 1: two scale of friction velocities
                                 - 2: scalable wall functions */
  int           iwallt;       /* exchange coefficient correlation
                                 - 0: not use by default
                                 - 1: exchange coefficient computed with a
                                      correlation */
  int           ilogpo;       /* wall function with
                                 - 0: a power lay (deprecated)
                                 - 1: a log lay */
  double        ypluli;       /* limit value of y+ for the viscous
                                 sublayer */

} cs_wall_functions_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Pointer to wall functions descriptor structure */

extern const cs_wall_functions_t *cs_glob_wall_functions;

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_wall_functions_velocity.
 *----------------------------------------------------------------------------*/

void CS_PROCF (wallfunctions, WALLFUNCTIONS)
(
 const cs_int_t   *const iwallf,
 const cs_lnum_t  *const ifac,
 const cs_real_t  *const viscosity,
 const cs_real_t  *const t_visc,
 const cs_real_t  *const vel,
 const cs_real_t  *const y,
 const cs_real_t  *const kinetic_en,
       cs_int_t        *iuntur,
       cs_lnum_t        *nsubla,
       cs_lnum_t        *nlogla,
       cs_real_t        *ustar,
       cs_real_t        *uk,
       cs_real_t        *yplus,
       cs_real_t        *ypup,
       cs_real_t        *cofimp,
       cs_real_t        *dplus
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_wall_functions_scalar.
 *----------------------------------------------------------------------------*/

void CS_PROCF (hturbp, HTURBP)
(
 const cs_real_t  *const prl,
 const cs_real_t  *const prt,
 const cs_real_t  *const ckarm,
 const cs_real_t  *const yplus,
 const cs_real_t  *const dplus,
       cs_real_t        *htur,
       cs_real_t        *yplim
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*! \brief  Compute the friction velocity and \f$y^+\f$ / \f$u^+\f$.

*/
/*-------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 ______________________________________________________________________________*/
/*!
 * \param[in]     iwallf        wall function type
 * \param[in]     ifac          face number
 * \param[in]     l_visc        kinematic viscosity
 * \param[in]     t_visc        turbulent kinematic viscosity
 * \param[in]     vel           wall projected
 *                              cell center velocity
 * \param[in]     y             wall distance
 * \param[in]     kinetic_en    turbulente kinetic energy
 * \param[in]     iuntur        indicator:
 *                              0 in the viscous sublayer
 * \param[in]     nsubla        counter of cell in the viscous
 *                              sublayer
 * \param[in]     nlogla        counter of cell in the log-layer
 * \param[out]    ustar         friction velocity
 * \param[out]    uk            friction velocity
 * \param[out]    yplus         non-dimension wall distance
 * \param[out]    ypup          yplus projected vel ratio
 * \param[out]    cofimp        \f$\frac{|U_F|}{|U_I^p|}\f$ to ensure a good
 *                              turbulence production
 * \param[out]    dplus         dimensionless shift to the wall
 *                              for scalable wall functions
 */
/*-------------------------------------------------------------------------------*/

void
cs_wall_functions_velocity(int          iwallf,
                           cs_lnum_t    ifac,
                           cs_real_t    l_visc,
                           cs_real_t    t_visc,
                           cs_real_t    vel,
                           cs_real_t    y,
                           cs_real_t    kinetic_en,
                           int         *iuntur,
                           cs_lnum_t   *nsubla,
                           cs_lnum_t   *nlogla,
                           cs_real_t   *ustar,
                           cs_real_t   *uk,
                           cs_real_t   *yplus,
                           cs_real_t   *ypup,
                           cs_real_t   *cofimp,
                           cs_real_t   *dplus);

/*-------------------------------------------------------------------------------*/

/*! \brief Compute the correction of the exchange coefficient between the fluid and
  the wall for a turbulent flow.

  This is function of the dimensionless
  distance to the wall \f$ y^+ = \dfrac{\centip \centf u_\star}{\nu}\f$.

  Then the return coefficient reads:
  \f[
  h_{tur} = Pr \dfrac{y^+}{T^+}
  \f]

  This coefficient is computed thanks to a similarity model between
  dynamic viscous sub-layer and themal sub-layer.

  \f$ T^+ \f$ is computed as follows:

  - For a laminar Prandtl number smaller than 0.1 (such as liquid metals),
    the standard model with two sub-layers (Prandtl-Taylor) is used.

  - For a laminar Prandtl number larger than 0.1 (such as liquids and gaz),
    a model with three sub-layers (Arpaci-Larsen) is used.

  The final exchange coefficient is:
  \f[
  h = \dfrac{K}{\centip \centf} h_{tur}
  \f]

*/
/*-------------------------------------------------------------------------------
  Arguments
 ______________________________________________________________________________.
   mode           name          role                                           !
 ______________________________________________________________________________*/
/*!
 * \param[in]     prl           laminar Prandtl number
 * \param[in]     prt           turbulent Prandtl number
 * \param[in]     ckarm         Von Karman constant
 * \param[in]     yplus         dimensionless distance to the wall
 * \param[in]     dplus         dimensionless distance for scalable
 *                              wall functions
 * \param[out]    htur          corrected exchange coefficient
 * \param[out]    yplim         value of the limit for \f$ y^+ \f$
 */
/*-------------------------------------------------------------------------------*/

void
cs_wall_functions_scalar(double  prl,
                         double  prt,
                         double  ckarm,
                         double  yplus,
                         double  dplus,
                         double  *htur,
                         double  *yplim);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_WALL_FUNCTIONS_H__ */
