#ifndef __CS_TURBULENCE_BC_H__
#define __CS_TURBULENCE_BC_H__

/*============================================================================
 * Base turbulence boundary conditions.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * Initialize turbulence model boundary condition ids.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_model_init_bc_ids(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory allocations for turbulence boundary conditions ids.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_model_free_bc_ids(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of \f$ u^\star \f$, \f$ k \f$ and \f$\varepsilon \f$
 *        from a diameter \f$ D_H \f$ and the reference velocity \f$ U_{ref} \f$
 *        for a circular duct flow with smooth wall
 *        (use for inlet boundary conditions).
 *
 * Both \f$ u^\star \f$ and\f$ (k,\varepsilon )\f$ are returned, so that
 * the user may compute other values of \f$ k \f$ and \f$ \varepsilon \f$
 * with \f$ u^\star \f$.
 *
 * We use the laws from Idel'Cik, i.e.
 * the head loss coefficient \f$ \lambda \f$ is defined by:
 * \f[ |\dfrac{\Delta P}{\Delta x}| =
 *                        \dfrac{\lambda}{D_H} \frac{1}{2} \rho U_{ref}^2 \f]
 *
 * then  the relation reads \f$u^\star = U_{ref} \sqrt{\dfrac{\lambda}{8}}\f$.
 * \f$\lambda \f$ depends on the hydraulic Reynolds number
 * \f$ Re = \dfrac{U_{ref} D_H}{ \nu} \f$ and is given by:
 *  - for \f$ Re < 2000 \f$
 *      \f[ \lambda = \dfrac{64}{Re} \f]
 *
 *  - for \f$ Re > 4000 \f$
 *      \f[ \lambda = \dfrac{1}{( 1.8 \log_{10}(Re)-1.64 )^2} \f]
 *
 *  - for \f$ 2000 < Re < 4000 \f$, we complete by a straight line
 *      \f[ \lambda = 0.021377 + 5.3115. 10^{-6} Re \f]
 *
 *  From \f$ u^\star \f$, we can estimate \f$ k \f$ and \f$ \varepsilon\f$
 *  from the well known formulae of developped turbulence
 *
 * \f[ k = \dfrac{u^{\star 2}}{\sqrt{C_\mu}} \f]
 * \f[ \varepsilon = \dfrac{ u^{\star 3}}{(\kappa D_H /10)} \f]
 *
 * \param[in]     uref2         square of the reference flow velocity
 * \param[in]     dh            hydraulic diameter \f$ D_H \f$
 * \param[in]     rho           mass density \f$ \rho \f$
 * \param[in]     mu            dynamic viscosity \f$ \nu \f$
 * \param[out]    ustar2        square of friction speed
 * \param[out]    k             calculated turbulent intensity \f$ k \f$
 * \param[out]    eps           calculated turbulent dissipation
 *                              \f$ \varepsilon \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_ke_hyd_diam(double   uref2,
                             double   dh,
                             double   rho,
                             double   mu,
                             double  *ustar2,
                             double  *k,
                             double  *eps);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of \f$ k \f$ and \f$\varepsilon\f$
 *        from a diameter \f$ D_H \f$, a turbulent intensity \f$ I \f$
 *        and the reference velocity \f$ U_{ref} \f$
 *        for a circular duct flow with smooth wall
 *        (for inlet boundary conditions).
 *
 * \f[ k = 1.5 I {U_{ref}}^2 \f]
 * \f[ \varepsilon = 10 \dfrac{{C_\mu}^{0.75} k^{1.5}}{ \kappa D_H} \f]
 *
 * \param[in]     uref2         square of the reference flow velocity
 * \param[in]     t_intensity   turbulent intensity \f$ I \f$
 * \param[in]     dh            hydraulic diameter \f$ D_H \f$
 * \param[out]    k             calculated turbulent intensity \f$ k \f$
 * \param[out]    eps           calculated turbulent dissipation
 *                               \f$ \varepsilon \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_ke_turb_intensity(double   uref2,
                                   double   t_intensity,
                                   double   dh,
                                   double  *k,
                                   double  *eps);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set inlet boundary condition values for turbulence variables based
 *        on a diameter \f$ D_H \f$ and the reference velocity \f$ U_{ref} \f$
 *        for a circular duct flow with smooth wall
 *        (use for inlet boundary conditions).
 *
 * We use the laws from Idel'Cik, i.e.
 * the head loss coefficient \f$ \lambda \f$ is defined by:
 * \f[ |\dfrac{\Delta P}{\Delta x}| =
 *                        \dfrac{\lambda}{D_H} \frac{1}{2} \rho U_{ref}^2 \f]
 *
 * then  the relation reads \f$u^\star = U_{ref} \sqrt{\dfrac{\lambda}{8}}\f$.
 * \f$\lambda \f$ depends on the hydraulic Reynolds number
 * \f$ Re = \dfrac{U_{ref} D_H}{ \nu} \f$ and is given by:
 *  - for \f$ Re < 2000 \f$
 *      \f[ \lambda = \dfrac{64}{Re} \f]
 *
 *  - for \f$ Re > 4000 \f$
 *      \f[ \lambda = \dfrac{1}{( 1.8 \log_{10}(Re)-1.64 )^2} \f]
 *
 *  - for \f$ 2000 < Re < 4000 \f$, we complete by a straight line
 *      \f[ \lambda = 0.021377 + 5.3115. 10^{-6} Re \f]
 *
 *  From \f$ u^\star \f$, we can estimate \f$ k \f$ and \f$ \varepsilon\f$
 *  from the well known formulae of developped turbulence
 *
 * \param[in]     face_id    boundary face id
 * \param[in]     uref2      square of the reference flow velocity
 * \param[in]     dh         hydraulic diameter \f$ D_H \f$
 * \param[in]     rho        mass density \f$ \rho \f$
 * \param[in]     mu         dynamic viscosity \f$ \nu \f$
 * \param[out]    rcodcl     boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_inlet_hyd_diam(cs_lnum_t   face_id,
                                double      uref2,
                                double      dh,
                                double      rho,
                                double      mu,
                                double     *rcodcl);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set inlet boundary condition values for turbulence variables based
 *        on a diameter \f$ D_H \f$, a turbulent intensity \f$ I \f$
 *        and the reference velocity \f$ U_{ref} \f$
 *        for a circular duct flow with smooth wall.
 *
 * \param[in]     face_id       boundary face id
 * \param[in]     uref2         square of the reference flow velocity
 * \param[in]     t_intensity   turbulent intensity \f$ I \f$
 * \param[in]     dh            hydraulic diameter \f$ D_H \f$
 * \param[out]    rcodcl        boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_inlet_turb_intensity(cs_lnum_t   face_id,
                                      double      uref2,
                                      double      t_intensity,
                                      double      dh,
                                      double     *rcodcl);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set inlet boundary condition values for turbulence variables based
 *        on given k and epsilon values.
 *
 * \param[in]     face_id    boundary face id
 * \param[in]     k          turbulent kinetic energy
 * \param[in]     eps        turbulent dissipation
 * \param[out]    rcodcl     boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_inlet_k_eps(cs_lnum_t   face_id,
                             double      k,
                             double      eps,
                             double     *rcodcl);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set inlet boundary condition values for turbulence variables based
 *        on given k and epsilon values only if not already initialized.
 *
 * \param[in]     face_id    boundary face id
 * \param[in]     k          turbulent kinetic energy
 * \param[in]     eps        turbulent dissipation
 * \param[out]    rcodcl     boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_set_uninit_inlet_k_eps(cs_lnum_t   face_id,
                                        double      k,
                                        double      eps,
                                        double     *rcodcl);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute matrix \f$\tens{\alpha}\f$ used in the computation of the
 * Reynolds stress tensor boundary conditions.
 *
 * We note \f$\tens{R}_g\f$ the Reynolds Stress tensor in the global reference
 * frame (mesh reference frame) and \f$\tens{R}_l\f$ the Reynolds stress
 * tensor in the local reference frame (reference frame associated to the
 * boundary face).
 *
 * \f$\tens{P}_{lg}\f$ is the change of basis orthogonal matrix from local
 * to global reference frame.

 * \f$\tens{\alpha}\f$ is a 6 by 6 matrix such that:
 * \f[
 * \vect{R}_{g,\fib} = \tens{\alpha} \vect{R}_{g,\centip} + \vect{R}_{g}^*
 * \f]
 * where symetric tensors \f$\tens{R}_g\f$ have been unfolded as follows:
 * \f[
 * \vect{R}_g = \transpose{\left(R_{g,11},R_{g,22},R_{g,33},
 *                              R_{g,12},R_{g,13},R_{g,23}\right)}
 * \f].
 *
 * \f$\tens{\alpha}\f$ is defined so that \f$ \tens{R}_{g,\fib} \f$ is computed
 * as a function of \f$\tens{R}_{g,\centip}\f$ as follows:
 * \f[
 * \tens{R}_{g,\fib}=\tens{P}_{lg}\tens{R}_{l,\fib}\transpose{\tens{P}_{lg}}
 * \f]
 *
 * with
 * \f[
 * \tens{R}_{l,\fib} =
 * \begin{bmatrix}
 * R_{l,11,\centip}   &   u^* u_k        & c R_{l,13,\centip}\\
 *   u^* u_k          & R_{l,22,\centip} & 0                 \\
 * c R_{l,13,\centip} & 0                & R_{l,33,\centip}
 * \end{bmatrix} +
 * \underbrace{\begin{bmatrix}
 *                 0  &   u^* u_k        & 0                 \\
 *   u^* u_k          & 0                & 0                 \\
 * 0                  & 0                & 0
 * \end{bmatrix}}_{\tens{R}_l^*}
 * \f]
 *
 * and
 * \f$\tens{R}_{l,\centip}=\transpose{\tens{P}_{lg}}\tens{R}_{g,\centip}
 *                       \tens{P}_{lg}\f$.
 *
 * Constant c is chosen depending on the type of the boundary face:
 * \f$c = 0\f$ at a wall face, \f$c = 1\f$ at a symmetry face.
 *
 * \param[in]      is_sym  Constant c in description above
 *                         (1 at a symmetry face, 0 at a wall face)
 * \param[in]      p_lg    change of basis matrix (local to global)
 * \param[out]     alpha   transformation matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_bc_rij_transform(int        is_sym,
                               cs_real_t  p_lg[3][3],
                               cs_real_t  alpha[][6]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_BC_H__ */
