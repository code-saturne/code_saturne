#ifndef __CS_VOF_H__
#define __CS_VOF_H__

/*============================================================================
 * Functions associated to VOF model
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
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*!
   @addtogroup vof
   @{

   \defgroup vof_masks Masks used to specify Volume of Fluid models
   @{
 */

/*! Volume of Fluid model */
#define CS_VOF_ENABLED (1 << 0)

/*! Free surface model */
#define CS_VOF_FREE_SURFACE (1 << 1)

/*! Mass transfer Merkle model for vaporization / condensation (cavitation) */
#define CS_VOF_MERKLE_MASS_TRANSFER (1 << 2)

/*!
    @}
    @}
*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* VOF model parameters */
/*----------------------*/

typedef struct {

  unsigned      vof_model;   /* VoF model (sum of masks defining Volume of
                                Fluid model and submodels */

  double        rho1;        /* density */

  double        rho2;

  double        mu1;         /* viscosity */

  double        mu2;

  double        sigma_s;     /* surface tension */

  int           idrift;      /* drift velocity model */

  double        cdrift;      /* C_gamma constante (drift flux factor)*/

  double        kdrift;

} cs_vof_parameters_t;

/* Cavitation parameters */
/*-----------------------*/

typedef struct {

  cs_real_t        presat;  /* reference saturation pressure */
  cs_real_t        uinf;    /* reference velocity */
  cs_real_t        linf;    /* reference length scale */
  cs_real_t        cdest;   /* constant of condensation model (Merkle) */
  cs_real_t        cprod;   /* constant of vaporization model (Merkle) */
  int              icvevm;  /* eddy-viscosity correction indicator */
  cs_real_t        mcav;    /* eddy-viscosity correction cstt (Reboud) */
  int              itscvi;  /* eddy-viscosity correction indicator */

} cs_cavitation_parameters_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/* pointer to VOF model parameters structure */

extern const cs_vof_parameters_t *cs_glob_vof_parameters;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Provide access to VOF structure.
 */
/*----------------------------------------------------------------------------*/

cs_vof_parameters_t *
cs_get_glob_vof_parameters(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the mixture density, mixture dynamic viscosity given fluid
 *         volume fractions and the reference density and dynamic viscosity
 *         \f$ \rho_l, \mu_l \f$ (liquid), \f$ \rho_v, \mu_v \f$ (gas).
 *
 * Computation is done as follows on cells:
 * \f[
 * \rho_\celli = \alpha_\celli \rho_v + (1-\alpha_\celli) \rho_l,
 * \f]
 * \f[
 * \mu_\celli = \alpha_\celli \mu_v + (1-\alpha_\celli) \mu_l,
 * \f]
 *
 * A similar linear formula is followed on boundary using fluid volume fraction
 * value on the boundary.
 *
 * \param[in]  m  pointer to mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_compute_linear_rho_mu(const cs_mesh_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the mixture density, mixture dynamic viscosity and mixture
 *        mass flux given the volumetric flux, the volume fraction and the
 *        reference density and dynamic viscosity \f$ \rho_l, \mu_l \f$
 *        (liquid), \f$ \rho_v, \mu_v \f$ (gas).
 *
 * For the computation of mixture density, mixture dynamic viscosity, see
 * \ref cs_vof_compute_linear_rho_mu.
 *
 * Computation of mass flux is as follows:
 * \f[
 * \left( \rho\vect{u}\cdot\vect{S} \right)_\ij = \\ \left\lbrace
 * \begin{array}{ll}
 *   \rho_\celli (\vect{u}\cdot\vect{S})_\ij
 *  &\text{ if } (\vect{u}\cdot\vect{S})_\ij>0, \\
 *   \rho_\cellj (\vect{u}\cdot\vect{S})_\ij
 *  &\text{ otherwise },
 * \end{array} \right.
 * \f]
 * \f[
 * \left( \rho\vect{u}\cdot\vect{S} \right)_\ib = \\ \left\lbrace
 * \begin{array}{ll}
 *   \rho_\celli (\vect{u}\cdot\vect{S})_\ib
 *  &\text{ if } (\vect{u}\cdot\vect{S})_\ib>0, \\
 *   \rho_b (\vect{u}\cdot\vect{S})_\ib
 *  &\text{ otherwise }.
 * \end{array} \right.
 * \f]
 *
 * \param[in]  m  pointer to mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_update_phys_prop(const cs_mesh_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface tension momentum source term following the CSF
 * model of Brackbill et al. (1992).
 *
 * \param[in]   m    pointer to mesh structure
 * \param[in]   mq   pointer to mesh quantities structure
 * \param[out]  stf  surface tension momentum source term
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_surface_tension(const cs_mesh_t             *m,
                       const cs_mesh_quantities_t  *mq,
                       cs_real_3_t                  stf[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write in main log the global mixture mass budget:
 * \f[
 * \sum_i\left(
 * |\Omega_i|\dfrac{\alpha_i^n - \alpha_i^{n-1}}{\Delta t} +
 * \sum_{j\in\Face{\celli}}\left(\rho\vect{u}\vect{S}\right)_{ij}^n
 * \right).
 * \f]
 *
 * \param[in]  m   pointer to mesh structure
 * \param[in]  mq  pointer to mesh quantities structure
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_log_mass_budget(const cs_mesh_t             *m,
                       const cs_mesh_quantities_t  *mq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the flux of the drift velocity \f$ \vect u _d \f$,
 *        by using the flux of the standard velocity \f$ \vect u \f$
 *        following the approach described by
 *        Suraj S Deshpande et al 2012 Comput. Sci. Disc. 5 014016.
 *        It is activated with the option idrift = 1
 *
 * Using the notation:
 * \f[
 * \begin{cases}
 * \left ( \vect u ^{n+1} . \vect S \right ) _{\face} = \Dot{m}_{\face}\\
 * \left ( \vect u _d^{n+1} . \vect S \right ) _{\face} = \Dot{m^d}_{\face}
 * \end{cases}
 * \f]
 * The drift flux is computed as:
 * \f[
 * \Dot{m^d}_{\face} = min \left ( C_{\gamma} \dfrac{\Dot{m}_{\face}}
 * {\vect S_{\face}}, \underset{\face'}{max} \left [ \dfrac{\Dot{m}_{\face'}}
 * {\vect S_{\face'}} \right ] \right ) \left ( \vect n . \vect S \right )
 * _{\face}
 * \f]
 * Where \f$ C_{\gamma} \f$ is the drift flux factor defined with the variable
 * \ref cdrift, \f$ \vect n _{\face} \f$ the normal vector to the interface.
 * The gradient is computed using a centered scheme:
 * \f[
 * {\vect n _{\face}} = \dfrac{\left ( \grad \alpha \right ) _{\face}}
 * {\norm {\left ( \grad \alpha \right ) _{\face} + \delta}},
 * \text{ with: }
 * \left ( \grad \alpha \right ) _{\face _{\celli \cellj}} = \dfrac{\left (
 * \grad \alpha \right ) _\celli + \left ( \grad \alpha \right ) _\cellj}{2},
 * \text{ and: }
 * \delta = 10^{-8} / \overline{\vol \celli} ^{1/3}
 * \f]
 *
 * \param[in]   m    pointer to mesh structure
 * \param[in]   mq   pointer to mesh quantities structure
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_deshpande_drift_flux(const cs_mesh_t             *m,
                            const cs_mesh_quantities_t  *mq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a
 * standard transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \alpha_\celli^{n+1} \left( 1 - \alpha_\cellj^{n+1} \right) \left(
 *        \dot{m}_\fij^{d} \right)^{+} + \alpha_\cellj^{n+1} \left( 1 -
 *        \alpha_\celli^{n+1} \right) \left( \dot{m}_\fij^{d} \right)^{-}
 *       \right)
 * \f]
 *
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least squares gradient
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 by neighboring gradients
 *                               - = 1 by the mean gradient
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coefficient for the computation of
 *                               the gradient
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_drift_term(int                        imrgra,
                  int                        nswrgp,
                  int                        imligp,
                  int                        iwarnp,
                  cs_real_t                  epsrgp,
                  cs_real_t                  climgp,
                  cs_real_t        *restrict pvar,
                  const cs_real_t  *restrict pvara,
                  cs_real_t        *restrict rhs);

/*----------------------------------------------------------------------------
 *!
 * \brief Provide access to cavitation parameters structure.
 */
/*----------------------------------------------------------------------------*/

cs_cavitation_parameters_t *
cs_get_glob_cavitation_parameters(void);

/*----------------------------------------------------------------------------*/
/*
 * \param[in]     iterns        Navier-Stokes iteration number
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_solve_void_fraction(int  iterns);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_VOF_H__ */
