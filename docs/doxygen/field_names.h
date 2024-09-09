/*============================================================================
 * Doxygen documentation for specific field names
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "cs_equation_param.h"
#include "cs_parameters.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file field_names.h
        General field names
*/

/*----------------------------------------------------------------------------*/

/*!
 * \defgroup field_names Field names
 */
/*!@{*/

/*!
  \var algo

  Fields used to check algorithm behavior.

  Fields of the form "algo:<type>_<variable_name>" are reserved
  by the code to allow visualization of some intermediate computational
  field values.

  If the user creates such a field, and the code calls the matching
  operator, this field will be updated automatically.

  Current reserved fields are:


  - <tt> algo:predicted_velocity </tt>
    Velocity field after the prediction step, cell-based field of dimension 3.
  - <tt> algo:predicted_velocity_divergence </tt>
    Divergence of the velocity field after the prediction step, cell-based field of dimension 1.
  - <tt> algo:pressure_gradient </tt>
    Pressure gradient, cell-based field of dimension 3.
  - <tt> algo:gradient_ </tt> supplemented with a transported field name
    Field gradient, cell-based field of dimension 3x the dimension of the corresponding field.
  - <tt> algo:turbulent_flux_production_ </tt> supplemented with a transported field name
    Field corresponding to the production of u'T' if scalar taken is temperature, cell-based field of dimension 3.
  - <tt> algo:velocity_gradient</tt>
    Velocity gradient, cell-based field of dimension 3x3.
  - <tt> algo:pressure_gradient_increment </tt>
    Gradient of the pressure increment solved in the correction step, cell-based field of dimension 3.
  - <tt> algo:rij_divergence </tt>
    Divergence of the Reynolds stress in the momentum equation for  Reynolds stress RANS models, cell-based field of dimension 3.
  - <tt> algo:rij_production </tt>
    So called production term in Reynolds stress RANS models, cell-based field of dimension 6.
  - <tt> algo:rij_pressure_strain_correlation </tt>
    So called pressure strain correlation term in Reynolds stress RANS models, cell-based field of dimension 6.
  - <tt> algo:rij_buoyancy </tt>
    So called buoyancy term in Reynolds stress RANS models, cell-based field of dimension 6.
  - <tt> algo:k_production </tt>
    So called production term in k-epsilon RANS models, cell-based field of dimension 1.
  - <tt> algo:k_buoyancy </tt>
    So called buoyancy term in k-epsilon RANS models, cell-based field of dimension 1.
  - <tt> algo:turbulent_flux_divergence </tt>
    Divergence of the turbulent flux in the energy equation, cell-based field of dimension 1.

  - <tt> algo:grad_clip_factor_<variable_name> </tt>
    Least-squares gradient clip factor if gradient clipping is activated
    (scalar field on cells).
  - <tt> algo:grad_b_iter_<variable_name> </tt>
    Number of iterations for convergence of fixed-point algorithm for
    least-squares boundary gradient of vector and tensor fields.
    (scalar field on boundary faces).
  - <tt> algo:b_velocity_inout </tt>
    Boundary field to visualize entering velocity at outlets.
*/
char *algo;

/*!
  \var est_error_pre_2

  Error estimator for Navier-Stokes: prediction.

  After the velocity prediction step (yielding \f$\vect{u}^*\f$), the
  estimator \f$\eta^{\,pred}_{\,i,k}(\vect{u}^*)\f$, local variable calculated
  at every cell \f$ \Omega_i \f$, is created from
  \f$\vect{\mathcal R}^{\,pred}(\vect{u}^*)\f$,
  which represents the residual of the equation solved during this step:
  \f$\vect{u}\f$ and \f$ P \f$:
  \f{eqnarray*}{
    \vect{\mathcal R}^{\,pred}(\vect{u}^*)
        & = & \rho^n \dfrac{\vect{u}^*-\vect{u}^n}{\Delta t}
            + \rho^n \vect{u}^n \cdot \gradt (\vect{u}^*)
            - \divv \left((\mu+\mu_t)^n \gradt (\vect{u}^*) \right)
            + \grad(P^n)
    \\  & - & \text{rest of the right-hand member }
             (\vect{u}^n, P^n, \text{other variables}^n)
  \f}
  - By definition:
    \f$ \eta^{\,pred}_{\,i,k}(\vect{u}^*)= {|\Omega_i|}^{\,(k-2)/2}\ ||\vect{\mathcal R}^{\,pred}(\vect{u}^*)|| _{{IL}^{2}(\Omega_i)} \f$
  - The first family, k=1, suppresses the
    volume \f$ |\Omega_i| \f$ which intrinsicly appears  with the norm
    \f$ {IL}^{2}(\Omega_i) \f$.
  - The second family, k=2, exactly represents the norm
    \f$ {IL}^{2}(\Omega_i) \f$. The size of the cell therefore
    appears in its calculation and induces a weighting effect.
  - \f$ \eta^{\,pred}_{\,i,k}(\vect{u}^*)\f$  is ideally equal to zero when the
    reconstruction methods are perfect and the associated system is solved exactly.
*/
char *est_error_pre_2;

/*!
  \var est_error_der_2

  Error estimator for Navier-Stokes: drift.

  The estimator \f$\eta^{\,der}_{\,i,k}(\vect{u}^{\,n+1})\f$ is based on the
  following quantity (intrinsic to the code):
  \f{eqnarray*}{
  \eta^{\,der}_{\,i,k}(\vect{u}^{\,n+1})
     &=& {|\Omega_i|}^{(k-2)/2}
        || \divs (\text{corrected mass flow after the pressure step})
        - \Gamma||_{{L}^{2}(\Omega_i)}
  \\ &=& {|\Omega_i|}^{(1-k)/2}
       | \divs (\text{corrected mass flow after the pressure step})- \Gamma|
  \f}
  - Ideally, it is equal to zero when the Poisson equation related to the
    pressure is solved exactly.
*/
char *est_error_der_2;

/*!
  \var est_error_cor_2

  Error estimator for Navier-Stokes: correction.

  The estimator \f$ \eta^{\,corr}_{\,i,k}(\vect{u}^{\,n+1})\f$ comes directly
  from the mass flow calculated with the updated velocity field:
  \f{eqnarray*}{
  \eta^{\,corr}_{\,i,k}(\vect{u}^{\,n+1})=
  |\Omega_i|^{\,\delta_{\,2,k}}\ |div (\rho^n \vect{u}^{n+1}) - \Gamma|
  \f}
  - The velocities \f$\vect{u}^{n+1}\f$ are taken at the cell centers,
    the divergence is calculated after projection on the faces.
    \f$ \,\delta_{\,2,k}\f$ represents the Kronecker symbol.
  - The first family, k=1, is the absolute raw value of the divergence of the mass flow
    minus the mass source term.
    The second family, $k=2$, represents a physical property and allows to evaluate
    the difference in \f$kg.s^{\,-1}\f$.
  - Ideally, it is equal to zero when the Poisson equation is solved exactly and
    the projection from the mass flux at the faces to the velocity at the cell
    centers is made in a set of  functions with null divergence.
*/
char *est_error_cor_2;

/*!
  \var est_error_tot_2

  Error estimator for Navier-Stokes: total.

  The estimator \f$ \eta^{\,tot}_{\,i,k}(\vect{u}^{\,n+1})\f$, local variable
  calculated at every cell \f$\Omega_i\f$, is based on the quantity
  \f$\vect{\mathcal R}^{\,tot}(\vect{u}^{\,n+1})\f$, which represents the
  residual of the equation using the updated values of
  \f$\vect{u}\f$ and \f$P\f$:
  \f{eqnarray*}{
    \vect{\mathcal R}^{\,pred}(\vect{u}^*)
        & = & \rho^n \dfrac{\vect{u}^*-\vect{u}^n}{\Delta t}
            + \rho^n \vect{u}^n \cdot \gradt (\vect{u}^*)
            - \divv \left((\mu+\mu_t)^n \gradt (\vect{u}^*) \right)
            + \grad(P^n)
    \\  & - & \text{rest of the right-hand member }
             (\vect{u}^n, P^n, \text{other variables}^n)
  \f}
  - By definition:
    \f$ \eta^{\,tot}_{\,i,k}(\vect{u}^{\,n+1})= {|\Omega_i|}^{\,(k-2)/2}\ ||\vect{\mathcal R}^{\,tot}(\vect{u}^{\,n+1})||
    _{{I\hspace{-.25em}L}^{2}(\Omega_i)} \f$
  - The mass flux in the convective term is recalculated from \f$\vect{u}^{n+1}\f$
    expressed at the cell centers (and not taken from the updated mass flow at the
    faces).
  - As for the prediction estimator:
    - The first family, k=1, suppresses the
      volume \f$ |\Omega_i| \f$ which intrinsicly appears  with the norm
      \f$ {IL}^{2}(\Omega_i) \f$.
    - The second family, k=2, exactly represents the norm
      \f$ {IL}^{2}(\Omega_i) \f$. The size of the cell therefore
      appears in its calculation and induces a weighting effect.
*/
char *est_error_tot_2;

/*!@}*/

/*----------------------------------------------------------------------------*/
