/*============================================================================
 * Doxygen documentation for specific field names
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

#include "cs_equation_param.h"
#include "cs_parameters.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file field names.h
        General field names
*/

/*----------------------------------------------------------------------------*/

/*!
 * \defgroup field_names Field names
 */
/*!@{*/

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
