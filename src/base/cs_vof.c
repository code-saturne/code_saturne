/*============================================================================
 * Functions associated to VOF model
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_boundary_conditions.h"
#include "cs_convection_diffusion.h"
#include "cs_equation.h"
#include "cs_face_viscosity.h"
#include "cs_divergence.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_log.h"
#include "cs_physical_constants.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_rotation.h"
#include "cs_sles_default.h"
#include "cs_turbomachinery.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_vof.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_vof.c
        VOF model data.
*/

/*----------------------------------------------------------------------------*/

/*!

  \defgroup homogeneous_mixture Homogeneous mixture modelling

  @addtogroup homogeneous_mixture
  @{

  \defgroup vof VOF model for free surface flow or dispersed flow

  @addtogroup vof
  @{

  \defgroup vof_mixture_properties Mixture properties

  @addtogroup vof_mixture_properties
  @{

  \struct cs_vof_parameters_t

  \brief VOF model parameters. Void fraction variable tracks fluid 2.

  Members of this structure are publicly accessible to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_vof_parameters_t::vof_model
        Volume of Fluid model - sum of masks defining VoF model and submodels.
        See defined masks in \ref vof_masks.

  \var  cs_vof_parameters_t::rho1
        reference density of fluid 1 (kg/m3).
        By convention, liquid phase for cavitation model.

  \var  cs_vof_parameters_t::rho2
        reference density of fluid 2 (kg/m3).
        By convention, gas phase for cavitation model.

  \var  cs_vof_parameters_t::mu1
        reference molecular viscosity of fluid 1 (kg/(m s))

  \var  cs_vof_parameters_t::mu2
        reference molecular viscosity of fluid 2 (kg/(m s))

  \var  cs_vof_parameters_t::sigmaS
        reference surface tension (N/m)

  \var  cs_vof_parameters_t::idrift
        drift velocity model
            - 0: drift model disable
            - 1: field inner_drift_velocity_flux is used
                 Deshpande's model
            - 2: field drift_velocity is used
                 User-defined drift velocity field

  \var  cs_vof_parameters_t::cdrift
        Flux factor parameter.
        In case of drift flux, factor of the local flux compared to the
        global max flux.

  \var  cs_vof_parameters_t::kdrift
        Turbulent like diffusion effect (m2/s).
        In case of drift velocity, factor of a volume fraction gradient

  @}

  \defgroup cavitation Cavitation model

  @addtogroup cavitation
  @{

  \struct cs_cavitation_parameters_t

  \brief Cavitation model parameters.

  Members of this structure are publicly accessible to allow for concise
  syntax, as it is expected to be used in many places.

  \defgroup cav_source_term Vaporization/condensation model

  @addtogroup cav_source_term
  @{

  \var  cs_cavitation_parameters_t::presat
        Reference saturation pressure (kg/(m s2)).

  \var  cs_cavitation_parameters_t::uinf
        Reference velocity of the flow (m/s).

  \var  cs_cavitation_parameters_t::linf
        Reference length scale of the flow (m).

  \var  cs_cavitation_parameters_t::cdest
        Constant Cdest of the condensation source term (Merkle model).

  \var  cs_cavitation_parameters_t::cprod
        Constant Cprod of the vaporization source term (Merkle model).

  @}

  \defgroup cav_turbulence Interaction with turbulence

  @addtogroup cav_turbulence
  @{

  \var  cs_cavitation_parameters_t::icvevm
        Activation of the eddy-viscosity correction (Reboud correction).
            - 1: activated
            - 0: desactivated

  \var  cs_cavitation_parameters_t::mcav
        Constant mcav of the eddy-viscosity correction (Reboud correction).

  @}

  \defgroup cav_numerics Numerical parameters

  @addtogroup cav_numerics
  @{

  \var  cs_cavitation_parameters_t::itscvi
        Implicitation in pressure of the vaporization/condensation model
           - 1: activated
           - 0: desactivated

  @}

  @}

  @}

  @}

*/
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static variables
 *============================================================================*/

static cs_vof_parameters_t  _vof_parameters =
{
  .vof_model     = 0,
  .rho1          = 1.e3,
  .rho2          = 1.,
  .mu1           = 1.e-3,
  .mu2           = 1.e-5,
  .sigmaS        = 0.,
  .idrift        = 0,
  .cdrift        = 1.,
  .kdrift        = 0.
};

static cs_cavitation_parameters_t  _cavit_parameters =
{
  .presat =  2.e3,
  .uinf   =  -1e13,
  .linf   =  1.e-1,
  .cdest  =  5.e1,
  .cprod  =  1.e4,
  .icvevm =  1,
  .mcav   =  1.e1,
  .itscvi =  1
};

/*============================================================================
 * Global variables
 *============================================================================*/

const cs_vof_parameters_t *cs_glob_vof_parameters = &_vof_parameters;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_vof_get_pointers(unsigned **ivofmt,
                      double   **rho1,
                      double   **rho2,
                      double   **mu1,
                      double   **mu2,
                      double   **sigmaS,
                      int      **idrift,
                      double   **cdrift,
                      double   **kdrift);

void
cs_f_vof_compute_linear_rho_mu(void);

void
cs_f_vof_update_phys_prop(void);

void
cs_f_vof_surface_tension(cs_real_3_t stf[]);

void
cs_f_vof_log_mass_budget(void);

void
cs_f_vof_deshpande_drift_flux(void);

void
cs_f_vof_update_drift_flux(void);

void
cs_f_cavitation_get_pointers(double **presat,
                             double **uinf,
                             double **linf,
                             double **cdest,
                             double **cprod,
                             int    **icvevm,
                             double **mcav,
                             int    **itscvi);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointer to VOF model indicator and parameters
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ivofmt --> pointer to cs_glob_vof_parameters->vof_model
 *   rho1   --> pointer to cs_glob_vof_parameters->rho1
 *   rho2   --> pointer to cs_glob_vof_parameters->rho2
 *   mu1    --> pointer to cs_glob_vof_parameters->mu1
 *   mu2    --> pointer to cs_glob_vof_parameters->mu2
 *   sigmaS --> pointer to cs_glob_vof_parameters->sigmaS
 *   idrift --> pointer to cs_glob_vof_parameters->idrift
 *   cdrift --> pointer to cs_glob_vof_parameters->cdrift
 *   kdrift --> pointer to cs_glob_vof_parameters->kdrift
 *----------------------------------------------------------------------------*/

void
cs_f_vof_get_pointers(unsigned **ivofmt,
                      double   **rho1,
                      double   **rho2,
                      double   **mu1,
                      double   **mu2,
                      double   **sigmaS,
                      int      **idrift,
                      double   **cdrift,
                      double   **kdrift)
{
  *ivofmt = &(_vof_parameters.vof_model);
  *rho1   = &(_vof_parameters.rho1);
  *rho2   = &(_vof_parameters.rho2);
  *mu1    = &(_vof_parameters.mu1);
  *mu2    = &(_vof_parameters.mu2);
  *sigmaS = &(_vof_parameters.sigmaS);
  *idrift = &(_vof_parameters.idrift);
  *cdrift = &(_vof_parameters.cdrift);
  *kdrift = &(_vof_parameters.kdrift);
}

/*----------------------------------------------------------------------------
 * wrapper vof functions, intended for use by Fortran wrapper only.
 *----------------------------------------------------------------------------*/

void
cs_f_vof_compute_linear_rho_mu(void)
{
  cs_vof_compute_linear_rho_mu(cs_glob_mesh);
}

void
cs_f_vof_update_phys_prop(void)
{
  cs_vof_update_phys_prop(cs_glob_mesh);
}

void
cs_f_vof_surface_tension(cs_real_3_t stf[])
{
  cs_vof_surface_tension(cs_glob_mesh, cs_glob_mesh_quantities, stf);
}

void
cs_f_vof_log_mass_budget(void)
{
  cs_vof_log_mass_budget(cs_glob_mesh, cs_glob_mesh_quantities);
}

void
cs_f_vof_deshpande_drift_flux(void)
{
  cs_vof_deshpande_drift_flux(cs_glob_mesh, cs_glob_mesh_quantities);
}

/*----------------------------------------------------------------------------
 * Get pointer to cavitation model indicator and parameters
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   presat --> pointer to cs_glob_cavitation_parameters->presat
 *   uinf   --> pointer to cs_glob_cavitation_parameters->uinf
 *   linf   --> pointer to cs_glob_cavitation_parameters->linf
 *   cdest  --> pointer to cs_glob_cavitation_parameters->cdest
 *   cprod  --> pointer to cs_glob_cavitation_parameters->cprod
 *   icvevm --> pointer to cs_glob_cavitation_parameters->icvevm
 *   mcav   --> pointer to cs_glob_cavitation_parameters->mcav
 *   itscvi --> pointer to cs_glob_cavitation_parameters->itscvi
 *----------------------------------------------------------------------------*/

void
cs_f_cavitation_get_pointers(double **presat,
                             double **uinf,
                             double **linf,
                             double **cdest,
                             double **cprod,
                             int    **icvevm,
                             double **mcav,
                             int    **itscvi)
{
  *presat = &(_cavit_parameters.presat);
  *uinf   = &(_cavit_parameters.uinf);
  *linf   = &(_cavit_parameters.linf);
  *cdest  = &(_cavit_parameters.cdest);
  *cprod  = &(_cavit_parameters.cprod);
  *icvevm = &(_cavit_parameters.icvevm);
  *mcav   = &(_cavit_parameters.mcav);
  *itscvi = &(_cavit_parameters.itscvi);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Provide access to VOF structure.
 */
/*----------------------------------------------------------------------------*/

cs_vof_parameters_t *
cs_get_glob_vof_parameters(void)
{
  return &_vof_parameters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Smoothing a variable after several double-projections
 * cells->faces->cells.
 *
 * \param[in]     m        pointer to mesh structure
 * \param[in]     mq       pointer to mesh quantities structure
 * \param[in]     coefa    boundary condition array for the variable
 * \param[in]     coefb    boundary condition array for the variable
 * \param[in,out] pvar     diffused variable
 */
/*----------------------------------------------------------------------------*/

static void
_smoothe(const cs_mesh_t              *m,
         const cs_mesh_quantities_t   *mq,
         cs_real_t                    *restrict coefa,
         cs_real_t                    *restrict coefb,
         cs_real_t                    *restrict pvar)
{
  int n_cells_ext = m->n_cells_with_ghosts;
  int n_cells     = m->n_cells;
  int n_i_faces   = m->n_i_faces;
  int n_b_faces   = m->n_b_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)m->i_face_cells;

  cs_real_t *restrict dist   = mq->i_dist;
  cs_real_t *restrict pond   = mq->weight;
  cs_real_t *restrict volume = mq->cell_vol;
  cs_real_t *restrict surfn  = mq->i_face_surf;

  cs_real_3_t *restrict surfac = (cs_real_3_t *restrict )mq->i_face_normal;
  cs_real_3_t *restrict diipf  = (cs_real_3_t *restrict )mq->diipf;
  cs_real_3_t *restrict djjpf  = (cs_real_3_t *restrict )mq->djjpf;

  double d_tau = 0.1; /* Sharpening interface on 5 cells (0.1 for 3 cells) */
  /* User Intialization Triple line model */
  /*   alpha_p = 0 - Surface hydrophobe
       alpha_p = 1 - Surface hydrophile
       B_s         - Paramètre de pénalité */
  double B_s = 0.;
  double alpha_p = 0.263544;
  /* User Intialization Triple line model */

  const cs_equation_param_t *eqp_volf
    = cs_field_get_equation_param_const(CS_F_(void_f));

  cs_real_t *viscf, *xam, *dam, *rtpdp, *smbdp;
  BFT_MALLOC(viscf, n_i_faces, cs_real_t);
  BFT_MALLOC(xam, n_i_faces, cs_real_t);
  BFT_MALLOC(dam, n_cells_ext, cs_real_t);
  BFT_MALLOC(rtpdp, n_cells_ext, cs_real_t);
  BFT_MALLOC(smbdp, n_cells_ext, cs_real_t);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    rtpdp[c_id] = 0.;
    smbdp[c_id] = pvar[c_id] * volume[c_id];
  }

  /* PREPARE SYSTEM TO SOLVE */
  /* Compute the gradient of "diffused void fraction" */
  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp_volf->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_gradient_scalar("pvar_grad",
                     gradient_type,
                     halo_type,
                     1,     /* inc */
                     true,  /* iccocg */
                     eqp_volf->nswrgr,
                     0,
                     0,
                     1,     /* w_stride */
                     eqp_volf->verbosity,
                     eqp_volf->imligr,
                     eqp_volf->epsrgr,
                     eqp_volf->climgr,
                     NULL,
                     coefa,
                     coefb,
                     pvar,
                     pond,
                     NULL, /* internal coupling */
                     grad);

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t ii = i_face_cells[f_id][0];
    cs_lnum_t jj = i_face_cells[f_id][1];

    cs_real_t taille = 0.5 * (  pow(volume[ii], cs_math_1ov3)
                              + pow(volume[jj], cs_math_1ov3));
    cs_real_t visco = taille * taille * d_tau;

    cs_real_3_t distxyz;
    for (int i = 0; i < 3; i++)
      distxyz[i] =  dist[f_id] * surfac[f_id][i] / surfn[f_id]
                  + diipf[f_id][i] + djjpf[f_id][i];

    viscf[f_id] = visco * surfn[f_id] / cs_math_3_norm(distxyz);

    cs_real_t *gradi = (cs_real_t *)grad + 3 * ii;
    cs_real_t *gradj = (cs_real_t *)grad + 3 * jj;

    cs_real_t reconstr
      = viscf[f_id] * (  cs_math_3_dot_product(diipf[f_id], gradi)
                       - cs_math_3_dot_product(djjpf[f_id], gradj));

    smbdp[ii] -= reconstr;
    smbdp[jj] += reconstr;
  }

  /* Initialization */
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    dam[c_id] = volume[c_id];

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    xam[f_id] = 0.;

  /* Extra-diagonal terms computation */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    xam[f_id] -= viscf[f_id];

  /* Extra-diagonal terms contribution to the diagonal
     (without boundary contribution (0 flux)) */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t ii = i_face_cells[f_id][0];
    cs_lnum_t jj = i_face_cells[f_id][1];

    dam[ii] -= xam[f_id];
    dam[jj] -= xam[f_id];
  }

  /* Slight re-inforcement of the diagonal if we are not on a
     Dirichlet condition (move the eigenvalues spectrum) */
  cs_real_t epsi = 1.e-7;

  for (int c_id = 0; c_id < n_cells; c_id++)
    dam[c_id] *= (1. + epsi);

  /* SOLVE SYSTEM */
  /* Linear system initialization is in nc_sles_default */
  /* Linear system resolution options */
  cs_real_t precision = 1.e-5; /* Solver precision */
  int n_equiv_iter = 0;     /* Number of equivalent iterative solver iterations */
  cs_real_t residue;        /* Residue */

  /* Linear system resolution */
  /* Get the residue normalization */
  cs_real_t rnorm = cs_dot(n_cells, smbdp, smbdp);

  cs_parall_sum(1, CS_DOUBLE, &rnorm);
  rnorm = sqrt(rnorm); /* Residue normalization */

  /* Triple line model (WARNING: model terms are added after
     residue normalization ?) */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t ii = b_face_cells[f_id];

    cs_real_t is_wall = 0.;
    if (   cs_glob_bc_type[f_id] == CS_SMOOTHWALL
        || cs_glob_bc_type[f_id] == CS_ROUGHWALL)
      is_wall = 1.;

    smbdp[ii] += B_s * volume[ii] * alpha_p * is_wall;
    dam[ii] += volume[ii] * B_s * is_wall;
  }

  cs_sles_solve_native(-1,
                       "ITM_diffusion_equation",
                       true,                   /* symmetric */
                       1, 1,                   /* diag/extradiag block size */
                       dam, xam,
                       precision,
                       rnorm,
                       &n_equiv_iter,
                       &residue,
                       smbdp,
                       rtpdp);

  /* Destroy multigrid structure */
  cs_sles_free_native(-1, "ITM_diffusion_equation");

  /* Increment the distance to the wall */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    pvar[c_id] = rtpdp[c_id];

  BFT_FREE(viscf);
  BFT_FREE(xam);
  BFT_FREE(dam);
  BFT_FREE(rtpdp);
  BFT_FREE(smbdp);
  BFT_FREE(grad);
}

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
cs_vof_compute_linear_rho_mu(const cs_mesh_t  *m)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_t *b_face_cells = m->b_face_cells;

  cs_real_t *cvar_voidf = CS_F_(void_f)->val;
  cs_real_t *a_voidf = CS_F_(void_f)->bc_coeffs->a;
  cs_real_t *b_voidf = CS_F_(void_f)->bc_coeffs->b;

  cs_real_t *cpro_rom = CS_F_(rho)->val;
  cs_real_t *bpro_rom = CS_F_(rho_b)->val;

  cs_real_t *cpro_viscl = CS_F_(mu)->val;

  const cs_real_t rho1 = _vof_parameters.rho1;
  const cs_real_t rho2 = _vof_parameters.rho2;
  const cs_real_t mu1 = _vof_parameters.mu1;
  const cs_real_t mu2 = _vof_parameters.mu2;

  /* Update mixture density and viscosity on cells */

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t vf = cvar_voidf[c_id];
    cpro_rom[c_id]   = rho2*vf + rho1*(1. - vf);
    cpro_viscl[c_id] =  mu2*vf +  mu1*(1. - vf);
  }

  cs_halo_type_t halo_type = m->halo_type;
  cs_field_synchronize(CS_F_(rho), halo_type);
  cs_field_synchronize(CS_F_(mu), halo_type);

  /* Update mixture density on boundary faces */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];
    cs_real_t vf = a_voidf[f_id] + b_voidf[f_id]*cvar_voidf[c_id];

    bpro_rom[f_id] = rho2*vf + rho1*(1. - vf);
  }
}

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
cs_vof_update_phys_prop(const cs_mesh_t  *m)
{
  /* update rho and mu with linear laws */
  cs_vof_compute_linear_rho_mu(m);

  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_real_t rho1 = _vof_parameters.rho1;
  const cs_real_t rho2 = _vof_parameters.rho2;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const int kiflux = cs_field_key_id("inner_flux_id");
  const int kbflux = cs_field_key_id("boundary_flux_id");

  const cs_real_t *restrict i_voidflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(void_f), kiflux))->val;
  const cs_real_t *restrict b_voidflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(void_f), kbflux))->val;

  const cs_real_t *restrict i_volflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(void_f), kimasf))->val;
  const cs_real_t *restrict b_volflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(void_f), kbmasf))->val;

  cs_real_t *restrict i_massflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kimasf))->val;
  cs_real_t *restrict b_massflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf))->val;

  cs_real_t drho = rho2 - rho1;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    i_massflux[f_id] += drho * i_voidflux[f_id] + rho1*i_volflux[f_id];
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    b_massflux[f_id] += drho * b_voidflux[f_id] + rho1*b_volflux[f_id];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write in main log the global mixture mass budget:
 * \f[
 * \sum_\celli\left(
 * |\Omega_\celli|\dfrac{\alpha_\celli^n - \alpha_\celli^{n-1}}{\Delta t} +
 * \sum_{\cellj\in\Face{\celli}}\left(\rho\vect{u}\vect{S}\right)_{\ij}^n
 * \right).
 * \f]
 *
 * \param[in]  m   pointer to mesh structure
 * \param[in]  mq  pointer to mesh quantities structure
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_log_mass_budget(const cs_mesh_t             *m,
                       const cs_mesh_quantities_t  *mq)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_with_ghosts = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const cs_real_t *restrict cell_f_vol = mq->cell_f_vol;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)mq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)mq->b_face_cog;

  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)mq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)mq->b_f_face_normal;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");

  cs_real_t *restrict i_massflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kimasf))->val;
  cs_real_t *restrict b_massflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf))->val;

  cs_real_t *cpro_rom = CS_F_(rho)->val;
  cs_real_t *cproa_rom = CS_F_(rho)->val_pre;
  cs_real_t *bpro_rom = CS_F_(rho_b)->val;

  int icorio = cs_glob_physical_constants->icorio;
  cs_turbomachinery_model_t iturbo = cs_turbomachinery_get_model();

  cs_real_t *i_massflux_abs = NULL, *b_massflux_abs = NULL;

  if (icorio == 1 || iturbo > CS_TURBOMACHINERY_NONE) {

    cs_lnum_t cr_step = 0;
    const int *cell_rotor_num = NULL;

    const int _corio_rotor_num[] = {1};

    if (iturbo > CS_TURBOMACHINERY_NONE) {
      cr_step = 1;
      cell_rotor_num = cs_turbomachinery_get_cell_rotor_num();
    }
    else
      cell_rotor_num = _corio_rotor_num;

    BFT_MALLOC(i_massflux_abs, n_i_faces, cs_real_t);
    BFT_MALLOC(b_massflux_abs, n_b_faces, cs_real_t);

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      i_massflux_abs[f_id] = i_massflux[f_id];

      cs_lnum_t  c_id_i = i_face_cells[f_id][0];
      cs_lnum_t  c_id_j = i_face_cells[f_id][1];
      int rot_ce_i = cell_rotor_num[cr_step*c_id_i];
      int rot_ce_j = cell_rotor_num[cr_step*c_id_j];

      if (rot_ce_i != 0 || rot_ce_j != 0) {
        cs_real_t rhofac = 0.5*(cpro_rom[c_id_i] + cpro_rom[c_id_j]);

        cs_real_t vr1[3], vr2[3];
        cs_rotation_velocity(cs_glob_rotation + rot_ce_i,
                             i_face_cog[f_id],
                             vr1);
        cs_rotation_velocity(cs_glob_rotation + rot_ce_i,
                             i_face_cog[f_id],
                             vr2);
        cs_real_t vr[] = {0.5*(vr1[0]+vr2[0]),
                          0.5*(vr1[1]+vr2[1]),
                          0.5*(vr1[2]+vr2[2])};

        i_massflux_abs[f_id] +=
          rhofac * cs_math_3_dot_product(i_f_face_normal[f_id], vr);
      }
    }

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      b_massflux_abs[f_id] = b_massflux[f_id];

      cs_lnum_t  c_id = b_face_cells[f_id];
      int rot_ce_i = cell_rotor_num[cr_step*c_id];

      if (rot_ce_i != 0) {
        cs_real_t vr[3];
        cs_rotation_velocity(cs_glob_rotation + rot_ce_i, b_face_cog[f_id], vr);

        b_massflux[f_id] +=
          bpro_rom[f_id] * cs_math_3_dot_product(b_f_face_normal[f_id], vr);
      }
    }

    /* massflux point to absolute ones now */
    i_massflux = i_massflux_abs;
    b_massflux = b_massflux_abs;
  }

  /* (Absolute) Mass flux divergence */

  cs_real_t *divro;
  BFT_MALLOC(divro, n_cells_with_ghosts, cs_real_t);
  cs_divergence(m,
                1, /* initialize to 0 */
                i_massflux,
                b_massflux,
                divro);

  if (icorio == 1 || iturbo > CS_TURBOMACHINERY_NONE) {
    BFT_FREE(i_massflux_abs);
    BFT_FREE(b_massflux_abs);
  }
  i_massflux = NULL, b_massflux = NULL;

  /* Unsteady term  and mass budget */

  cs_real_t glob_m_budget = 0.;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t tinsro =  cell_f_vol[c_id]
                      * (cpro_rom[c_id]-cproa_rom[c_id]) / CS_F_(dt)->val[c_id];

    glob_m_budget += tinsro + divro[c_id];
  }

  cs_parall_sum(1, CS_DOUBLE, &glob_m_budget);

  if (cs_log_default_is_active())
    cs_log_printf(CS_LOG_DEFAULT,
                  _("   ** VOF model, mass balance: %12.4e\n\n"),
                  glob_m_budget);

  BFT_FREE(divro);
}

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
                       cs_real_3_t                  stf[])
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)m->i_face_cells;

  cs_real_t *pond = mq->weight;

  cs_real_3_t *restrict surfac = (cs_real_3_t *restrict )mq->i_face_normal;
  cs_real_3_t *restrict dofij = (cs_real_3_t *restrict)mq->dofij;

  const cs_equation_param_t *eqp_volf
    = cs_field_get_equation_param_const(CS_F_(void_f));

  cs_real_t *curv, *coefa, *coefb, *pvar;
  cs_real_3_t *coefa_vec;
  cs_real_33_t *coefb_vec;

  const cs_real_t cpro_surftens = _vof_parameters.sigmaS;

  BFT_MALLOC(curv, n_cells_ext, cs_real_t);
  BFT_MALLOC(pvar, n_cells_ext, cs_real_t);
  BFT_MALLOC(coefa, n_b_faces, cs_real_t);
  BFT_MALLOC(coefb, n_b_faces, cs_real_t);
  BFT_MALLOC(coefa_vec, n_b_faces, cs_real_3_t);
  BFT_MALLOC(coefb_vec, n_b_faces, cs_real_33_t);

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    coefa[face_id] = 0.;
    coefb[face_id] = 1.;

    for (int i = 0; i < 3; i++) {
      coefa_vec[face_id][i] = 0.;

      for (int j = 0; j < 3; j++)
        coefb_vec[face_id][i][j] = 0.;
      coefb_vec[face_id][i][i] = 1.;
    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    pvar[c_id] = 1. - CS_F_(void_f)->val[c_id];
    pvar[c_id] = cs_math_fmin(cs_math_fmax(pvar[c_id], 0.), 1.);
  }

  /* Void fraction diffusion solving */
  int ncycles = 5;
  for (int i = 0; i < ncycles; i++)
    _smoothe(m, mq, coefa, coefb, pvar);

  /* Compute the gradient of "diffused void fraction" */
  cs_real_3_t *surfxyz_unnormed;
  BFT_MALLOC(surfxyz_unnormed, n_cells_ext, cs_real_3_t);

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp_volf->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_gradient_scalar("diff_void_grad",
                     gradient_type,
                     halo_type,
                     1,     /* inc */
                     true,  /* iccocg */
                     eqp_volf->nswrgr,
                     0,
                     0,
                     1,     /* w_stride */
                     eqp_volf->verbosity,
                     eqp_volf->imligr,
                     eqp_volf->epsrgr,
                     eqp_volf->climgr,
                     NULL,
                     coefa,
                     coefb,
                     pvar,
                     pond,
                     NULL, /* internal coupling */
                     surfxyz_unnormed);


  /* Compute the norm of grad(alpha_diffu) */
  cs_real_3_t *surfxyz_norm;
  BFT_MALLOC(surfxyz_norm, n_cells_ext, cs_real_3_t);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    double snorm = cs_math_3_norm(surfxyz_unnormed[c_id]);
    if (snorm < 1.e-10)
      snorm = 1.;

    double unsnorm = 1. / snorm;

    for (int i = 0; i < 3; i++)
      surfxyz_norm[c_id][i] = surfxyz_unnormed[c_id][i] * unsnorm;
  }

  /* Curvature Computation */
  cs_real_33_t *gradnxyz;
  BFT_MALLOC(gradnxyz, n_cells_ext, cs_real_33_t);

  cs_gradient_vector("grad_diff_void_grad",
                     gradient_type,
                     halo_type,
                     1,     /* inc */
                     eqp_volf->nswrgr,
                     eqp_volf->verbosity,
                     eqp_volf->imligr,
                     eqp_volf->epsrgr,
                     eqp_volf->climgr,
                     (const cs_real_3_t *)coefa_vec,
                     (const cs_real_33_t *)coefb_vec,
                     surfxyz_norm,
                     pond,
                     NULL,
                     gradnxyz);

  /* Reconstructions for curvature computation */
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    curv[c_id] = 0.;

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    cs_real_3_t gradf;

    for (int k = 0; k < 3; k++)
      gradf[k] =         pond[face_id]  * surfxyz_norm[ii][k]
                 + (1. - pond[face_id]) * surfxyz_norm[jj][k];

    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        gradf[k] += 0.5 * dofij[face_id][l]
                        * (gradnxyz[ii][k][l] + gradnxyz[jj][k][l]);

    double flux = cs_math_3_dot_product(gradf, surfac[face_id]);
    curv[ii] += flux;
    curv[jj] -= flux;
  }

  /* Compute volumic surface tension */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    for (int i = 0; i < 3; i++) {
      stf[c_id][i] = -cpro_surftens * surfxyz_unnormed[c_id][i] * curv[c_id];
    }

  BFT_FREE(surfxyz_norm);
  BFT_FREE(surfxyz_unnormed);
  BFT_FREE(gradnxyz);
  BFT_FREE(curv);
  BFT_FREE(pvar);
  BFT_FREE(coefa);
  BFT_FREE(coefb);
  BFT_FREE(coefa_vec);
  BFT_FREE(coefb_vec);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a relative velocity \f$ \vect u _d \f$ directly at internal
 *        faces (drift flux), following the approach described by
 *        Suraj S. Deshpande et al 2012 Comput. Sci. Disc. 5 014016.
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
 * {\vect S_{\face'}} \right ] \right ) \left ( \vect n \cdot \vect S \right )
 * _{\face}
 * \f]
 * Where \f$ C_{\gamma} \f$ is the drift flux factor defined with the variable
 * \ref cdrift, \f$ \vect n _{\face} \f$ the normal vector to the interface.
 * The gradient is computed using a centered scheme:
 * \f[
 * \vect {n} _{\face} = \dfrac{\left ( \grad \alpha \right ) _{\face}}
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
                            const cs_mesh_quantities_t  *mq)
{
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_gnum_t n_g_cells = m->n_g_cells;
  const cs_lnum_t n_cells_with_ghosts = m->n_cells_with_ghosts;

  const cs_real_t tot_vol = mq->tot_vol;
  const cs_real_t *i_face_surf = (const cs_real_t *)mq->i_face_surf;
  const cs_real_3_t *i_face_normal = (const cs_real_3_t *)mq->i_face_normal;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)m->i_face_cells;

  /* Constant parameter */
  const cs_real_t cdrift = _vof_parameters.cdrift;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const cs_real_t *restrict i_volflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(void_f), kimasf))->val;

  cs_field_t *idriftflux = NULL;
  idriftflux = cs_field_by_name_try("inner_drift_velocity_flux");

  /* Check if field exists */
  if (idriftflux == NULL)
    bft_error(__FILE__, __LINE__, 0,_("error drift velocity not defined\n"));
  cs_real_t *cpro_idriftf = idriftflux->val;

  cs_real_3_t *voidf_grad;
  BFT_MALLOC(voidf_grad, n_cells_with_ghosts, cs_real_3_t);
  /* Compute the gradient of the void fraction */
  cs_field_gradient_scalar(CS_F_(void_f),
                           true,           // use_previous_t
                           1,              // inc
                           true,           // _recompute_cocg
                           voidf_grad);

  /* Stabilization factor */
  cs_real_t delta = pow(10,-8)/pow(tot_vol/n_g_cells,(1./3.));

  /* Compute the max of flux/Surf over the entire domain*/
  cs_real_t maxfluxsurf = 0.;
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    if (maxfluxsurf < CS_ABS(i_volflux[f_id])/i_face_surf[f_id])
      maxfluxsurf = CS_ABS(i_volflux[f_id])/i_face_surf[f_id];
  }
  cs_parall_max(1, CS_DOUBLE, &maxfluxsurf);

  /* Compute the relative velocity at internal faces */
  cs_real_3_t gradface, normalface;
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t cell_id1 = i_face_cells[f_id][0];
    cs_lnum_t cell_id2 = i_face_cells[f_id][1];
    cs_real_t fluxfactor =
      CS_MIN(cdrift*CS_ABS(i_volflux[f_id])/i_face_surf[f_id], maxfluxsurf);

    for (int idim = 0; idim < 3; idim++)
      gradface[idim] = (  voidf_grad[cell_id1][idim]
                        + voidf_grad[cell_id2][idim])/2.;

    cs_real_t normgrad = sqrt(pow(gradface[0],2)+
                              pow(gradface[1],2)+
                              pow(gradface[2],2));

    for (int idim = 0; idim < 3; idim++)
      normalface[idim] = gradface[idim] / (normgrad+delta);

    cpro_idriftf[f_id] = fluxfactor*(normalface[0]*i_face_normal[f_id][0]+
                                     normalface[1]*i_face_normal[f_id][1]+
                                     normalface[2]*i_face_normal[f_id][2]);
  }

  BFT_FREE(voidf_grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the divergence of the drift velocity term in the volume
 *        fraction equation.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \alpha_\celli^{n+1} \left( 1 - \alpha_\cellj^{n+1} \right) \left(
 *        \dot{m}_\fij^{d} \right)^{+} + \alpha_\cellj^{n+1} \left( 1 -
 *        \alpha_\celli^{n+1} \right) \left( \dot{m}_\fij^{d} \right)^{-}
 *       \right)
 * \f]
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
                  cs_real_t        *restrict rhs)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != NULL)
    cs_sync_scalar_halo(m, pvar);
  else if (pvara == NULL)
    pvara = (const cs_real_t *restrict)pvar;

  const cs_real_t  *restrict _pvar = (pvar != NULL) ? pvar : pvara;

  /*======================================================================
    Computation of the drift flux
    ======================================================================*/

  cs_field_t *vr = cs_field_by_name_try("drift_velocity");
  cs_field_t *idriftflux = cs_field_by_name_try("inner_drift_velocity_flux");
  cs_field_t *bdriftflux = cs_field_by_name_try("boundary_drift_velocity_flux");

  if (_vof_parameters.idrift == 1) {

    // FIXME Handle boundary terms bdriftflux
    cs_vof_deshpande_drift_flux(cs_glob_mesh, cs_glob_mesh_quantities);

  } else {

    const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
    int f_id, itypfl, iflmb0, init, inc;

    cs_real_3_t *coefav;
    cs_real_33_t *coefbv;

    /* Check if field exist */
    if (idriftflux == NULL)
      bft_error(__FILE__, __LINE__, 0,_("error drift velocity not defined\n"));

    cs_real_3_t *cpro_vr = (cs_real_3_t *)vr->val;
    cs_real_t *cpro_idriftf = idriftflux->val;
    cs_real_t *cpro_bdriftf = bdriftflux->val;

    BFT_MALLOC(coefav, n_b_faces, cs_real_3_t);
    BFT_MALLOC(coefbv, n_b_faces, cs_real_33_t);

    f_id = -1;
    itypfl = 0;
    iflmb0 = 1;
    init = 1;
    inc = 1;

    /* Boundary coefficients */
    for (cs_lnum_t ifac = 0 ; ifac < n_b_faces ; ifac++) {
      for (int ii = 0 ; ii < 3 ; ii++) {
        coefav[ifac][ii] = 0.;
        for (int jj = 0 ; jj < 3 ; jj++) {
          coefbv[ifac][ii][jj] = 0.;
        }
        coefbv[ifac][ii][ii] = 1.;
      }
    }

    cs_mass_flux(m,
                 fvq,
                 f_id,
                 itypfl,
                 iflmb0,
                 init,
                 inc,
                 imrgra,
                 nswrgp,
                 imligp,
                 iwarnp,
                 epsrgp,
                 climgp,
                 NULL, /* rom */
                 NULL, /* romb */
                 (const cs_real_3_t *)cpro_vr,
                 (const cs_real_3_t *)coefav,
                 (const cs_real_33_t *)coefbv,
                 cpro_idriftf,
                 cpro_bdriftf);

    BFT_FREE(coefav);
    BFT_FREE(coefbv);

  }

  /*======================================================================
    Contribution from interior faces
    ======================================================================*/

  const int kiflux = cs_field_key_id("inner_flux_id");
  int i_flux_id = cs_field_get_key_int(CS_F_(void_f), kiflux);
  cs_field_t *i_flux = cs_field_by_id(i_flux_id);

  if (n_cells_ext>n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      rhs[cell_id] = 0.;
    }
  }

  for (int g_id = 0; g_id < n_i_groups; g_id++) {
#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {
      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        cs_real_t irvf = 0.;
        if (idriftflux != NULL)
          irvf = idriftflux->val[face_id];

        cs_real_2_t fluxij = {0.,0.};

        cs_i_conv_flux(1,
                       1.,
                       0,
                       _pvar[ii],
                       _pvar[jj],
                       _pvar[ii]*(1.-_pvar[jj]),
                       _pvar[ii]*(1.-_pvar[jj]),
                       _pvar[jj]*(1.-_pvar[ii]),
                       _pvar[jj]*(1.-_pvar[ii]),
                       irvf,
                       1.,
                       1.,
                       fluxij);

        const cs_real_t kdrift = _vof_parameters.kdrift;
        cs_i_diff_flux(1,
                       1.,
                       _pvar[ii],
                       _pvar[jj],
                       _pvar[ii],
                       _pvar[jj],
                       kdrift*(2.-_pvar[ii]-_pvar[jj])
                        / 2.*i_face_surf[face_id]/i_dist[face_id],
                       fluxij);

        rhs[ii] -= fluxij[0];
        rhs[jj] += fluxij[1];
        /* store void fraction convection flux contribution */
        i_flux->val[face_id] += fluxij[0];
      }
    }
  }
}

/*----------------------------------------------------------------------------
 *!
 * \brief Provide access to cavitation parameters structure.
 */
/*----------------------------------------------------------------------------*/

cs_cavitation_parameters_t *
cs_get_glob_cavitation_parameters(void)
{
  return &_cavit_parameters;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
