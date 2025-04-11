/*============================================================================
 * Functions associated to VOF model
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"

#include "base/cs_array.h"
#include "base/cs_base.h"
#include "alge/cs_blas.h"
#include "base/cs_boundary_conditions.h"
#include "alge/cs_convection_diffusion.h"
#include "alge/cs_convection_diffusion_priv.h"
#include "cs_dispatch.h"
#include "alge/cs_divergence.h"
#include "cdo/cs_equation.h"
#include "base/cs_equation_iterative_solve.h"
#include "alge/cs_face_viscosity.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_log.h"
#include "base/cs_log_iteration.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_location.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_physical_constants.h"
#include "base/cs_prototypes.h"
#include "base/cs_reducers.h"
#include "base/cs_rotation.h"
#include "alge/cs_sles_default.h"
#include "base/cs_turbomachinery.h"
#include "turb/cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_vof.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_vof.cpp
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

  \var  cs_vof_parameters_t::sigma_s
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
  .sigma_s        = 0.,
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

/* Contact angle choice */
static cs_vof_contact_angle_t _contact_angle_choice = CS_VOF_CONTACT_ANGLE_OFF;

/*============================================================================
 * Global variables
 *============================================================================*/

const cs_vof_parameters_t *cs_glob_vof_parameters = &_vof_parameters;
const cs_cavitation_parameters_t *cs_glob_cavitation_parameters
  = &_cavit_parameters;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_vof_compute_linear_rho_mu(void);

void
cs_f_vof_deshpande_drift_flux(void);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * wrapper vof functions, intended for use by Fortran wrapper only.
 *----------------------------------------------------------------------------*/

void
cs_f_vof_compute_linear_rho_mu(void)
{
  cs_vof_compute_linear_rho_mu(cs_glob_mesh);
}

void
cs_f_vof_deshpande_drift_flux(void)
{
  cs_vof_deshpande_drift_flux(cs_glob_mesh, cs_glob_mesh_quantities);
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Smoothing a variable after several double-projections
 * cells->faces->cells.
 *
 * \param[in]     m          pointer to mesh structure
 * \param[in]     mq         pointer to mesh quantities structure
 * \param[in]     bc_coeffs  boundary condition structure for the variable
 * \param[in,out] pvar       diffused variable
 */
/*----------------------------------------------------------------------------*/

static void
_smoothe(const cs_mesh_t              *m,
         const cs_mesh_quantities_t   *mq,
         const cs_field_bc_coeffs_t   *bc_coeffs,
         cs_real_t                    *restrict pvar)
{
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_2_t *i_face_cells = m->i_face_cells;

  const cs_real_t *restrict dist = mq->i_dist;
  const cs_real_t *restrict cell_vol = mq->cell_vol;
  const cs_real_t *restrict i_face_surf = mq->i_face_surf;

  const cs_nreal_3_t *restrict i_face_u_normal = mq->i_face_u_normal;
  const cs_rreal_3_t *restrict diipf = mq->diipf;
  const cs_rreal_3_t *restrict djjpf = mq->djjpf;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;
  const cs_lnum_t *cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *cell_b_faces = ma->cell_b_faces;

  const int *bc_type = cs_glob_bc_type;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);

  const cs_real_t d_tau = 0.1; /* Sharpening interface on 5 cells
                                  (0.1 for 3 cells) */
  /* User Intialization Triple line model */
  /*   alpha_p = 0 - hydrophobic surface
       alpha_p = 1 - hydrophilic surface
       b_s         - penalty parameter */
  const cs_real_t b_s = 0.;
  const cs_real_t alpha_p = 0.263544;
  /* User Intialization Triple line model */

  const cs_equation_param_t *eqp_volf
    = cs_field_get_equation_param_const(CS_F_(void_f));

  cs_real_t *xam, *dam, *smbdp;
  CS_MALLOC_HD(xam, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(dam, n_cells, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(smbdp, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Prepare system to solve */
  /* Compute the gradient of "diffused void fraction" */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp_volf->imrgra,
                             &gradient_type,
                             &halo_type);
  cs_real_t *gweight = nullptr;
  int w_stride = 1;
  if (eqp_volf->iwgrec == 1 && eqp_volf->idiff > 0) {
    int key_id = cs_field_key_id("gradient_weighting_id");
    int diff_id = cs_field_get_key_int(CS_F_(void_f), key_id);
    if (diff_id > -1) {
      cs_field_t *weight_f = cs_field_by_id(diff_id);
      gweight = weight_f->val;
      w_stride = weight_f->dim;
    }
  }

  cs_real_3_t *grad;
  CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  cs_gradient_scalar("pvar_grad",
                     gradient_type,
                     halo_type,
                     1,     /* inc */
                     eqp_volf->nswrgr,
                     0,
                     w_stride,
                     eqp_volf->verbosity,
                     (cs_gradient_limit_t)eqp_volf->imligr,
                     eqp_volf->epsrgr,
                     eqp_volf->climgr,
                     nullptr,
                     bc_coeffs,
                     pvar,
                     gweight,
                     nullptr, /* internal coupling */
                     grad);

  constexpr cs_real_t c_1ov3 = 1. / 3.;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    smbdp[c_id] = 0.;
  });

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
    cs_lnum_t ii = i_face_cells[f_id][0];
    cs_lnum_t jj = i_face_cells[f_id][1];

    cs_real_t taille = 0.5 * (  pow(cell_vol[ii], c_1ov3)
                              + pow(cell_vol[jj], c_1ov3));
    cs_real_t visco = taille * taille * d_tau;

    cs_real_t distxyz[3];
    for (cs_lnum_t i = 0; i < 3; i++) {
      distxyz[i] =   dist[f_id] * i_face_u_normal[f_id][i]
                   + diipf[f_id][i] + djjpf[f_id][i];
    }

    cs_real_t visc_f = visco * i_face_surf[f_id] / cs_math_3_norm(distxyz);

    /* Extra-diagonal terms computation */
    xam[f_id] = -visc_f;

    cs_real_t reconstr
      = visc_f * (  cs_math_3_dot_product(diipf[f_id], grad[ii])
                  - cs_math_3_dot_product(djjpf[f_id], grad[jj]));

    if (ii < n_cells)
      cs_dispatch_sum(&smbdp[ii], -reconstr, i_sum_type);
    if (jj < n_cells)
      cs_dispatch_sum(&smbdp[jj], reconstr, i_sum_type);
  });

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* Initialization */
    smbdp[c_id] += pvar[c_id] * cell_vol[c_id];

    cs_real_t dam_c = cell_vol[c_id];

    /* Loop on interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

    /* Extra-diagonal terms contribution to the diagonal
       (without boundary contribution (0 flux)) */

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t f_id = cell_i_faces[cidx];
      dam_c -= xam[f_id];
    }

    /* Slight re-inforcement of the diagonal if we are not on a
       Dirichlet condition (move the eigenvalues spectrum) */
    constexpr cs_real_t epsi = 1.e-7;
    dam_c *= (1. + epsi);

    /* Triple line model (WARNING: model terms are added after
       residual normalization ?) */

    /* Loop on boundary faces */

    if (b_s > 0.) {
      const cs_lnum_t s_id_b = cell_b_faces_idx[c_id];
      const cs_lnum_t e_id_b = cell_b_faces_idx[c_id + 1];

      for (cs_lnum_t cidx = s_id_b; cidx < e_id_b; cidx++) {
        const cs_lnum_t f_id = cell_b_faces[cidx];

        if (   bc_type[f_id] == CS_SMOOTHWALL
            || bc_type[f_id] == CS_ROUGHWALL) {
          smbdp[c_id] += b_s * cell_vol[c_id] * alpha_p;
          dam_c += cell_vol[c_id] * b_s;
        }
      }
    }

    dam[c_id] = dam_c;
    pvar[c_id] = 0.; // Set to 0 before solving.
  });

  ctx.wait();

  CS_FREE_HD(grad);

  /* FIXME: the following synchronization should not
     be needed, unless perhaps when using a native matrix format */
  cs_halo_sync_var(m->halo, CS_HALO_STANDARD, smbdp);

  /* SOLVE SYSTEM
     Linear system initialization is in nc_sles_default
     Linear system resolution options */

  cs_real_t precision = 1.e-5; /* Solver precision */
  int n_equiv_iter = 0;        /* Number of equivalent iterative
                                  solver iterations */
  cs_real_t residual;          /* Residual */

  /* Linear system resolution */
  /* Get the residual normalization */
  cs_real_t rnorm = cs_dot(n_cells, smbdp, smbdp);

  cs_parall_sum(1, CS_REAL_TYPE, &rnorm);
  rnorm = sqrt(rnorm); /* Residual normalization */

  cs_sles_solve_native(-1,
                       "ITM_diffusion_equation",
                       true,                   /* symmetric */
                       1, 1,                   /* diag/extradiag block size */
                       dam, xam,
                       precision,
                       rnorm,
                       &n_equiv_iter,
                       &residual,
                       smbdp,
                       pvar);                  /* Distance to the wall */

  /* Free solver setup and associated matrix. */
  cs_sles_free_native(-1, "ITM_diffusion_equation");

  ctx.wait();

  cs_halo_sync_var(m->halo, CS_HALO_STANDARD, pvar);

  CS_FREE_HD(xam);
  CS_FREE_HD(dam);
  CS_FREE_HD(smbdp);
}

END_C_DECLS // Temporary used here for the templated function

/*----------------------------------------------------------------------------*/
/*!
 * \brief Contact angle based correction for the boundary condition.
 */
/*----------------------------------------------------------------------------*/

template <cs_vof_contact_angle_t choice>
static void
_contact_angle_correction
(
  const cs_mesh_t            *m,         /*!<[in] pointer to mesh structure */
  const cs_mesh_quantities_t *mq,        /*!<[in] pointer to mesh quantities */
  cs_dispatch_context        &ctx,       /*!<[in] reference to dispatch context */
  cs_real_3_t                *surf_norm  /*!<[in,out] Corrected norm of grad(alpha_diffu) */
)
{
  /* Sanity check */
  if (choice < CS_VOF_CONTACT_ANGLE_STATIC)
    return;

  /* Get values or local pointers */
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  cs_real_3_t *b_face_u_normal = mq->b_face_u_normal;
  cs_real_t *restrict volume = mq->cell_vol;

  const cs_real_t cpro_surftens = _vof_parameters.sigma_s;
  const cs_real_t mu1 = _vof_parameters.mu1;

  const cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;

  cs_real_t *alpha_g = CS_F_(void_f)->val;

  // D. Legendre, Computers & Fluids 113 (2015) 2–13 (for these two values)
  constexpr cs_real_t slip_length = 1.e-9;
  constexpr cs_real_t angle = 60.;

  const cs_real_t pi_inv = 180. / cs_math_pi;
  const cs_real_t theta_micro = cs_math_pi * angle / 180.;

  /* Static contact angle only depends on theta_micro, hence we compute it
   * outside of the loop on boundary faces.
   */
  const cs_real_t static_contact_angle
    = (cs_math_pow3(theta_micro) / 9.)
    - 0.00183985 * pow(theta_micro, 4.5)
    + 1.845823e-06 * pow(theta_micro, 12.258487);

  const int *bc_type = cs_glob_bc_type;

  /* Dispatch class */
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);


  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

    if (   bc_type[face_id] == CS_SMOOTHWALL
        || bc_type[face_id] == CS_ROUGHWALL)
    {
      cs_lnum_t cell_id = b_face_cells[face_id];

      if ((1. - alpha_g[cell_id]) * alpha_g[cell_id] > 0.01) {

        // Tangential velocity (tangential projection to the wall)
        cs_real_3_t nt;
        cs_math_3_orthogonal_projection(b_face_u_normal[face_id],
                                        vel[cell_id],
                                        nt);

        /* Compute theta depending on the chosen model. */
        cs_real_t theta = theta_micro;

        if (choice == CS_VOF_CONTACT_ANGLE_DYN) {
          cs_real_t utau = cs_math_3_norm(nt);

          cs_real_t capillary_number = mu1 * utau / cpro_surftens;

          cs_real_t apparent_length = cbrt(volume[cell_id]) / 2.0;

          cs_real_t theta_macro =  static_contact_angle
                                 + capillary_number
                                 * log(apparent_length / slip_length);

          theta = (  cbrt(9.0 * theta_macro)
                   + 0.0727387 * theta_macro
                   - 0.0515388 * cs_math_pow2(theta_macro)
                   + 0.00341336 * cs_math_pow3(theta_macro)) * pi_inv;
        }

        cs_math_3_normalize(nt, nt);

        cs_real_t ctheta = cos(theta);
        cs_real_t stheta = sin(theta);
        for (int i = 0; i < 3; i++)
          cs_dispatch_sum(&surf_norm[cell_id][i],
                          - b_face_u_normal[face_id][i] * ctheta + nt[i] * stheta,
                          b_sum_type);
      }
    }
  });

  ctx.wait();
}

BEGIN_C_DECLS

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Provide write access to VOF structure.
 */
/*----------------------------------------------------------------------------*/

cs_vof_parameters_t *
cs_get_glob_vof_parameters(void)
{
  return &_vof_parameters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create VoF fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_field_create(void)
{
  if (cs_glob_vof_parameters->vof_model < 1)
    return;

  // Key ids for logging and visualization
  const int keylog = cs_field_key_id("log");
  const int keyvis = cs_field_key_id("post_vis");

  // Key id for mass flux
  const int k_imasf = cs_field_key_id("inner_mass_flux_id");
  const int k_bmasf = cs_field_key_id("boundary_mass_flux_id");

  // Key id for flux id
  const int kiflux = cs_field_key_id("inner_flux_id");
  const int kbflux = cs_field_key_id("boundary_flux_id");

  /* Interior faces*/

  cs_field_t *f_ivf = cs_field_create("inner_volume_flux",
                                      CS_FIELD_EXTENSIVE | CS_FIELD_PROPERTY,
                                      CS_MESH_LOCATION_INTERIOR_FACES,
                                      1,
                                      false);
  cs_field_set_key_int(CS_F_(void_f), k_imasf, f_ivf->id);
  cs_field_set_key_int(f_ivf, keyvis, 0);
  cs_field_set_key_int(f_ivf, keylog, 0);

  cs_field_t *f_ivff = cs_field_create("inner_void_fraction_flux",
                                       CS_FIELD_EXTENSIVE | CS_FIELD_PROPERTY,
                                       CS_MESH_LOCATION_INTERIOR_FACES,
                                       1,
                                       false);
  cs_field_set_key_int(CS_F_(void_f), kiflux, f_ivff->id);
  cs_field_set_key_int(f_ivff, keyvis, 0);
  cs_field_set_key_int(f_ivff, keylog, 0);

  /* Boundary faces*/

  cs_field_t *f_bvf = cs_field_create("boundary_volume_flux",
                                      CS_FIELD_EXTENSIVE | CS_FIELD_PROPERTY,
                                      CS_MESH_LOCATION_BOUNDARY_FACES,
                                      1,
                                      false);
  cs_field_set_key_int(CS_F_(void_f), k_bmasf, f_bvf->id);
  cs_field_set_key_int(f_bvf, keyvis, 0);
  cs_field_set_key_int(f_bvf, keylog, 0);

  cs_field_t *f_bvff = cs_field_create("boundary_void_fraction_flux",
                                       CS_FIELD_EXTENSIVE | CS_FIELD_PROPERTY,
                                       CS_MESH_LOCATION_BOUNDARY_FACES,
                                       1,
                                       false);
  cs_field_set_key_int(CS_F_(void_f), kbflux, f_bvff->id);
  cs_field_set_key_int(f_bvff, keyvis, 0);
  cs_field_set_key_int(f_bvff, keylog, 0);

  /* Drift */

  if (cs_glob_vof_parameters->idrift > 0) {

    cs_field_t *f_idvf
      = cs_field_create("inner_drift_velocity_flux",
                        CS_FIELD_EXTENSIVE | CS_FIELD_PROPERTY,
                        CS_MESH_LOCATION_INTERIOR_FACES,
                        1,
                        false);
    cs_field_set_key_int(f_idvf, keyvis, 0);
    cs_field_set_key_int(f_idvf, keylog, 0);

    cs_field_t *f_bdvf
      = cs_field_create("boundary_drift_velocity_flux",
                        CS_FIELD_EXTENSIVE | CS_FIELD_PROPERTY,
                        CS_MESH_LOCATION_BOUNDARY_FACES,
                        1,
                        false);
    cs_field_set_key_int(f_bdvf, keyvis, 0);
    cs_field_set_key_int(f_bdvf, keylog, 0);

  }

  /* Cavitation model

     Liquid-vapor mass transfer term for cavitating flows
     and its part implicit in pressure. */

  if (cs_glob_vof_parameters->vof_model & CS_VOF_MERKLE_MASS_TRANSFER) {

    cs_field_t *f_gamma = cs_field_create("model:cavitation_gamma",
                                          CS_FIELD_INTENSIVE,
                                          CS_MESH_LOCATION_CELLS,
                                          1,
                                          false);
    cs_field_set_key_int(f_gamma, keyvis, 0);
    cs_field_set_key_int(f_gamma, keylog, 0);

    cs_field_t *f_dgdp = cs_field_create("model:cavitation_st_dgdp",
                                          CS_FIELD_INTENSIVE,
                                          CS_MESH_LOCATION_CELLS,
                                          1,
                                          false);
    cs_field_set_key_int(f_dgdp, keyvis, 0);
    cs_field_set_key_int(f_dgdp, keylog, 0);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log setup of VoF model.
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_log_setup(void)
{
  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Homogeneous mixture model VoF\n"
                  "-----------------------------\n\n"));

  auto vp = cs_glob_vof_parameters;
  if (vp->vof_model <= 0) {
    cs_log_printf(CS_LOG_SETUP,
                  _("  Disabled (%d)\n"), vp->vof_model);
    return;
  }

  cs_log_printf(CS_LOG_SETUP,
                _("  Enabled (%d)\n"), vp->vof_model);

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n"
       "  Fluid 1:\n"
       "    ro1:         %14.5e (Reference density)\n"
       "    mu1:         %14.5e (Ref. molecular dyn. visc.)\n"
       "  Fluid 2:\n"
       "    rho2:        %14.5e (Reference density)\n"
       "    mu2:         %14.5e (Ref. molecular dyn. visc.)\n"
       "  Surface tension:\n"
       "    sigma:       %14.5e (surface tension)\n"
       "  Drift velocity:\n"
       "    idrift:          %10d (0: disabled; > 0: enabled)\n"
       "    kdrift       %14.5e (Diffusion effect coeff.)\n"
       "    cdrift       %14.5e (Diffusion flux coeff.)\n"),
     vp->rho1, vp->mu1, vp->rho2, vp->mu2, vp->sigma_s,
     vp->idrift, vp->cdrift, vp->kdrift);

  /* cavitation model */

  if (vp->vof_model & CS_VOF_MERKLE_MASS_TRANSFER) {

    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    "Cavitation model\n"
                    "----------------\n\n"));

    const cs_cavitation_parameters_t *cvp
      = cs_get_glob_cavitation_parameters();

    cs_log_printf
      (CS_LOG_SETUP,
     _("\n"
       "  Liquid phase: fluid 1:\n"
       "  Gas phase: fluid 2:\n"
       "  Vaporization/condensation model (Merkle):\n"
       "    presat:      %14.5e (saturation pressure)\n"
       "    linf:        %14.5e (reference lenght scale)\n"
       "    uinf:        %14.5e (reference velocity)\n"),
       cvp->presat, cvp->linf, cvp->uinf);

    int itytur = cs_glob_turb_model->itytur;

    if (itytur == 2 || itytur == 5 || itytur == 6 || itytur == 7)
      cs_log_printf
        (CS_LOG_SETUP,
         _("\n"
           "  Eddy-viscosity correction (Rebound correction):\n"
           "    icvevm:          %10d (activated (1) or not (0)\n"
           "    mcav:        %14.5e (mcav constant)\n"),
         cvp->icvevm, cvp->mcav);
  }
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

  cs_dispatch_context ctx;

  /* Update mixture density and viscosity on cells */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  c_id) {
    cs_real_t vf = cvar_voidf[c_id];
    cpro_rom[c_id]   = rho2*vf + rho1*(1. - vf);
    cpro_viscl[c_id] =  mu2*vf +  mu1*(1. - vf);
  });

  /* Update mixture density on boundary faces */

  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  f_id) {
    cs_lnum_t c_id = b_face_cells[f_id];
    cs_real_t vf = a_voidf[f_id] + b_voidf[f_id]*cvar_voidf[c_id];

    bpro_rom[f_id] = rho2*vf + rho1*(1. - vf);
  });

  ctx.wait();

  cs_halo_type_t halo_type = m->halo_type;
  cs_field_synchronize(CS_F_(rho), halo_type);
  cs_field_synchronize(CS_F_(mu), halo_type);
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

  cs_dispatch_context ctx;

  ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    i_massflux[f_id] += drho * i_voidflux[f_id] + rho1 * i_volflux[f_id];
  });

  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    b_massflux[f_id] += drho * b_voidflux[f_id] + rho1 * b_volflux[f_id];
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write in main log the global mixture mass budget:
 * \f[
 * \sum_\celli\left(
 * |\Omega_\celli|\dfrac{\rho_\celli^n - \rho_\celli^{n-1}}{\Delta t} +
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

  const cs_lnum_2_t *i_face_cells = m->i_face_cells;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const cs_real_t *restrict cell_f_vol = mq->cell_vol;
  const cs_real_3_t *restrict i_face_cog = mq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = mq->b_face_cog;

  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *)mq->i_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *)mq->b_face_normal;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");

  cs_real_t *restrict i_massflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kimasf))->val;
  cs_real_t *restrict b_massflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf))->val;

  cs_real_t *cpro_rom = CS_F_(rho)->val;
  cs_real_t *cproa_rom = CS_F_(rho)->val_pre;
  cs_real_t *bpro_rom = CS_F_(rho_b)->val;
  const cs_real_t *dt = CS_F_(dt)->val;

  int icorio = cs_glob_physical_constants->icorio;
  cs_turbomachinery_model_t iturbo = cs_turbomachinery_get_model();

  cs_real_t *i_massflux_abs = nullptr, *b_massflux_abs = nullptr;

  cs_dispatch_context ctx;

  if (icorio == 1 || iturbo > CS_TURBOMACHINERY_NONE) {

    cs_lnum_t cr_step = 0;
    const int *cell_rotor_num = nullptr;

    const int _corio_rotor_num[] = {1};

    const cs_rotation_t *r = cs_glob_rotation;

    if (iturbo > CS_TURBOMACHINERY_NONE) {
      cr_step = 1;
      cell_rotor_num = cs_turbomachinery_get_cell_rotor_num();
    }
    else
      cell_rotor_num = _corio_rotor_num;

    CS_MALLOC_HD(i_massflux_abs, n_i_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(b_massflux_abs, n_b_faces, cs_real_t, cs_alloc_mode);

    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      i_massflux_abs[f_id] = i_massflux[f_id];

      cs_lnum_t c_id_i = i_face_cells[f_id][0];
      cs_lnum_t c_id_j = i_face_cells[f_id][1];
      const cs_lnum_t id_i = cr_step*c_id_i;
      const cs_lnum_t id_j = cr_step*c_id_j;
      cs_lnum_t rot_ce_i = cell_rotor_num[id_i];
      cs_lnum_t rot_ce_j = cell_rotor_num[id_j];

      if (rot_ce_i != 0 || rot_ce_j != 0) {
        cs_real_t rhofac = 0.5 * (cpro_rom[c_id_i] + cpro_rom[c_id_j]);

        cs_real_t vr1[3], vr2[3];
        /* cs_rotation_velocity is actually a host function.
           the code will crash on GPU */
        cs_rotation_velocity(r + rot_ce_i,
                             i_face_cog[f_id],
                             vr1);

        cs_rotation_velocity(r + rot_ce_i,
                             i_face_cog[f_id],
                             vr2);

        cs_real_t vr[3] = {0.5 * (vr1[0] + vr2[0]),
                           0.5 * (vr1[1] + vr2[1]),
                           0.5 * (vr1[2] + vr2[2])};

        i_massflux_abs[f_id]
          += rhofac * cs_math_3_dot_product(i_f_face_normal[f_id],vr);
      }
    });

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      b_massflux_abs[f_id] = b_massflux[f_id];

      cs_lnum_t c_id = b_face_cells[f_id];
      const cs_lnum_t id = cr_step*c_id;
      cs_lnum_t rot_ce_i = cell_rotor_num[id];

      if (rot_ce_i != 0) {
        cs_real_t vr[3];
        cs_rotation_velocity(r + rot_ce_i, b_face_cog[f_id], vr);

        b_massflux_abs[f_id]
          += bpro_rom[f_id] * cs_math_3_dot_product(b_f_face_normal[f_id], vr);

      }
    });

    ctx.wait();

    /* massflux point to absolute ones now */
    i_massflux = i_massflux_abs;
    b_massflux = b_massflux_abs;
  }

  if (icorio == 1 || iturbo > CS_TURBOMACHINERY_NONE) {
    CS_FREE_HD(i_massflux_abs);
    CS_FREE_HD(b_massflux_abs);
  }

  /* Unsteady term and mass budget */

  if (cs_log_default_is_active()) {

    /* (Absolute) Mass flux divergence */

    cs_real_t *divro;
    CS_MALLOC_HD(divro, n_cells_with_ghosts, cs_real_t, cs_alloc_mode);
    cs_divergence(m,
                  1, /* initialize to 0 */
                  i_massflux,
                  b_massflux,
                  divro);

    double glob_m_budget = 0.;

    ctx.parallel_for_reduce_sum
      (n_cells, glob_m_budget, [=] CS_F_HOST_DEVICE
       (cs_lnum_t c_id,
        CS_DISPATCH_REDUCER_TYPE(double) &sum) {

      cs_real_t tinsro = cell_f_vol[c_id] * (cpro_rom[c_id]-cproa_rom[c_id]) / dt[c_id];

      sum += (tinsro + divro[c_id]);
    });

    ctx.wait();

    cs_parall_sum(1, CS_DOUBLE, &glob_m_budget);

    cs_log_printf(CS_LOG_DEFAULT,
                  _("   ** VOF model, mass balance: %12.4e\n\n"),
                  glob_m_budget);

    CS_FREE_HD(divro);
  }

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

  cs_mesh_adjacencies_update_cell_i_faces();

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *c2c = ma->cell_cells;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const short int *c2f_sgn = ma->cell_i_faces_sgn;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;

  const cs_real_t *restrict pond = mq->weight;
  const cs_real_3_t *restrict i_face_u_normal
    = (const cs_real_3_t *)mq->i_face_u_normal;
  const cs_real_3_t *restrict dofij = (const cs_real_3_t *)mq->dofij;
  const cs_real_t *restrict i_face_surf = (const cs_real_t *)mq->i_face_surf;

  const cs_equation_param_t *eqp_volf
    = cs_field_get_equation_param_const(CS_F_(void_f));

  const cs_real_t cpro_surftens = _vof_parameters.sigma_s;

  cs_real_t *pvar;
  CS_MALLOC_HD(pvar, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Boundary condition */

  cs_field_bc_coeffs_t bc_coeffs_loc;
  cs_field_bc_coeffs_init(&bc_coeffs_loc);
  CS_MALLOC_HD(bc_coeffs_loc.a, n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(bc_coeffs_loc.b, n_b_faces, cs_real_t, cs_alloc_mode);
  cs_real_t *coefa = bc_coeffs_loc.a;
  cs_real_t *coefb = bc_coeffs_loc.b;

  cs_field_bc_coeffs_t bc_coeffs_v_loc;
  cs_field_bc_coeffs_init(&bc_coeffs_v_loc);
  CS_MALLOC_HD(bc_coeffs_v_loc.a, 3*n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(bc_coeffs_v_loc.b, 9*n_b_faces, cs_real_t, cs_alloc_mode);
  cs_real_3_t  *coefa_vec = (cs_real_3_t  *)bc_coeffs_v_loc.a;
  cs_real_33_t *coefb_vec = (cs_real_33_t *)bc_coeffs_v_loc.b;

  cs_dispatch_context ctx;

  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
    coefa[face_id] = 0.;
    coefb[face_id] = 1.;

    for (cs_lnum_t i = 0; i < 3; i++) {
      coefa_vec[face_id][i] = 0.;

      for (cs_lnum_t j = 0; j < 3; j++)
        coefb_vec[face_id][i][j] = 0.;
      coefb_vec[face_id][i][i] = 1.;
    }
  });

  cs_real_t *cvar_voidf = CS_F_(void_f)->val;
  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  c_id) {
    pvar[c_id] = cs::min(cs::max(1.-cvar_voidf[c_id], 0.), 1.);
  });

  ctx.wait();

  cs_halo_sync_var(m->halo, CS_HALO_STANDARD, pvar);

  /* Void fraction diffusion solving */
  int ncycles = 5;
  for (cs_lnum_t i = 0; i < ncycles; i++) {
    _smoothe(m, mq, &bc_coeffs_loc, pvar);
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  c_id) {
      pvar[c_id] = (pvar[c_id] <= 0.001) ? 0. : pvar[c_id];
      pvar[c_id] = (pvar[c_id] >= 0.999) ? 1. : pvar[c_id];
    });

    ctx.wait();
  }

  /* Compute the gradient of "diffused void fraction" */
  cs_real_3_t *surfxyz_unnormed;
  CS_MALLOC_HD(surfxyz_unnormed, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp_volf->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_real_t *gweight = nullptr;
  int w_stride = 1;
  if (eqp_volf->iwgrec == 1 && eqp_volf->idiff > 0) {
    int key_id = cs_field_key_id("gradient_weighting_id");
    int diff_id = cs_field_get_key_int(CS_F_(void_f), key_id);
    if (diff_id > -1) {
      cs_field_t *weight_f = cs_field_by_id(diff_id);
      gweight = weight_f->val;
      w_stride = weight_f->dim;
    }
  }

  cs_gradient_scalar("diff_void_grad",
                     gradient_type,
                     halo_type,
                     1,     /* inc */
                     eqp_volf->nswrgr,
                     0,
                     w_stride,
                     eqp_volf->verbosity,
                     (cs_gradient_limit_t)eqp_volf->imligr,
                     eqp_volf->epsrgr,
                     eqp_volf->climgr,
                     nullptr,
                     &bc_coeffs_loc,
                     pvar,
                     gweight,
                     nullptr, /* internal coupling */
                     surfxyz_unnormed);

  /* Compute the norm of grad(alpha_diffu) */
  cs_real_3_t *surfxyz_norm;
  CS_MALLOC_HD(surfxyz_norm, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  c_id) {
    cs_real_t snorm = cs_math_3_norm(surfxyz_unnormed[c_id])+1.e-8;
    cs_real_t unsnorm = 1. / snorm;

    for (cs_lnum_t i = 0; i < 3; i++)
      surfxyz_norm[c_id][i] = surfxyz_unnormed[c_id][i] * unsnorm;
  });

  ctx.wait();

  /* Correction of surfxyz_norm at walls (Contact angle) */
  if (_contact_angle_choice == CS_VOF_CONTACT_ANGLE_STATIC) {
    _contact_angle_correction<CS_VOF_CONTACT_ANGLE_STATIC>(m,
                                                           mq,
                                                           ctx,
                                                           surfxyz_norm);
  }
  else if (_contact_angle_choice == CS_VOF_CONTACT_ANGLE_DYN) {
    _contact_angle_correction<CS_VOF_CONTACT_ANGLE_DYN>(m,
                                                        mq,
                                                        ctx,
                                                        surfxyz_norm);
  }

  /* Curvature Computation */
  cs_real_33_t *gradnxyz;
  CS_MALLOC_HD(gradnxyz, n_cells_ext, cs_real_33_t, cs_alloc_mode);

  cs_gradient_vector("grad_diff_void_grad",
                     gradient_type,
                     halo_type,
                     1,     /* inc */
                     eqp_volf->nswrgr,
                     eqp_volf->verbosity,
                     (cs_gradient_limit_t)eqp_volf->imligr,
                     eqp_volf->epsrgr,
                     eqp_volf->climgr,
                     &bc_coeffs_v_loc,
                     surfxyz_norm,
                     gweight,
                     nullptr,
                     gradnxyz);

  /* Reconstructions for curvature computation */
  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  c_id) {

    /* Initialization */
    cs_real_t curv = 0.;

    /* Loop on interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t face_id = cell_i_faces[cidx];
      const cs_lnum_t c_id_adj = c2c[cidx];
      const short int f_sgn = c2f_sgn[cidx];
      cs_real_t gradf[3];

      for (cs_lnum_t k = 0; k < 3; k++) {
        gradf[k] =         pond[face_id]  * surfxyz_norm[c_id][k]
                   + (1. - pond[face_id]) * surfxyz_norm[c_id_adj][k];

        for (cs_lnum_t l = 0; l < 3; l++)
          gradf[k] +=   0.5 * dofij[face_id][l]
                      * (gradnxyz[c_id][k][l] + gradnxyz[c_id_adj][k][l]);
      }

      cs_real_t flux =   f_sgn * i_face_surf[face_id]
                       * cs_math_3_dot_product(gradf, i_face_u_normal[face_id]);

      curv += flux;
    }

    /* Compute volumic surface tension */
    for (cs_lnum_t i = 0; i < 3; i++) {
      stf[c_id][i] = -cpro_surftens * surfxyz_unnormed[c_id][i] * curv;
    }

  });

  ctx.wait();

  CS_FREE_HD(surfxyz_norm);
  CS_FREE_HD(surfxyz_unnormed);
  CS_FREE_HD(gradnxyz);
  CS_FREE_HD(pvar);
  CS_FREE_HD(coefa);
  CS_FREE_HD(coefb);
  CS_FREE_HD(coefa_vec);
  CS_FREE_HD(coefb_vec);
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
  const cs_real_3_t *i_face_u_normal = (const cs_real_3_t *)mq->i_face_u_normal;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)m->i_face_cells;

  cs_dispatch_context ctx;

  /* Constant parameter */
  const cs_real_t cdrift = _vof_parameters.cdrift;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const cs_real_t *restrict i_volflux
    = cs_field_by_id(cs_field_get_key_int(CS_F_(void_f), kimasf))->val;

  cs_real_t *ipro_idriftf = cs_field_by_name("inner_drift_velocity_flux")->val;

  /* Check if field exists */
  if (ipro_idriftf == nullptr)
    bft_error(__FILE__, __LINE__, 0,_("error drift velocity not defined\n"));

  cs_real_3_t *voidf_grad;
  CS_MALLOC_HD(voidf_grad, n_cells_with_ghosts, cs_real_3_t, cs_alloc_mode);

  /* Compute the gradient of the void fraction */
  cs_field_gradient_scalar(CS_F_(void_f),
                           true,           // use_previous_t
                           1,              // inc
                           voidf_grad);

  /* Stabilization factor */
  const cs_real_t delta = 1.e-8/pow(tot_vol/n_g_cells,(1./3.));

  /* Compute the max of flux/Surf over the entire domain */
  cs_real_t maxfluxsurf;
  struct cs_reduce_max1r reducer;

  ctx.parallel_for_reduce
    (n_i_faces, maxfluxsurf, reducer,
     [=] CS_F_HOST_DEVICE (cs_lnum_t f_id, cs_real_t &res) {

    res = cs::abs(i_volflux[f_id])/i_face_surf[f_id];
  });

  ctx.wait();

  cs_parall_max(1, CS_REAL_TYPE, &maxfluxsurf);

  /* Compute the relative velocity at internal faces */

  ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    cs_lnum_t cell_id1 = i_face_cells[f_id][0];
    cs_lnum_t cell_id2 = i_face_cells[f_id][1];

    cs_real_t fluxfactor
      = cs::min(cdrift*std::abs(i_volflux[f_id])/i_face_surf[f_id],
                maxfluxsurf);

    cs_real_t gradface[3], normalface[3];

    for (cs_lnum_t idim = 0; idim < 3; idim++)
      gradface[idim] = 0.5 * (  voidf_grad[cell_id1][idim]
                              + voidf_grad[cell_id2][idim]);

    cs_real_t normgrad = cs_math_3_norm(gradface);

    for (cs_lnum_t idim = 0; idim < 3; idim++)
      normalface[idim] = gradface[idim] / (normgrad + delta);

    ipro_idriftf[f_id]
      =   i_face_surf[f_id] * fluxfactor
        * cs_math_3_dot_product(normalface, i_face_u_normal[f_id]);

  });

  ctx.wait();

  CS_FREE_HD(voidf_grad);
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
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t sum_type = ctx.get_parallel_for_i_faces_sum_type(m);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr) {
    if (m->halo != nullptr) {
      cs_halo_type_t halo_type = CS_HALO_STANDARD;
      cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;
      cs_gradient_type_by_imrgra(imrgra,
                                 &gradient_type,
                                 &halo_type);
      cs_halo_sync_var(m->halo, halo_type, pvar);
    }
  }
  else if (pvara == nullptr)
    pvara = (const cs_real_t *)pvar;

  const cs_real_t  *restrict _pvar = (pvar != nullptr) ? pvar : pvara;

  /*======================================================================
    Computation of the drift flux
    ======================================================================*/

  cs_field_t *vr = cs_field_by_name_try("drift_velocity");
  cs_field_t *idriftflux = cs_field_by_name_try("inner_drift_velocity_flux");
  cs_field_t *bdriftflux = cs_field_by_name_try("boundary_drift_velocity_flux");
  cs_real_t *ipro_idriftf = idriftflux->val;

  if (_vof_parameters.idrift == 1) {

    // FIXME Handle boundary terms bdriftflux
    cs_vof_deshpande_drift_flux(cs_glob_mesh, cs_glob_mesh_quantities);

  } else {

    const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

    /* Check if field exist */
    if (idriftflux == nullptr)
      bft_error(__FILE__, __LINE__, 0,_("error drift velocity not defined\n"));

    cs_real_3_t *cpro_vr = (cs_real_3_t *)vr->val;
    cs_real_t *cpro_bdriftf = bdriftflux->val;

    cs_field_bc_coeffs_t bc_coeffs_v_loc;
    cs_field_bc_coeffs_init(&bc_coeffs_v_loc);
    CS_MALLOC_HD(bc_coeffs_v_loc.a, 3*n_b_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(bc_coeffs_v_loc.b, 9*n_b_faces, cs_real_t, cs_alloc_mode);

    cs_real_3_t *coefav = (cs_real_3_t *)bc_coeffs_v_loc.a;
    cs_real_33_t *coefbv = (cs_real_33_t *)bc_coeffs_v_loc.b;

    /* Boundary coefficients */
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      for (cs_lnum_t ii = 0 ; ii < 3 ; ii++) {
        coefav[face_id][ii] = 0.;

        for (cs_lnum_t jj = 0 ; jj < 3 ; jj++) {
          coefbv[face_id][ii][jj] = 0.;
        }
        coefbv[face_id][ii][ii] = 1.;
      }
    });

    ctx.wait();

    int f_id = -1;
    int itypfl = 0;
    int iflmb0 = 1;
    int init = 1;
    int inc = 1;

    cs_mass_flux(m,
                 fvq,
                 f_id,
                 itypfl,
                 iflmb0,
                 init,
                 inc,
                 imrgra,
                 nswrgp,
                 static_cast<cs_gradient_limit_t>(imligp),
                 iwarnp,
                 epsrgp,
                 climgp,
                 nullptr, /* rom */
                 nullptr, /* romb */
                 (const cs_real_3_t *)cpro_vr,
                 &bc_coeffs_v_loc,
                 ipro_idriftf,
                 cpro_bdriftf);

    CS_FREE_HD(coefav);
    CS_FREE_HD(coefbv);
  }

  /*======================================================================
    Contribution from interior faces
    ======================================================================*/

  const int kiflux = cs_field_key_id("inner_flux_id");
  int i_flux_id = cs_field_get_key_int(CS_F_(void_f), kiflux);
  cs_real_t *i_flux = cs_field_by_id(i_flux_id)->val;

  const cs_real_t kdrift = _vof_parameters.kdrift;

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
    cs_lnum_t c_id1 = i_face_cells[face_id][0];
    cs_lnum_t c_id2 = i_face_cells[face_id][1];

    cs_real_t irvf = 0.;
    if (idriftflux != nullptr)
      irvf = ipro_idriftf[face_id];

    cs_real_t fluxij[2] = {0., 0.};

    cs_i_conv_flux(1,
                   1.,
                   0,
                   _pvar[c_id1],
                   _pvar[c_id2],
                   _pvar[c_id1]*(1.-_pvar[c_id2]),
                   _pvar[c_id1]*(1.-_pvar[c_id2]),
                   _pvar[c_id2]*(1.-_pvar[c_id1]),
                   _pvar[c_id2]*(1.-_pvar[c_id1]),
                   irvf,
                   1.,
                   1.,
                   fluxij);

    cs_i_diff_flux(1,
                   1.,
                   _pvar[c_id1],
                   _pvar[c_id2],
                   _pvar[c_id1],
                   _pvar[c_id2],
                   kdrift*(2.-_pvar[c_id1]-_pvar[c_id2])
                    / 2.*i_face_surf[face_id]/i_dist[face_id],
                   fluxij);

    if (c_id1 < n_cells)
      cs_dispatch_sum(&rhs[c_id1], -fluxij[0], sum_type);
    if (c_id2 < n_cells)
      cs_dispatch_sum(&rhs[c_id2],  fluxij[1], sum_type);

    /* store void fraction convection flux contribution */
    i_flux[face_id] += fluxij[0];
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the void fraction \f$ \alpha \f$ for the Volume of Fluid
 *        method (and hence for cavitating flows).
 *
 * This function solves:
 * \f[
 * \dfrac{\alpha^n - \alpha^{n-1}}{\Delta t}
 *     + \divs \left( \alpha^n \vect{u}^n \right)
 *     + \divs \left( \left[ \alpha^n
 *                           \left( 1 - \alpha^{n} \right)
 *                    \right] \vect{u^r}^n \right)
 *     = \dfrac{\Gamma_V \left( \alpha^{n-1}, p^n \right)}{\rho_v}
 * \f]
 * with \f$ \Gamma_V \f$ the eventual vaporization source term (Merkle model) in
 * case the cavitation model is enabled, \f$ \rho_v \f$ the reference gas
 * density and \f$ \vect{u^r} \f$ the drift velocity for the compressed
 * interface.
 *
 * \param[in]     iterns        Navier-Stokes iteration number
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_solve_void_fraction(int  iterns)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_real_t *cell_f_vol = fvq->cell_vol;

  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t rho1 = _vof_parameters.rho1;
  const cs_real_t rho2 = _vof_parameters.rho2;
  const cs_cavitation_parameters_t *cavitation_parameters
    = cs_get_glob_cavitation_parameters();

  cs_dispatch_context ctx;

  /* Initialization
     -------------- */

  const int i_vof_mass_transfer
    = (_vof_parameters.vof_model & CS_VOF_MERKLE_MASS_TRANSFER);
  const int isno2t = cs_glob_time_scheme->isno2t;

  cs_field_t *volf2 = cs_field_by_name("void_fraction");
  cs_equation_param_t *eqp_vol = cs_field_get_equation_param(volf2);
  int istat = eqp_vol->istat;

  cs_real_t *cvar_voidf = volf2->val;
  cs_real_t *cvara_voidf = volf2->val_pre;

  /* Implicitation in pressure of the
     vaporization/condensation model (cavitation) */

  cs_field_t *f_pr = CS_F_(p);
  cs_real_t *cvar_pr = nullptr;
  cs_real_t *cvara_pr = nullptr;
  if (i_vof_mass_transfer != 0 && cavitation_parameters->itscvi == 1) {
    cvar_pr = f_pr->val;
    cvara_pr = f_pr->val_pre;
  }

  /* Allocate temporary arrays */

  cs_real_t *i_visc, *b_visc, *smbrs, *rovsdt;
  CS_MALLOC_HD(rovsdt, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(i_visc, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_visc, n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(smbrs, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Allocate work arrays */
  cs_real_t *dpvar;
  CS_MALLOC_HD(dpvar, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Boundary conditions */
  cs_field_bc_coeffs_t *bc_coeffs_vof = volf2->bc_coeffs;

  /* Physical quantities */

  const int iflmas
    = cs_field_get_key_int(volf2, cs_field_key_id("inner_mass_flux_id"));
  cs_real_t *i_mass_flux_volf = cs_field_by_id(iflmas)->val;

  const int iflmab
    = cs_field_get_key_int(volf2, cs_field_key_id("boundary_mass_flux_id"));
  cs_real_t *b_mass_flux_volf = cs_field_by_id(iflmab)->val;

  /* Key id for clipping */
  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");

  /* Theta related to explicit source terms */

  cs_real_t *c_st_voidf = nullptr;
  const cs_real_t thets = cs_glob_time_scheme->thetsn;
  if (isno2t > 0) {
    const int kstprv = cs_field_key_id("source_term_prev_id");
    const int istprv = cs_field_get_key_int(volf2, kstprv);

    c_st_voidf = cs_field_by_id(istprv)->val;
  }

  /* Theta related to void fraction */
  //const cs_real_t thetv = eqp_vol->theta; /* UNUSED */

  /*  Arbitrary initialization :
   * - RHS at 0
   * - No diffusion for void fraction (visc = 1)
   *   convection flux at 0
   */

  const int kiflux = cs_field_key_id("inner_flux_id");
  const int kbflux = cs_field_key_id("boundary_flux_id");

  const int icflux_id = cs_field_get_key_int(volf2, kiflux);
  const int bcflux_id = cs_field_get_key_int(volf2, kbflux);

  cs_real_t *icflux = cs_field_by_id(icflux_id)->val;
  cs_real_t *bcflux = cs_field_by_id(bcflux_id)->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    smbrs[c_id] = 0.;
  });

  ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
    i_visc[face_id] = 1.;
    icflux[face_id] = 0.;
  });

  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
    b_visc[face_id] = 1.;
    bcflux[face_id] = 0.;
  });

  ctx.wait();

  /* Preliminary computations
     ------------------------ */

  /* Update the cavitation source term with pressure increment
     if it has been implicited in pressure at correction step,
     in order to ensure mass conservation. */

  cs_real_t *gamcav = nullptr, *dgdpca = nullptr;
  if (i_vof_mass_transfer != 0) {
    gamcav = cs_field_by_name("model:cavitation_gamma")->val;
    dgdpca = cs_field_by_name("model:cavitation_st_dgdp")->val;
  }

  if (i_vof_mass_transfer != 0 && cavitation_parameters->itscvi == 1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  c_id) {
      gamcav[c_id] += dgdpca[c_id] * (cvar_pr[c_id] - cvara_pr[c_id]);
    });

    ctx.wait();
  }

  /* Compute the limiting time step to satisfy min/max principle.
     Only if a source term is accounted for. */

  cs_real_t dtmaxg = HUGE_VAL;

  if (i_vof_mass_transfer != 0) {

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t dtmaxl;
      if (gamcav[c_id] < 0.0)
        dtmaxl = - rho2 * cvara_voidf[c_id] / gamcav[c_id];
      else
        dtmaxl = rho1 * (1.0 - cvara_voidf[c_id]) / gamcav[c_id];

      dtmaxg = cs::min(dtmaxl, dtmaxg);

    }

    cs_parall_min(1, CS_REAL_TYPE, &dtmaxg);

    if (dt[0] > dtmaxg) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("@\n"
           "@ @@ Warning: Void fraction resolution\n"
           "@    ========\n"
           "@  The current time step is too large to ensure the min/max\n"
           "@     principle on void fraction.\n"
           "@\n"
           "@  The current time step is %13.5e while\n"
           "@     the maximum admissible value is %13.5e\n"
           "@\n"
           "@  Clipping on void fraction should occur and\n"
           "@     mass conservation is lost.\n"
           "@\n"),
         dt[0], dtmaxg);
    }

  }

  /* Visualization of divu (only for advanced analysis purpose,
     not used in the source terms hereafter) */

  cs_field_t *vel_div = cs_field_by_name_try("algo:velocity_divergence");
  cs_real_t *divu = nullptr, *_divu = nullptr;
  if (vel_div != nullptr)
    divu = vel_div->val;
  else {
    /* Allocation */
    CS_MALLOC_HD(_divu, n_cells_ext, cs_real_t, cs_alloc_mode);
    divu = _divu;
  }

  cs_divergence(mesh,
                1, /* init */
                i_mass_flux_volf,
                b_mass_flux_volf,
                divu);

  CS_FREE_HD(_divu);

  /* Construct the system to solve
     ----------------------------- */

  /* Source terms */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* Cavitation source term (explicit) */
    if (i_vof_mass_transfer != 0)
      smbrs[c_id] += cell_f_vol[c_id] * gamcav[c_id] / rho2;

    /* Source terms assembly for cs_equation_iterative_solve_scalar */

    /* If source terms are extrapolated over time */
    if (isno2t > 0) {
      const cs_real_t tsexp = c_st_voidf[c_id];
      c_st_voidf[c_id] = smbrs[c_id];
      smbrs[c_id] = - thets * tsexp + (1.0 + thets) * smbrs[c_id];
    }

    /* Source term linked with the non-conservative form of convection term
       in cs_equation_iterative_solve_scalar (always implicited)
       FIXME set imasac per variable? Here it could be set to 0
       (or Gamma(1/rho2 - 1/rho1) for cavitation) and divu not used.
       Note: we prefer to be not perfectly conservative (up to the precision
       of the pressure solver, but that allows to fullfill the min/max principle
       on alpha */

    if (i_vof_mass_transfer != 0) {
      /* Should be for the conservative form */

      smbrs[c_id] += - cell_f_vol[c_id] * gamcav[c_id]
                     * (1.0/rho2 - 1.0/rho1) * cvara_voidf[c_id];

      rovsdt[c_id] = cell_f_vol[c_id] * gamcav[c_id] * (1.0/rho2 - 1.0/rho1);
    }
    else {
      rovsdt[c_id] = 0;
    }

    /* Unteady term */
    if (istat > 0) {
      rovsdt[c_id] += cell_f_vol[c_id] / dt[c_id];
    }

  });

  ctx.wait();

  if (cs_glob_vof_parameters->idrift > 0) {

    cs_vof_drift_term(eqp_vol->imrgra,
                      eqp_vol->nswrgr,
                      eqp_vol->imligr,
                      eqp_vol->verbosity,
                      eqp_vol->epsrgr,
                      eqp_vol->climgr,
                      cvar_voidf,
                      cvara_voidf,
                      smbrs);
  }

  /* Solving
     ------- */

  /* Solving void fraction */
  int iescap = 0;
  int imucpp = 0;
  /* All boundary convective flux with upwind */
  int icvflb = 0;
  cs_real_t normp = -1.0;

  cs_equation_param_t eqp_loc = *eqp_vol;

  eqp_loc.icoupl = -1;
  eqp_loc.idifft = -1;
  eqp_loc.iwgrec = 0; /* Warning, may be overwritten if a field */
  eqp_loc.blend_st = 0; /* Warning, may be overwritten if a field */

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     iterns,
                                     volf2->id,
                                     nullptr,
                                     iescap,
                                     imucpp,
                                     normp,
                                     &eqp_loc,
                                     cvara_voidf,
                                     cvara_voidf,
                                     bc_coeffs_vof,
                                     i_mass_flux_volf, b_mass_flux_volf,
                                     i_visc, b_visc,
                                     i_visc, b_visc,
                                     nullptr, /* viscel */
                                     nullptr, nullptr, /* weighf, weighb */
                                     icvflb,
                                     nullptr, /* icvfli */
                                     rovsdt,
                                     smbrs,
                                     cvar_voidf,
                                     dpvar,
                                     nullptr,  /* xcpp */
                                     nullptr); /* eswork */

  /* Clipping: only if min/max principle is not satisfied for cavitation
     ------------------------------------------------------------------- */

  if (  (i_vof_mass_transfer != 0 && dt[0] > dtmaxg)
      || i_vof_mass_transfer == 0) {

    /* Compute min and max */
    struct cs_data_2r rd;
    struct cs_reduce_min1r_max1r reducer;

    ctx.parallel_for_reduce
      (n_cells, rd, reducer,
       [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_data_2r &res) {
      res.r[0] = cvar_voidf[c_id];
      res.r[1] = cvar_voidf[c_id];
    });

    ctx.wait();

    /* Get the min and max clipping */

    const cs_real_t scminp = cs_field_get_key_double(volf2, kscmin);
    const cs_real_t scmaxp = cs_field_get_key_double(volf2, kscmax);
    const int kclipp = cs_field_key_id("clipping_id");

    /* Postprocess clippings ? */

    const int clip_voidf_id = cs_field_get_key_int(volf2, kclipp);
    cs_real_t *voidf_clipped = nullptr;

    if (clip_voidf_id >= 0) {
      voidf_clipped = cs_field_by_id(clip_voidf_id)->val;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        voidf_clipped[c_id] = 0.;
      });
      ctx.wait();
    }

    struct cs_data_2i rd_sum;
    rd_sum.i[0] = 0;
    rd_sum.i[1] = 0;
    struct cs_reduce_sum2i reducer_sum;

    if (scmaxp > scminp) {

      ctx.parallel_for_reduce
        (n_cells, rd_sum, reducer_sum,
         [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_data_2i &res) {

        res.i[0] = 0, res.i[1] = 0;

        if (cvar_voidf[c_id] > scmaxp) {
          res.i[0] = 1;

          if (clip_voidf_id >= 0)
            voidf_clipped[c_id] = cvar_voidf[c_id] - scmaxp;

          cvar_voidf[c_id] = scmaxp;
        }

        if (cvar_voidf[c_id] < scminp) {
          res.i[1] = 1;

          if (clip_voidf_id >= 0)
            voidf_clipped[c_id] = cvar_voidf[c_id] - scminp;

          cvar_voidf[c_id] = scminp;
        }
      });

      ctx.wait();
    }

    cs_log_iteration_clipping_field(volf2->id,
                                    rd_sum.i[1],
                                    rd_sum.i[0],
                                    &rd.r[0],
                                    &rd.r[1],
                                    &rd_sum.i[1],
                                    &rd_sum.i[0]);

  }

  /* Free memory */

  CS_FREE_HD(rovsdt);
  CS_FREE_HD(i_visc);
  CS_FREE_HD(b_visc);
  CS_FREE_HD(smbrs);
  CS_FREE_HD(dpvar);
}

/*----------------------------------------------------------------------------
 *!
 * \brief Provide write access to cavitation parameters structure.
 */
/*----------------------------------------------------------------------------*/

cs_cavitation_parameters_t *
cs_get_glob_cavitation_parameters(void)
{
  return &_cavit_parameters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the vaporization source term
 * \f$ \Gamma_V \left(\alpha, p\right) = m^+ + m^- \f$ using the
 * Merkle model:
 *
 * \f[
 * m^+ = -\dfrac{C_{prod} \rho_l \min \left( p-p_V,0 \right)\alpha(1-\alpha)}
 *              {0.5\rho_lu_\infty^2t_\infty},
 * \f]
 * \f[
 * m^- = -\dfrac{C_{dest} \rho_v \max \left( p-p_V,0 \right)\alpha(1-\alpha)}
 *              {0.5\rho_lu_\infty^2t_\infty},
 * \f]
 * with \f$ C_{prod}, C_{dest} \f$ empirical constants,
 * \f$ t_\infty=l_\infty/u_\infty \f$ a reference time scale and \f$p_V\f$
 * the reference saturation pressure.
 * \f$ l_\infty \f$, \f$ u_\infty \f$ and \f$p_V\f$ may be provided by
 * the user (user function).
 *
 * Note that the r.h.s. of the void fraction transport equation is
 * \f$ \Gamma_V/\rho_v \f$.
 *
 * \param[in]  pressure  Pressure array
 * \param[in]  voidf     Void fraction array
 */
/*----------------------------------------------------------------------------*/

void
cs_cavitation_compute_source_term(const cs_real_t  pressure[],
                                  const cs_real_t  voidf[])
{
  if (!(cs_glob_vof_parameters->vof_model & CS_VOF_MERKLE_MASS_TRANSFER))
    return;

  const cs_vof_parameters_t *vof = cs_glob_vof_parameters;
  const cs_cavitation_parameters_t *cav = cs_glob_cavitation_parameters;

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t n_cells = mesh->n_cells;

  cs_dispatch_context ctx;

  // Merkle model

  cs_real_t *gamcav = cs_field_by_name("model:cavitation_gamma")->val;
  cs_real_t *dgdpca = cs_field_by_name("model:cavitation_st_dgdp")->val;

  cs_real_t tinf = cav->linf / cav->uinf;

  cs_real_t cond =   (cav->cdest * vof->rho2)
                   / (0.5 * vof->rho1 * cav->uinf * cav->uinf * tinf);
  cs_real_t cvap =   (cav->cprod * vof->rho1)
                   / (0.5 * vof->rho1 * cav->uinf * cav->uinf * tinf);

  const cs_real_t presat = cav->presat;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    const cs_real_t w = voidf[c_id] * (1. - voidf[c_id]);
    cs_real_t condens = -cond * cs::max(0., pressure[c_id] - presat) * w;
    cs_real_t vaporis = -cvap * cs::min(0., pressure[c_id] - presat) * w;

    gamcav[c_id] = condens + vaporis;

    if (gamcav[c_id] < 0)
      dgdpca[c_id] = - cond * voidf[c_id] * (1. - voidf[c_id]);
    else
      dgdpca[c_id] = - cvap * voidf[c_id] * (1. - voidf[c_id]);
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the contact angle mode
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_contact_angle_set
(
  const cs_vof_contact_angle_t choice /*!<[in] Contact angle model to set. */
)
{
  /* Sanity check */
  assert(choice < CS_VOF_N_CONTACT_ANGLE_TYPES);

  _contact_angle_choice = choice;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
