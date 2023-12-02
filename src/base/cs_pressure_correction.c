/*============================================================================
 * Pressure correction.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_air_props.h"
#include "cs_cdo_headers.h"
#include "cs_ale.h"
#include "cs_array.h"
#include "cs_atmo.h"
#include "cs_balance.h"
#include "cs_blas.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_cf_thermo.h"
#include "cs_convection_diffusion.h"
#include "cs_divergence.h"
#include "cs_equation_iterative_solve.h"
#include "cs_face_viscosity.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_halo.h"
#include "cs_lagr.h"
#include "cs_log.h"
#include "cs_matrix_building.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_porous_model.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_sat_coupling.h"
#include "cs_sles_default.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"
#include "cs_volume_mass_injection.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_pressure_correction.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_pressure_correction.c
        Pressure correction step.
*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_PRESSURE_CORRECTION_CDO_DBG  0

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_prcdo[] =
  " Stop execution. The structure related to the pressure correction is"
  " empty.\n Please check your settings.\n";

/*============================================================================
 * Global variables
 *============================================================================*/

static bool cs_pressure_correction_cdo_active = false;

static cs_pressure_correction_cdo_t *cs_pressure_correction_cdo = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Poisson equation resolution for hydrostatic pressure:
 *         \f$ \divs ( \grad P ) = \divs ( f ) \f$
 *
 * \param[in]      m                pointer to glob mesh
 * \param[in]      mq               pointer to glob mesh quantiites
 * \param[out]     indhyd           indicateur de mise a jour de cvar_hydro_pres
 * \param[in]      iterns           Navier-Stokes iteration number
 * \param[in]      frcxt            external force generating hydrostatic pressure
 * \param[in]      dfrcxt           external force increment
 *                                  generating hydrostatic pressure
 * \param[out]     cvar_hydro_pres  hydrostatic pressure increment
 * \param[in]      iflux            work array
 * \param[in]      bflux            work array
 * \param[in,out]  i_visc           work array
 * \param[in,out]  b_visc           work array
 * \param[in,out]  dam              work array
 * \param[in,out]  xam              work array
 * \param[in,out]  dphi             work array
 * \param[in,out]  rhs              work array
 */
/*----------------------------------------------------------------------------*/

static void
_hydrostatic_pressure_compute(const cs_mesh_t       *m,
                              cs_mesh_quantities_t  *mq,
                              int                   *indhyd,
                              int                    iterns,
                              const cs_real_t        frcxt[][3],
                              const cs_real_t        dfrcxt[][3],
                              cs_real_t              cvar_hydro_pres[],
                              cs_real_t              iflux[],
                              cs_real_t              bflux[],
                              cs_real_t              i_visc[],
                              cs_real_t              b_visc[],
                              cs_real_t              dam[],
                              cs_real_t              xam[],
                              cs_real_t              dphi[],
                              cs_real_t              rhs[])
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *distb = mq->b_dist;

  const int ksinfo = cs_field_key_id("solving_info");
  cs_field_t *f = cs_field_by_name_try("hydrostatic_pressure");
  cs_solving_info_t *sinfo = cs_field_get_key_struct_ptr(f, ksinfo);
  const cs_equation_param_t *eqp_pr = cs_field_get_equation_param_const(f);

  if (iterns < 2)
    sinfo->n_it = 0;

  cs_real_t *coefap = f->bc_coeffs->a;
  cs_real_t *coefbp = f->bc_coeffs->b;
  cs_real_t *cofafp = f->bc_coeffs->af;
  cs_real_t *cofbfp = f->bc_coeffs->bf;

  /* Check for variation of the hydrostatic pressure at outlet
   *
   * We check if the source term has changed. We exit directly
   * if we do not have standard outlets.
   * The precisiton for tests if more or less arbitrary. */

  int ical = 0;
  const cs_real_t precab = 1e2*cs_math_epzero;
  const cs_real_t precre = sqrt(cs_math_epzero);
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t rnrmf = cs_math_3_square_norm(frcxt[c_id]);
    const cs_real_t rnrmdf = cs_math_3_square_norm(dfrcxt[c_id]);
    if ((rnrmdf >= precre*rnrmf) && (rnrmdf >= precab))
      ical = 1;
  }

  cs_parall_sum(1, CS_INT_TYPE, &ical);
  if ((ical == 0) && (cs_glob_atmo_option->open_bcs_treatment == 0)) {
    *indhyd = 0;
    return;
  }

  const cs_time_step_t *ts = cs_glob_time_step;
  if (ts->nt_cur%cs_glob_log_frequency == 0 || eqp_pr->verbosity > 0)
    bft_printf("  Hydrostatic pressure computation:\n"
               "    updating the Dirichlets at the end"
               " (_hydrostatic_pressure_compute)\n");

  *indhyd = 1;

  cs_real_3_t *next_fext = NULL;
  BFT_MALLOC(next_fext, n_cells_ext, cs_real_3_t);

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const int c_act = cs_mesh_quantities_cell_is_active(mq, c_id);
    for (cs_lnum_t ii = 0; ii < 3; ii++)
      next_fext[c_id][ii] = frcxt[c_id][ii]*c_act + dfrcxt[c_id][ii];
  }

  cs_mesh_sync_var_vect((cs_real_t*)next_fext);

  /* Prepare matrix and boundary conditions
     -------------------------------------- */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    const cs_real_t qimp = 0;
    const cs_real_t hint = 1.0/distb[face_id];
    cs_boundary_conditions_set_neumann_scalar(&coefap[face_id],
                                              &cofafp[face_id],
                                              &coefbp[face_id],
                                              &cofbfp[face_id],
                                              qimp,
                                              hint);
  }

  cs_real_t *rovsdt = NULL, *viscce = NULL;

  BFT_MALLOC(viscce, n_cells_ext, cs_real_t);
  BFT_MALLOC(rovsdt, n_cells_ext, cs_real_t);

  cs_array_real_fill_zero(n_cells, rovsdt);
  cs_array_set_value_real(n_cells, 1, 1, viscce);

  cs_face_viscosity(m,
                    mq,
                    eqp_pr->imvisf,
                    viscce,
                    i_visc,
                    b_visc);

  cs_matrix_wrapper_scalar(eqp_pr->iconv,
                           eqp_pr->idiff,
                           0,
                           1,
                           eqp_pr->thetav,
                           0,
                           coefbp,
                           cofbfp,
                           rovsdt,
                           iflux,
                           bflux,
                           i_visc,
                           b_visc,
                           NULL,
                           dam,
                           xam);

  /* Compute right hand side
     ----------------------- */

  /* Compute div(f_ext^n+1) */
  cs_ext_force_flux(m,
                    mq,
                    1,
                    eqp_pr->nswrgr,
                    next_fext,
                    cofbfp,
                    iflux,
                    bflux,
                    i_visc,
                    b_visc,
                    viscce,
                    viscce,
                    viscce);

  cs_real_t *div_fext = NULL;
  BFT_MALLOC(div_fext, n_cells_ext, cs_real_t);

  cs_divergence(m, 1, iflux, bflux, div_fext);

  cs_real_t residual = 0;
  const cs_real_t rnorm = sqrt(cs_gdot(n_cells, div_fext, div_fext));

  /* Loops on non-orthogonalities (resolution)
     ----------------------------------------- */

  cs_array_real_fill_zero(n_cells, dphi);

  /* Reconstruction loop */

  for (int isweep = 1; isweep <= eqp_pr->nswrsm; isweep++) {

    /* Update the right hand side and update the residual
     *  rhs^{k+1} = - div(fext^n+1) - D(1, p_h^{k+1})
     *--------------------------------------------------- */

    cs_diffusion_potential(-1,              /* f_id */
                           m,
                           mq,
                           1,               /* init */
                           1,               /* inc  */
                           eqp_pr->imrgra,
                           eqp_pr->nswrgr,
                           eqp_pr->imligr,
                           1,               /* iphydp */
                           eqp_pr->iwgrec,
                           eqp_pr->verbosity,
                           eqp_pr->epsrgr,
                           eqp_pr->climgr,
                           next_fext,
                           cvar_hydro_pres,
                           coefap,
                           coefbp,
                           cofafp,
                           cofbfp,
                           i_visc,
                           b_visc,
                           viscce,
                           rhs);

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rhs[c_id] = -div_fext[c_id] - rhs[c_id];

    /* Convergence test */

    residual = sqrt(cs_gdot(n_cells, rhs, rhs));

    if (eqp_pr->verbosity > 1)
      bft_printf(_("hydrostatic_p: sweep = %d RHS residual = %10.14le"
                   " norm = %10.14le\n"),
                 isweep, residual, rnorm);

    if (isweep == 1)
      sinfo->rhs_norm = residual;

    if (residual <= eqp_pr->epsrsm*rnorm)
      break;

    /* Solving on the increment */

    int n_iter;
    cs_real_t ressol = residual;
    cs_array_real_fill_zero(n_cells, dphi);

    cs_sles_solve_native(f->id,
                         NULL,
                         true,  /* symmetric */
                         1,     /* diag_block_size */
                         1,     /* extra_diag_block_size */
                         dam,
                         xam,
                         eqp_pr->epsilo,
                         rnorm,
                         &n_iter,
                         &ressol,
                         rhs,
                         dphi);

    sinfo->n_it += n_iter;

    cs_axpy(n_cells, 1.0, dphi, cvar_hydro_pres);

  }

  if (eqp_pr->verbosity > 1)
    bft_printf("@\n"
               "@ @@ Warning: hydrostatic_p (hydrostatic pressure step)\n"
               "@    ========\n"
               "@  Maximum number of iterations %d reached\n"
               "@\n", eqp_pr->nswrsm);

  /* For logging */

  sinfo->res_norm = 0;
  if (fabs(rnorm) > 0)
    sinfo->res_norm = residual/rnorm;

  /* Free Memory and solver setup */

  cs_sles_free_native(f->id, NULL);

  BFT_FREE(viscce);
  BFT_FREE(rovsdt);
  BFT_FREE(div_fext);
  BFT_FREE(next_fext);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the pressure correction step of the Navier-Stokes equations
 *        for incompressible or slightly compressible flows.
 *
 * This function solves the following Poisson equation on the pressure:
 * \f[
 *     D \left( \Delta t, \delta p \right) =
 * \divs \left( \rho \vect{\widetilde{u}}\right)
 *     - \Gamma^n
 *     + \dfrac{\rho^n - \rho^{n-1}}{\Delta t}
 * \f]
 * The mass flux is then updated as follows:
 * \f[
 *  \dot{m}^{n+1}_\ij = \dot{m}^{n}_\ij
 *                    - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
 * \f]
 *
 * \remark:
 * - an iterative process is used to solve the Poisson equation.
 * - if the arak coefficient is set to 1, the the Rhie & Chow filter is
 *   activated.
 *
 * Please refer to the
 * <a href="../../theory.pdf#resopv"><b>resopv</b></a>
 * section of the theory guide for more information.
 *
 * \param[in]       iterns    Navier-Stokes iteration number
 * \param[in]       nfbpcd    number of faces with condensation source term
 * \param[in]       ncmast    number of cells with condensation source terms
 * \param[in]       ifbpcd    index of faces with condensation source term
 * \param[in]       ltmast    list of cells with condensation source terms
 *                            (1 to n numbering)
 * \param[in]       isostd    indicator of standard outlet and index
 *                            of the reference outlet face
 * \param[in]       vel       velocity
 * \param[in, out]  da_uu     velocity matrix
 * \param[in]       coefav    boundary condition array for the variable
 *                            (explicit part)
 * \param[in]       coefbv    boundary condition array for the variable
 *                            (implicit part)
 * \param[in]       coefa_dp  boundary conditions for the pressure increment
 * \param[in]       coefb_dp  boundary conditions for the pressure increment
 * \param[in]       spcond    variable value associated to the condensation
 *                            source term (for ivar=ipr, spcond is the
 *                            flow rate
 *                            \f$ \Gamma_{s,cond}^n \f$)
 * \param[in]       svcond    variable value associated to the condensation
 *                            source term (for ivar=ipr, svcond is the flow rate
 *                            \f$ \Gamma_{v, cond}^n \f$)
 * \param[in]       frcxt     external forces making hydrostatic pressure
 * \param[in]       dfrcxt    variation of the external forces
 *                            composing the hydrostatic pressure
 * \param[in]       i_visc    visc*surface/dist aux faces internes
 * \param[in]       b_visc    visc*surface/dist aux faces de bord
 */
/*----------------------------------------------------------------------------*/

static void
_pressure_correction_fv(int        iterns,
                        cs_lnum_t  nfbpcd,
                        cs_lnum_t  ncmast,
                        cs_lnum_t  ifbpcd[],
                        cs_lnum_t  ltmast[],
                        const int  isostd[],
                        cs_real_t  vel[restrict][3],
                        cs_real_t  da_uu[restrict][6],
                        cs_real_t  coefav[restrict][3],
                        cs_real_t  coefbv[restrict][3][3],
                        cs_real_t  coefa_dp[restrict],
                        cs_real_t  coefb_dp[restrict],
                        cs_real_t  spcond[restrict],
                        cs_real_t  svcond[restrict],
                        cs_real_t  frcxt[restrict][3],
                        cs_real_t  dfrcxt[restrict][3],
                        cs_real_t  i_visc[restrict],
                        cs_real_t  b_visc[restrict])
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_dist = (const cs_real_t *)fvq->b_dist;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_t *restrict i_f_face_surf = fvq->i_f_face_surf;
  const cs_real_t *restrict b_f_face_surf = fvq->b_f_face_surf;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;

  const int *bc_type = cs_glob_bc_type;

  cs_field_t *f_p = CS_F_(p);
  const cs_field_t *f_vel = CS_F_(vel);

  assert((cs_real_t *)vel == f_vel->val);

  const cs_equation_param_t *eqp_u
    = cs_field_get_equation_param_const(f_vel);
  const cs_equation_param_t *eqp_p
    = cs_field_get_equation_param_const(f_p);

  const cs_velocity_pressure_param_t  *vp_param
    = cs_glob_velocity_pressure_param;
  const cs_velocity_pressure_model_t  *vp_model
    = cs_glob_velocity_pressure_model;
  const cs_vof_parameters_t *vof_parameters = cs_glob_vof_parameters;
  const cs_cavitation_parameters_t *cavitation_parameters
    = cs_get_glob_cavitation_parameters();

  const int idilat = vp_model->idilat;
  const int idtvar = cs_glob_time_step_options->idtvar;

  const int compressible_flag
    = cs_glob_physical_model_flag[CS_COMPRESSIBLE];
  const int open_bcs_flag = cs_glob_atmo_option->open_bcs_treatment;

  const cs_fluid_properties_t  *fluid_props = cs_glob_fluid_properties;

  cs_real_t *restrict dt = CS_F_(dt)->val;

  cs_real_t *cvar_hydro_pres = NULL, *cvar_hydro_pres_prev = NULL;
  cs_real_t *c_visc = NULL;
  cs_real_6_t *vitenp = NULL;
  cs_real_t *taui = NULL, *taub = NULL;

  cs_field_t  *f_hp = cs_field_by_name_try("hydrostatic_pressure");
  if (f_hp != NULL) {
    cvar_hydro_pres = f_hp->vals[0];
    cvar_hydro_pres_prev = f_hp->vals[1];
  }

  /* Allocate temporary arrays */

  cs_real_t *dam, *xam, *rhs, *res;
  BFT_MALLOC(dam, n_cells_ext, cs_real_t);
  BFT_MALLOC(xam, m->n_i_faces, cs_real_t);
  CS_MALLOC_HD(rhs, n_cells_ext, cs_real_t, cs_alloc_mode);
  BFT_MALLOC(res, n_cells_ext, cs_real_t);

  cs_real_t *phia, *iflux, *bflux, *dphi;
  BFT_MALLOC(phia, n_cells_ext, cs_real_t);
  BFT_MALLOC(iflux, m->n_i_faces, cs_real_t);
  BFT_MALLOC(bflux, m->n_b_faces, cs_real_t);
  BFT_MALLOC(dphi, n_cells_ext, cs_real_t);

  cs_real_3_t *wrk;
  BFT_MALLOC(wrk, n_cells_ext, cs_real_3_t);
  cs_real_3_t *wrk2;
  BFT_MALLOC(wrk2, n_cells_ext, cs_real_3_t);

  cs_real_t *adxk = NULL, *adxkm1 = NULL, *dphim1 = NULL, *rhs0 = NULL;
  if (eqp_p->iswdyn > 0) {
    BFT_MALLOC(adxk, n_cells_ext, cs_real_t);
    BFT_MALLOC(adxkm1, n_cells_ext, cs_real_t);
    BFT_MALLOC(dphim1, n_cells_ext, cs_real_t);
    BFT_MALLOC(rhs0, n_cells_ext, cs_real_t);
  }

  cs_field_t *f_iddp = cs_field_by_name("pressure_increment");
  cs_real_t *phi = f_iddp->val;

  /* Diffusive flux Boundary conditions for delta P */
  cs_real_t *coefaf_dp = f_iddp->bc_coeffs->af;
  cs_real_t *coefbf_dp = f_iddp->bc_coeffs->bf;

  /* Associate pointers to pressure diffusion coefficients */
  c_visc = dt;
  if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
    vitenp = (cs_real_6_t *)(cs_field_by_name("dttens")->val);

  /* Index of the field */
  cs_solving_info_t *sinfo
    = cs_field_get_key_struct_ptr(f_p,
                                  cs_field_key_id("solving_info"));

  cs_field_t *f_weight = NULL;
  if (eqp_p->iwgrec == 1) {
    /* Weighting field for gradient */
    int kwgrec = cs_field_key_id("gradient_weighting_id");
    f_weight = cs_field_by_id(cs_field_get_key_int(f_p, kwgrec));
  }

  cs_real_t *cpro_divu = NULL, *_cpro_divu = NULL;
  cs_field_t *f_divu = cs_field_by_name_try("predicted_vel_divergence");
  if (f_divu != NULL)
    cpro_divu = f_divu->val;
  else {
    BFT_MALLOC(_cpro_divu, n_cells_ext, cs_real_t);
    cpro_divu = _cpro_divu;
  }

  /* Boundary conditions */

  cs_real_t *coefa_p = f_p->bc_coeffs->a;
  cs_real_t *coefb_p = f_p->bc_coeffs->b;
  cs_real_t *coefaf_p = f_p->bc_coeffs->af;
  cs_real_t *coefbf_p = f_p->bc_coeffs->bf;

  /* Physical quantities */

  cs_real_t *crom = NULL;
  const cs_real_t *brom = NULL;

  cs_real_t *crom_eos = CS_F_(rho)->val;
  const cs_real_t *croma = NULL;
  if (vp_param->icalhy == 1 || idilat > 1 || fluid_props->irovar) {
    croma = CS_F_(rho)->val_pre;
  }
  const cs_real_t *brom_eos = CS_F_(rho_b)->val;
  const cs_real_t *broma = NULL;
  if (fluid_props->irovar)
    broma = CS_F_(rho_b)->val_pre;

  /* Time-interpolated density */

  cs_real_t *cpro_rho_tc = NULL, *bpro_rho_tc = NULL;

  if (fluid_props->irovar && (   idilat > 1
                              || vof_parameters->vof_model > 0
                              || compressible_flag == 3)) {

    cs_real_t *cpro_rho_mass = cs_field_by_name("density_mass")->val;
    cs_real_t *bpro_rho_mass = cs_field_by_name("boundary_density_mass")->val;

    /* Staggered in time velocity and pressure */
    if (eqp_u->thetav < 1 && iterns > 1 && vp_param->itpcol == 0) {
      BFT_MALLOC(cpro_rho_tc, n_cells_ext, cs_real_t);
      BFT_MALLOC(bpro_rho_tc, m->n_b_faces, cs_real_t);

      const cs_real_t thetav =  eqp_u->thetav;
      const cs_real_t one_m_thetav = 1.0 - thetav;

      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
        cpro_rho_tc[c_id] =  thetav * cpro_rho_mass[c_id]
                            + one_m_thetav * croma[c_id];
      }

      crom = cpro_rho_tc;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        bpro_rho_tc[f_id] =   thetav * bpro_rho_mass[f_id]
                            + one_m_thetav * broma[f_id];
      }

      brom = bpro_rho_tc;
    }
    else {
      crom = cpro_rho_mass;
      brom = bpro_rho_mass;
    }

  }

  /* Weakly variable density algo. (idilat <=1) or constant density */

  else {
    crom = crom_eos;
    brom = brom_eos;
  }

  cs_real_t *cvar_pr = f_p->vals[0];
  cs_real_t *cvara_pr = f_p->vals[1];

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");

  cs_real_t *imasfl = cs_field_by_id(cs_field_get_key_int(f_p, kimasf))->val;
  cs_real_t *bmasfl = cs_field_by_id(cs_field_get_key_int(f_p, kbmasf))->val;

  /* Solving options */

  const int isym = (eqp_p->iconv > 0) ? 2 : 1;
  const bool symmetric = (isym == 1) ? true : false;

  /* Using to VOF algo,
     the pressure correction is done through the volumetric flux (that is
     the convective flux of the void fraction), not the mass flux */

  if (vof_parameters->vof_model > 0) {
    cs_field_t *f_vf = cs_field_by_name("void_fraction");
    imasfl = cs_field_by_id(cs_field_get_key_int(f_vf, kimasf))->val;
    bmasfl = cs_field_by_id(cs_field_get_key_int(f_vf, kbmasf))->val;
  }

  const int i_vof_mass_transfer = (  vof_parameters->vof_model
                                   & CS_VOF_MERKLE_MASS_TRANSFER);

  if ((eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) && vp_param->rcfact == 0) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t j = 0; j < 6; j++)
        da_uu[c_id][j] = vitenp[c_id][j];
    }
    cs_mesh_sync_var_sym_tens(da_uu);
  }

  /* Calculation of dt/rho */

  cs_real_t *xdtsro = NULL;
  cs_real_6_t *tpusro = NULL;

  if (vof_parameters->vof_model > 0 || idilat == 4) {

    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION) {
      BFT_MALLOC(xdtsro, n_cells_ext, cs_real_t);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        xdtsro[c_id] = dt[c_id]/crom[c_id];

      cs_halo_sync_var(m->halo, CS_HALO_STANDARD, xdtsro);
    }
    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) {
      BFT_MALLOC(tpusro, n_cells_ext, cs_real_6_t);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t drom = 1 / crom[c_id];
        for (cs_lnum_t j = 0; j < 6; j++)
          tpusro[c_id][j] = vitenp[c_id][j] * drom;
      }
      cs_mesh_sync_var_sym_tens(tpusro);

      vitenp = tpusro;
    }

    /* Associate pointers to pressure diffusion coefficient */

    c_visc = xdtsro;

  }

  if (vp_param->staggered == 1) {

    const cs_real_t *hli = cs_field_by_name("inner_face_head_loss")->val;
    const cs_real_t *hlb = cs_field_by_name("boundary_face_head_loss")->val;

    BFT_MALLOC(taui, n_i_faces, cs_real_t);
    BFT_MALLOC(taub, n_b_faces, cs_real_t);

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id_0 = i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = i_face_cells[f_id][1];
      cs_real_t dt_f = 0.5 * (dt[c_id_0] + dt[c_id_1]);
      taui[f_id] = dt_f / (1. + hli[f_id] * dt_f);
    }
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, taui);
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      cs_lnum_t c_id = b_face_cells[f_id];
      taub[f_id] = dt[c_id] / (1. + hlb[f_id] * dt[c_id]);
    }

  }

  /* Compute an approximated pressure increment if needed,
   * that is when there are buoyancy terms (gravity and variable density)
   * with a free outlet.
   * ==================================================================== */

  /* Standard initialization */
  cs_array_real_fill_zero(n_i_faces, iflux);

  if (vp_param->staggered == 0) {
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      coefa_dp[f_id] = 0.;
      coefaf_dp[f_id] = 0.;
      coefb_dp[f_id] = coefb_p[f_id];
      coefbf_dp[f_id] = coefbf_p[f_id];
      bflux[f_id] = 0.;
    }
  }
  else {
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      coefa_dp[f_id] = coefa_p[f_id];
      coefaf_dp[f_id] = coefaf_p[f_id];
      coefb_dp[f_id] = coefb_p[f_id];
      coefbf_dp[f_id] = coefbf_p[f_id];
      bflux[f_id] = 0.;
    }
  }

  /* Compute a pseudo hydrostatic pressure increment stored
     in cvar_hydro_pres(.) with Homogeneous Neumann BCs everywhere. */

  int indhyd = 0;

  if (vp_param->iphydr == 1 && vp_param->icalhy == 1) {

    int ifcsor = isostd[n_b_faces];
    cs_parall_max(1, CS_INT_TYPE, &ifcsor);

    /* This computation is needed only if there are outlet faces */

    if (ifcsor > -1 || open_bcs_flag != 0)
      _hydrostatic_pressure_compute(m, fvq,
                                    &indhyd, iterns,
                                    frcxt, dfrcxt, cvar_hydro_pres,
                                    iflux, bflux, i_visc, b_visc,
                                    dam, xam,
                                    dphi, rhs);
  }

  /* Compute the BCs for the pressure increment
     (first we set the BCs of a standard pressure increment,
     that are (A = 0, B_dp = B_p) for the gradient BCs
     Then the A_dp is set thank to the pre-computed hydrostatic pressure
     so that the pressure increment will be 0 on the reference outlet face. */

  if (vp_param->iphydr == 1 || vp_param->iifren == 1) {

    cs_real_t phydr0 = 0.;

    if (f_hp != NULL && indhyd == 1) {

      cs_lnum_t f_id_0 = isostd[n_b_faces] - 1;
      if (f_id_0 > -1) {
        cs_lnum_t c_id_0 = b_face_cells[f_id_0];
        cs_real_t d[3] = {b_face_cog[f_id_0][0] - cell_cen[c_id_0][0],
                          b_face_cog[f_id_0][1] - cell_cen[c_id_0][1],
                          b_face_cog[f_id_0][2] - cell_cen[c_id_0][2]};
        phydr0 =   cvar_hydro_pres[c_id_0]
                 + d[0] * (dfrcxt[c_id_0][0] + frcxt[c_id_0][0])
                 + d[1] * (dfrcxt[c_id_0][1] + frcxt[c_id_0][1])
                 + d[2] * (dfrcxt[c_id_0][2] + frcxt[c_id_0][2]);
      }

      cs_parall_sum(1, CS_REAL_TYPE, &phydr0);  /* > 0 only on one rank */

      /* Rescale cvar_hydro_pres so that it is 0 on the reference face */
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cvar_hydro_pres[c_id] -= phydr0;

    }

    /* If hydrostatic pressure increment or free entrance Inlet. */

    if (indhyd == 1 || vp_param->iifren == 1) {

      const int *auto_flag = cs_glob_bc_pm_info->iautom;
      const cs_real_t *b_head_loss
        = cs_boundary_conditions_get_b_head_loss(false);

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

        int iautof = 0;

        /* automatic inlet/outlet face for atmospheric flow */
        if (cs_glob_atmo_option->meteo_profile > 0)
          iautof = auto_flag[f_id];

        if (isostd[f_id] == 1 || (open_bcs_flag >= 1 && iautof >= 1)) {

          cs_lnum_t c_id = b_face_cells[f_id];

          cs_real_t hint = 0;

          /* Diffusive flux BCs */
          if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
            hint = c_visc[c_id]/b_dist[f_id];

          /* Symmetric tensor diffusivity */
          else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) {

            cs_real_t visci[3][3], dist[3];
            const cs_real_t *n = b_face_normal[f_id];

            visci[0][0] = vitenp[c_id][0];
            visci[1][1] = vitenp[c_id][1];
            visci[2][2] = vitenp[c_id][2];
            visci[1][0] = vitenp[c_id][3];
            visci[0][1] = vitenp[c_id][3];
            visci[2][1] = vitenp[c_id][4];
            visci[1][2] = vitenp[c_id][4];
            visci[2][0] = vitenp[c_id][5];
            visci[0][2] = vitenp[c_id][5];

            dist[0] = b_face_cog[f_id][0] - cell_cen[c_id][0];
            dist[1] = b_face_cog[f_id][1] - cell_cen[c_id][1];
            dist[2] = b_face_cog[f_id][2] - cell_cen[c_id][2];

            /* ||Ki.S||^2 */
            cs_real_t viscis =   cs_math_pow2(  visci[0][0]*n[0]
                                              + visci[1][0]*n[1]
                                              + visci[2][0]*n[2])
                               + cs_math_pow2(  visci[0][1]*n[0]
                                              + visci[1][1]*n[1]
                                              + visci[2][1]*n[2])
                               + cs_math_pow2(  visci[0][2]*n[0]
                                              + visci[1][2]*n[1]
                                              + visci[2][2]*n[2]);

            /* IF.Ki.S */
            cs_real_t fikis
              = (  cs_math_3_dot_product(dist, visci[0]) * n[0]
                 + cs_math_3_dot_product(dist, visci[1]) * n[1]
                 + cs_math_3_dot_product(dist, visci[2]) * n[2]);

            /* Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
               NB: eps =1.d-1 must be consistent with vitens.f90 */

            fikis = fmax(fikis, 1.e-1*sqrt(viscis)*b_dist[f_id]);

            hint = viscis/b_face_surf[f_id]/fikis;

          }

          if (indhyd == 1) {
            if (f_hp != NULL) {
              coefa_dp[f_id] =   cvar_hydro_pres[c_id]
                               - cvar_hydro_pres_prev[c_id];
            }
            coefa_dp[f_id] += cs_math_3_distance_dot_product(cell_cen[c_id],
                                                             b_face_cog[f_id],
                                                             dfrcxt[c_id]);
          }

          /* Free entrance boundary face (Bernoulli condition to link the
             pressure increment and the predicted velocity) */

          if (bc_type[f_id] == CS_FREE_INLET) {

            /* Boundary mass flux of the predicted velocity */
            cs_real_t bpmasf = cs_math_3_dot_product(vel[c_id],
                                                     b_face_normal[f_id]);

            /* Ingoing mass Flux, Bernoulli relation ship is used */
            if (bpmasf <= 0) {

              /* Head loss of the fluid outside the domain, between infinity and
                 the entrance */

              cs_real_t kpdc = b_head_loss[f_id];
              cs_real_t rho = brom[f_id];
              cs_real_t cfl =   -(bmasfl[f_id]/b_face_surf[f_id]*dt[c_id])
                              / (2.*rho*b_dist[f_id])*(1. + kpdc);

              cs_real_t pimp = coefa_dp[f_id] - cvar_pr[c_id]
                               -   0.5 * (1. + kpdc) * bmasfl[f_id]*bpmasf
                                 / cs_math_pow2(b_face_surf[f_id]);

              cs_boundary_conditions_set_convective_outlet_scalar
                (&(coefa_dp[f_id]), &(coefaf_dp[f_id]),
                 &(coefb_dp[f_id]), &(coefbf_dp[f_id]),
                 pimp, cfl, hint);

            }

            else
              coefaf_dp[f_id] = - hint*coefa_dp[f_id];

          }

          /* Outher boundary face types */

          else
            coefaf_dp[f_id] = - hint*coefa_dp[f_id];

        } /* if (isostd[f_id] == 1 || (open_bcs_flag >= 1 && iautof >= 1)) */

      } /* End of loop on boundary faces */

    }  /* if (indhyd == 1 || vp_param->iifren == 1) */

  } /* if (vp_param->iphydr == 1 && vp_param->iifren == 1) */

  /* Building the linear system to solve.
   * ==================================== */

  /* Implicit term */

  cs_real_t *rovsdt;
  BFT_MALLOC(rovsdt, n_cells_ext, cs_real_t);

  cs_array_real_fill_zero(n_cells, rovsdt);

  /* Compressible scheme implicit part;
     Getting the thermal parameters */

  int ieos = cs_glob_cf_model->ieos;
  int thermal_variable = cs_glob_thermal_model->thermal_variable;
  int kinetic_st = cs_glob_thermal_model->has_kinetic_st;

  const cs_real_t *temp = NULL;
  cs_real_t *xcpp = NULL, *dc2 = NULL;
  cs_real_t *cvar_th = NULL, *tempk = NULL;
  cs_real_t _coef = 0;
  cs_real_t *yw = NULL, *yv = NULL;

  if (idilat == 2) {

    /* Get the temperature */
    if (thermal_variable == CS_THERMAL_MODEL_TEMPERATURE)
      temp = CS_F_(t)->val;
    else {
      const cs_field_t *f_t = cs_field_by_name_try("temperature");
      if (f_t != NULL) {
        tempk = f_t->val;
        temp = f_t->val;
      }
    }

    if (temp != NULL) {

      cvar_th = CS_F_(t)->val;

      /* Allocation */
      BFT_MALLOC(dc2, n_cells_ext, cs_real_t);
      BFT_MALLOC(xcpp, n_cells_ext, cs_real_t);

      /* Theta scheme related term */
      _coef = 1. + 2. * (1. - eqp_u->thetav);

      /* Get cp */
      if (fluid_props->icp > 0) {
        cs_real_t *cpro_cp = CS_F_(cp)->val;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          xcpp[c_id] = cpro_cp[c_id];
      }
      else
        cs_array_set_value_real(n_cells, 1, 1., xcpp);

      /* Get mass fractions if needed */
      if (ieos == CS_EOS_MOIST_AIR) {
        yw = cs_field_by_name("yw")->val;
        yv = cs_field_by_name("yv")->val;
      }
      cs_real_t *cvar_fracm = NULL;

      /* Compute dc2 */
      cs_thermal_model_c_square(xcpp,
                                temp,
                                cvar_pr,
                                yv,
                                cvar_fracm,
                                yw,
                                dc2);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        rovsdt[c_id] += cell_f_vol[c_id] * _coef * dc2[c_id] / dt[c_id];

      BFT_FREE(xcpp);

    }
  }

  /* Implicit part of the cavitation source */
  if (i_vof_mass_transfer != 0 && cavitation_parameters->itscvi == 1) {
    cs_real_t *dgdpca = cs_get_cavitation_dgdp_st();
    cs_real_t rho1 = vof_parameters->rho1;
    cs_real_t rho2 = vof_parameters->rho2;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rovsdt[c_id] -= cell_f_vol[c_id] * dgdpca[c_id] * (1./rho2 - 1./rho1);
  }

  /* Strengthen the diagonal for Low Mach Algorithm */
  if (idilat == 3) {
    const cs_real_t epsdp = vp_param->epsdp;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rovsdt[c_id] += epsdp*cell_f_vol[c_id]/dt[c_id];
  }

  cs_real_t *c2 = NULL;

  if (compressible_flag == 3) {

    cs_real_t *cvar_fracv = NULL;
    cs_real_t *cvar_fracm = NULL;
    cs_real_t *cvar_frace = NULL;

    cs_real_t *cpro_cp = NULL;
    cs_real_t *cpro_cv = NULL;

    if (fluid_props->icp >= 0)
      cpro_cp = cs_field_by_id(fluid_props->icp)->val;

    if (fluid_props->icv >= 0)
      cpro_cv = cs_field_by_id(fluid_props->icv)->val;

    BFT_MALLOC(c2, n_cells_ext, cs_real_t);

    cs_cf_thermo_c_square(cpro_cp, cpro_cv, cvar_pr, crom,
                          cvar_fracv, cvar_fracm, cvar_frace, c2, n_cells);

    int istat = eqp_p->istat;
    if (istat != 0) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        rovsdt[c_id] += istat*(cell_f_vol[c_id]/(dt[c_id]*c2[c_id]));
    }

  }

  /* Face diffusivity */

  cs_real_t *weighb = NULL;
  cs_real_2_t *weighf = NULL;

  if (eqp_p->idiff >= 1) {

    /* Scalar diffusivity */
    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION) {

      cs_face_viscosity(m,
                        fvq,
                        eqp_p->imvisf,
                        c_visc,
                        i_visc,
                        b_visc);

      if (f_weight != NULL) {  /* Weighting for gradient */
        cs_array_real_copy(n_cells, c_visc, f_weight->val);
        cs_halo_sync_var(m->halo, CS_HALO_STANDARD, f_weight->val);
      }

    }

    /* Tensor diffusivity */
    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) {

      /* Allocate temporary arrays */
      BFT_MALLOC(weighf, n_i_faces, cs_real_2_t);
      BFT_MALLOC(weighb, n_b_faces, cs_real_t);

      cs_face_anisotropic_viscosity_scalar(m,
                                           fvq,
                                           vitenp,
                                           eqp_p->verbosity,
                                           weighf,
                                           weighb,
                                           i_visc,
                                           b_visc);

      if (f_weight != NULL) { /* Weighting for gradient */
        cs_array_real_copy(6*n_cells, (const cs_real_t *)vitenp, f_weight->val);
        cs_mesh_sync_var_sym_tens((cs_real_6_t *)(f_weight->val));
      }

    }

  }
  else {

    cs_array_real_fill_zero(n_i_faces, i_visc);
    cs_array_real_fill_zero(n_b_faces, b_visc);

  }

  if (vp_param->staggered == 1) {
    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id_0 = i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = i_face_cells[f_id][1];
      cs_real_t dtm = 0.5 * (dt[c_id_0] + dt[c_id_1]);
      i_visc[f_id] = taui[f_id] / dtm * i_visc[f_id];
    }
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      cs_lnum_t c_id = b_face_cells[f_id];
      b_visc[f_id] = taub[f_id] / dt[c_id] * b_visc[f_id];
    }
  }

  cs_matrix_wrapper_scalar(eqp_p->iconv, eqp_p->idiff, eqp_p->ndircl,
                           isym,
                           1.0, /* thetap */
                           0,   /* imucpp */
                           coefb_dp, coefbf_dp, rovsdt,
                           imasfl, bmasfl, i_visc, b_visc,
                           NULL, dam, xam);

  /* Mass flux initialization
     ------------------------ */

  /* Predicted mass flux and first Rhie and Chow component */

  cs_real_3_t  *gradp;
  BFT_MALLOC(gradp, n_cells_ext, cs_real_3_t);

  cs_gradient_porosity_balance(1);

  /* Pressure gradient.
     NB: for the VOF algo. the weighting is automatically done
     through the iwgrec variable calculation option.
  */

  cs_field_gradient_potential(f_p,
                              0,   /* iprev */
                              1,   /* inc, not by increment */
                              vp_param->iphydr,
                              frcxt,
                              gradp);

  const cs_real_t arak = vp_param->arak;

  /* Rhie and Chow filter */
  if (arak > 0.) {

    if (vp_param->iphydr == 1) {

      if (cs_glob_porous_model == 3) {

        cs_real_3_t *cpro_poro_div_duq
          = (cs_real_3_t *)(cs_field_by_name("poro_div_duq")->val);

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          for (cs_lnum_t j = 0; j < 3; j++) {
            wrk[c_id][j] =   gradp[c_id][j] - frcxt[c_id][j]
                           - cpro_poro_div_duq[c_id][j];
          }
      }

      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          for (cs_lnum_t j = 0; j < 3; j++) {
            wrk[c_id][j] = gradp[c_id][j] - frcxt[c_id][j];
          }
      }

    }
    else {

      cs_array_real_copy(3*n_cells, (const cs_real_t *)gradp, (cs_real_t *)wrk);

    }

    if (   (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
        && vp_param->rcfact == 0) {

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t arsrdt = (arak/crom[c_id]) * dt[c_id];
        for (cs_lnum_t j = 0; j < 3; j++)
          wrk[c_id][j] *= arsrdt;
      }

    }
    else if (   (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
             ||  vp_param->rcfact == 1) {

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t arsr = arak/crom[c_id];
        cs_real_t wrkn[3];

        wrkn[0] = arsr * (  da_uu[c_id][0]*wrk[c_id][0]
                          + da_uu[c_id][3]*wrk[c_id][1]
                          + da_uu[c_id][5]*wrk[c_id][2]);
        wrkn[1] = arsr * (  da_uu[c_id][3]*wrk[c_id][0]
                          + da_uu[c_id][1]*wrk[c_id][1]
                          + da_uu[c_id][4]*wrk[c_id][2]);
        wrkn[2] = arsr * (  da_uu[c_id][5]*wrk[c_id][0]
                          + da_uu[c_id][4]*wrk[c_id][1]
                          + da_uu[c_id][2]*wrk[c_id][2]);

        wrk[c_id][0] = wrkn[0];
        wrk[c_id][1] = wrkn[1];
        wrk[c_id][2] = wrkn[2];
      }

    }

  }

  else {
    cs_array_real_fill_zero(3*n_cells, (cs_real_t *)wrk);
  }

  if (idilat < 4) {
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        wrk[c_id][j] += vel[c_id][j];
        wrk2[c_id][j] = wrk[c_id][j];
      }
    }
  }

  /* Sync for parallelism and periodicity */
  cs_mesh_sync_var_vect((cs_real_t *)wrk);

  {
    /* BCs will be taken into account later if idilat >= 4 */
    int inc = (idilat >= 4) ? 0 : 1;
    int iflmb0 = (cs_glob_ale > CS_ALE_NONE) ? 0 : 1;
    int itypfl = (vof_parameters->vof_model > 0 || idilat == 4) ? 0 : 1;

    cs_mass_flux(m,
                 fvq,
                 -1, /* f_id */
                 itypfl,
                 iflmb0,
                 1, /* init */
                 inc,
                 eqp_u->imrgra,
                 eqp_u->nswrgr,
                 eqp_u->imligr,
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 wrk,
                 coefav, coefbv,
                 imasfl, bmasfl);
  }

  if (vp_param->staggered == 1) {

    const cs_time_step_t *ts = cs_glob_time_step;

    cs_real_t *imasfla
      = cs_field_by_id(cs_field_get_key_int(f_p, kimasf))->val_pre;
    cs_real_t *bmasfla
      = cs_field_by_id(cs_field_get_key_int(f_p, kbmasf))->val_pre;

    cs_real_t *sti = cs_field_by_name("inner_face_source_term")->val;

    if (ts->nt_cur == 1 && ts->nt_prev == 0) {
      cs_array_real_copy(n_i_faces, imasfl, imasfla);
      cs_array_real_copy(n_b_faces, bmasfl, bmasfla);
    }

    if (cs_glob_porous_model >= 1) {

      const cs_real_t *c_porosity = cs_field_by_name("porosity")->val;

      for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
        cs_lnum_t c_id_0 = i_face_cells[f_id][0];
        cs_lnum_t c_id_1 = i_face_cells[f_id][1];
        cs_real_t dtm = 0.5 * (dt[c_id_0]+dt[c_id_1]);
        cs_real_t porosf = fmin(c_porosity[c_id_0], c_porosity[c_id_1]);
        imasfl[f_id] =   taui[f_id] / dtm
                       * imasfla[f_id]+porosf*taui[f_id]*sti[f_id]
                       * i_f_face_surf[f_id];
      }

    }

    else { /* cs_glob_porous_model == 0) */

      for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
        cs_lnum_t c_id_0 = i_face_cells[f_id][0];
        cs_lnum_t c_id_1 = i_face_cells[f_id][1];
        cs_real_t dtm = 0.5 * (dt[c_id_0] + dt[c_id_1]);
        imasfl[f_id] =   taui[f_id] / dtm
                       * imasfla[f_id]+taui[f_id]*sti[f_id]
                       * i_f_face_surf[f_id];
      }

    } /* end of test on cs_glob_porous_model */

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      cs_lnum_t c_id = b_face_cells[f_id];
      if (bc_type[f_id] == CS_INLET)
        bmasfl[f_id] = taub[f_id] / dt[c_id] * bmasfl[f_id];
      else
        bmasfl[f_id] = taub[f_id] / dt[c_id] * bmasfla[f_id];
    }

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      if (bc_type[f_id] == CS_INLET) {
        cs_lnum_t c_id = b_face_cells[f_id];
        cs_real_t dimp =   -(1. - dt[c_id]/taub[f_id])
                         * bmasfl[f_id]/b_f_face_surf[f_id];
        cs_real_t hint = taub[f_id] / b_dist[f_id];

        cs_boundary_conditions_set_neumann_scalar
          (&(coefa_dp[f_id]), &(coefaf_dp[f_id]),
           &(coefb_dp[f_id]), &(coefbf_dp[f_id]),
           dimp, hint);
      }
    }

  } /* end if (vp_param->staggered == 1) */

  /* Project exterior forces to faces */

  if (vp_param->iphydr == 1) {

    cs_gradient_porosity_balance(0);

    /* Scalar diffusivity */
    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
      cs_ext_force_flux(m,
                        fvq,
                        0,  /* init */
                        eqp_p->nswrgr,
                        dfrcxt,
                        coefbf_dp,
                        imasfl,
                        bmasfl,
                        i_visc,
                        b_visc,
                        c_visc,
                        c_visc,
                        c_visc);

    /* Tensor diffusivity */
    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
      cs_ext_force_anisotropic_flux(m,
                                    fvq,
                                    0,  /* init */
                                    eqp_p->nswrgr,
                                    eqp_p->ircflu,
                                    dfrcxt,
                                    coefbf_dp,
                                    i_visc,
                                    b_visc,
                                    vitenp,
                                    weighf,
                                    imasfl,
                                    bmasfl);

  }

  cs_gradient_porosity_balance(1);

  /* Rhie and Chow filter
     ==================== */

  if (arak > 0.) {

    cs_real_t *ipro_visc = NULL, *bpro_visc = NULL;

    BFT_MALLOC(ipro_visc, n_i_faces, cs_real_t);
    BFT_MALLOC(bpro_visc, n_b_faces, cs_real_t);

    /* Scalar diffusivity */
    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION && vp_param->rcfact == 0) {

      cs_real_t *cpro_visc;
      BFT_MALLOC(cpro_visc, n_cells_ext, cs_real_t);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cpro_visc[c_id] = arak * c_visc[c_id];

      cs_face_viscosity(m,
                        fvq,
                        eqp_p->imvisf,
                        cpro_visc,
                        ipro_visc,
                        bpro_visc);

      /* We cancel the face viscosity for coupled faces so as not to modify the
         boundary mass flux in the case of a pressure Dirichlet:
         pressure correction and filter are canceled. */

      if (cs_sat_coupling_n_couplings() > 0) {
        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
          if (bc_type[f_id] == CS_COUPLED_FD)
            bpro_visc[f_id] = 0.;
        }
      }

      cs_face_diffusion_potential(-1,
                                  m,
                                  fvq,
                                  0,  /* init */
                                  1,  /* inc */
                                  eqp_p->imrgra,
                                  eqp_p->nswrgr,
                                  eqp_p->imligr,
                                  vp_param->iphydr,
                                  eqp_p->iwgrec,
                                  eqp_p->verbosity,
                                  eqp_p->epsrgr,
                                  eqp_p->climgr,
                                  frcxt,
                                  cvar_pr,
                                  coefa_p, coefb_p,
                                  coefaf_p, coefbf_p,
                                  ipro_visc, bpro_visc, cpro_visc,
                                  imasfl, bmasfl);

      /* Project source term to remove hydrostatic part from pressure */

      if (vp_param->iphydr == 1)
        cs_ext_force_flux(m,
                          fvq,
                          0,  /* init */
                          eqp_p->nswrgr,
                          frcxt,
                          coefbf_p,
                          imasfl, bmasfl,
                          ipro_visc, bpro_visc,
                          cpro_visc, cpro_visc, cpro_visc);

      BFT_FREE(cpro_visc);

    }

    /* Tensor diffusivity */
    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION || vp_param->rcfact == 1) {

      cs_real_6_t *cpro_vitenp;
      BFT_MALLOC(cpro_vitenp, n_cells_ext, cs_real_6_t);

      if (idilat == 4 || vof_parameters->vof_model > 0) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t arsr = arak / crom[c_id];
          for (cs_lnum_t j = 0; j < 6; j++)
            cpro_vitenp[c_id][j] = arsr * da_uu[c_id][j];
        }
      }
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          for (cs_lnum_t j = 0; j < 6; j++)
            cpro_vitenp[c_id][j] = arak * da_uu[c_id][j];
        }
      }

      cs_real_2_t *weighftp = NULL;
      cs_real_t *weighbtp = NULL;
      BFT_MALLOC(weighftp, n_i_faces, cs_real_2_t);
      BFT_MALLOC(weighbtp, n_b_faces, cs_real_t);

      /* A harmonic mean is used regardless of the imvisf option. */

      cs_face_anisotropic_viscosity_scalar(m,
                                           fvq,
                                           cpro_vitenp,
                                           eqp_p->verbosity,
                                           weighftp, weighbtp,
                                           ipro_visc, bpro_visc);

      /* We cancel the face viscosity for coupled faces so as not to modify the
         boundary mass flux in the case of a pressure Dirichlet:
         pressure correction and filter are canceled. */

      if (cs_sat_coupling_n_couplings() > 0) {
        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
          if (bc_type[f_id] == CS_COUPLED_FD)
            bpro_visc[f_id] = 0.;
        }
      }

      cs_face_anisotropic_diffusion_potential(-1,
                                              m,
                                              fvq,
                                              0,  /* init */
                                              1,  /* inc */
                                              eqp_p->imrgra,
                                              eqp_p->nswrgr,
                                              eqp_p->imligr,
                                              eqp_p->ircflu,
                                              vp_param->iphydr,
                                              eqp_p->iwgrec,
                                              eqp_p->verbosity,
                                              eqp_p->epsrgr,
                                              eqp_p->climgr,
                                              frcxt,
                                              cvar_pr,
                                              coefa_p, coefb_p,
                                              coefaf_p, coefbf_p,
                                              ipro_visc, bpro_visc, cpro_vitenp,
                                              weighftp, weighbtp,
                                              imasfl, bmasfl);

      /* Project source term to remove hydrostatic part from pressure */

      if (vp_param->iphydr == 1)
        cs_ext_force_anisotropic_flux(m,
                                      fvq,
                                      0,  /* init */
                                      eqp_p->nswrgr,
                                      eqp_p->ircflu,
                                      frcxt,
                                      coefbf_p,
                                      ipro_visc, bpro_visc,
                                      cpro_vitenp,
                                      weighftp,
                                      imasfl, bmasfl);

      BFT_FREE(cpro_vitenp);
      BFT_FREE(weighftp);
      BFT_FREE(weighbtp);
    }

    BFT_FREE(ipro_visc);
    BFT_FREE(bpro_visc);
  }

  /*
   * Solving (Loop over the non-orthogonalities)
   * =========================================== */

  int nswmpr = eqp_p->nswrsm;   /* Number of sweeps */

  /* Variables are set to 0 (or hydro pressure increment)
   *   phi        is the increment of the pressure
   *   dphi       is the increment of the increment between sweeps
   *   cpro_divu  is the initial divergence of the predicted mass flux */

  if (vp_param->staggered == 0) {
    if (f_hp != NULL) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        phi[c_id]  = cvar_hydro_pres[c_id] - cvar_hydro_pres_prev[c_id];
        dphi[c_id] = 0.;
        phia[c_id] = phi[c_id];
      }
    }
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        phi[c_id] = 0.;
        dphi[c_id] = 0.;
        phia[c_id] = phi[c_id];
      }
    }
  }
  else {
    if (f_hp != NULL) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        phi[c_id]  = cvara_pr[c_id] + cvar_hydro_pres[c_id]
                                    - cvar_hydro_pres_prev[c_id];
        dphi[c_id] = 0.;
        phia[c_id] = phi[c_id];
      }
    }
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        phi[c_id]  = cvara_pr[c_id];
        dphi[c_id] = 0.;
        phia[c_id] = phi[c_id];
      }
    }
  }

  /* Initial divergence */

  cs_divergence(m,
                1,  /* init */
                imasfl,
                bmasfl,
                cpro_divu);

  /* Weakly compressible algorithm: semi analytic scheme
   *   1. The RHS contains rho div(u*) and not div(rho u*)
   *   2. Add dilatation source term to rhs
   *   3. The mass flux is completed by rho u* . S
   */

  cs_real_t *velflx = NULL, *velflb = NULL;

  if (idilat >= 4) {

    cs_real_t *cpro_tsrho = cs_field_by_name("dila_st")->val;

    BFT_MALLOC(velflx, n_i_faces, cs_real_t);
    BFT_MALLOC(velflb, n_b_faces, cs_real_t);

    /* 1. The RHS contains rho div(u*) and not div(rho u*) */

    int init = 1;
    int iflmb0 = (cs_glob_ale > CS_ALE_NONE) ? 0 : 1;

    cs_mass_flux(m,
                 fvq,
                 f_vel->id,
                 0,  /* itypfl */
                 iflmb0,
                 init,
                 1,  /* inc */
                 eqp_u->imrgra,
                 eqp_u->nswrgr,
                 eqp_u->imligr,
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 vel,
                 coefav, coefbv,
                 velflx, velflb);

    cs_divergence(m, init, velflx, velflb, res);

    if (idilat == 4) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cpro_divu[c_id] += res[c_id];
    }
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cpro_divu[c_id] += res[c_id]*crom[c_id];
    }

    /* 2. Add the dilatation source term D(rho)/Dt */

    if (idilat == 4) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cpro_divu[c_id] += cpro_tsrho[c_id] / crom[c_id];
    }
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cpro_divu[c_id] += cpro_tsrho[c_id];
    }

    /* 3. The mass flux is completed by u*.S (idilat=4)
     *                                  rho u* . S (idilat=5)
     */

    int itypfl = (idilat == 4) ? 0 : 1;

    cs_mass_flux(m,
                 fvq,
                 f_vel->id,
                 itypfl,
                 iflmb0,
                 0,  /* init */
                 1,  /* inc */
                 eqp_u->imrgra,
                 eqp_u->nswrgr,
                 eqp_u->imligr,
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 vel,
                 coefav, coefbv,
                 imasfl, bmasfl);

  }

  /* Mass source terms adding for volumic flow rate */

  cs_lnum_t ncesmp = 0;
  cs_lnum_t *icetsm = NULL;
  int *itpsmp = NULL;
  cs_real_t *smcelp, *gamma = NULL;

  cs_volume_mass_injection_get_arrays(f_p, &ncesmp, &icetsm, &itpsmp,
                                      &smcelp, &gamma);

  if (ncesmp > 0) {
    for (cs_lnum_t c_idx = 0; c_idx < ncesmp; c_idx++) {
      cs_lnum_t c_id = icetsm[c_idx] - 1;
      cpro_divu[c_id] -= cell_f_vol[c_id] * smcelp[c_idx];
    }
  }

  /* Source term adding for condensation modelling */

  if (nfbpcd > 0) {
    const int var_id_key = cs_field_key_id("variable_id");
    const int ipr = cs_field_get_key_int(f_p, var_id_key);

    cs_real_t *_spcond = spcond + (ipr-1)*nfbpcd;

    for (cs_lnum_t f_idx = 0; f_idx < nfbpcd; f_idx++) {
      cs_lnum_t f_id = ifbpcd[f_idx];
      cs_lnum_t c_id = b_face_cells[f_id];
      cpro_divu[c_id] -= b_face_surf[f_id] * _spcond[f_idx];
    }
  }

  /* volume Gamma source for metal mass structures
     condensation modelling */

  if (ncmast > 0) {
    const int var_id_key = cs_field_key_id("variable_id");
    const int ipr = cs_field_get_key_int(f_p, var_id_key);

    cs_real_t *_svcond = svcond + (ipr-1)*ncmast;
    cs_real_t *surfbm = NULL;
    BFT_MALLOC(surfbm, ncmast, cs_real_t);

    cs_wall_condensation_volume_exchange_surf_at_cells(surfbm);

    for (cs_lnum_t c_idx = 0; c_idx < ncmast; c_idx++) {
      cs_lnum_t c_id = ltmast[c_idx];
      cpro_divu[c_id] -= surfbm[c_idx]* _svcond[c_idx];
    }

    BFT_FREE(surfbm);
  }

  /* Source term associated to the mass aggregation */

  if ((idilat == 2 || idilat == 3) && compressible_flag != 3) {
    if (ieos == CS_EOS_NONE) { // If no particular EOS is set
      if (vp_param->itpcol == 1 && eqp_u->thetav < 1.) {
        cs_real_t *imasfla
          = cs_field_by_id(cs_field_get_key_int(f_p, kimasf))->val_pre;
        cs_real_t *bmasfla
          = cs_field_by_id(cs_field_get_key_int(f_p, kbmasf))->val_pre;
        cs_real_t *divu_prev;
        BFT_MALLOC(divu_prev, n_cells_ext, cs_real_t);

        cs_divergence(m, 1, imasfla, bmasfla, divu_prev);

        cs_real_t thetav = eqp_u->thetav;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t drom = crom_eos[c_id] - croma[c_id];
          cpro_divu[c_id] +=  (1. + thetav) *drom *cell_f_vol[c_id] /dt[c_id]
            + thetav *divu_prev[c_id];
        }
        BFT_FREE(divu_prev);
      }
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t drom = crom_eos[c_id] - croma[c_id];
          cpro_divu[c_id] += drom *cell_f_vol[c_id] /dt[c_id];
        }
      }
    } else {/* compressible scheme explicit part*/
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t drom = crom_eos[c_id] - croma[c_id];
        cs_real_t drop = (-1 + _coef ) * (cvar_pr[c_id] - cvara_pr[c_id])
          * dc2[c_id];
        cpro_divu[c_id] += (drom + drop) *cell_f_vol[c_id] /dt[c_id];
      }
    }

  }

  /* Lagrangian source terms */

  if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
      && cs_glob_lagr_source_terms->ltsmas == 1) {

    const cs_lagr_source_terms_t  *lag_st = cs_glob_lagr_source_terms;

    cs_lnum_t itsmas = cs_glob_lagr_source_terms->itsmas;
    cs_real_t *lag_st_m = lag_st->st_val + (itsmas-1)*n_cells_ext;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      cpro_divu[c_id] -= lag_st_m[c_id];
  }

  /* Cavitation source term */

  if (i_vof_mass_transfer != 0) {
    cs_real_t *gamcav = cs_get_cavitation_gam();
    cs_real_t rho1 = vof_parameters->rho1;
    cs_real_t rho2 = vof_parameters->rho2;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      cpro_divu[c_id] -= cell_f_vol[c_id]*gamcav[c_id]*(1./rho2 - 1./rho1);
  }

  /* Norm residual
   * Historical norm for the pressure step:
   *   div(rho u* + dt gradP^(n))-Gamma
   *   i.e.  RHS of the pressure + div(dt gradP^n) (otherwise there is a risk
   *   a 0 norm at steady states...). Represents terms that pressure has to
   *   balance. */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t dt_d_rho = dt[c_id] / crom[c_id];
    for (cs_lnum_t j = 0; j < 3; j++) {
      wrk[c_id][j] = dt_d_rho * gradp[c_id][j];
    }
  }

  /* Parallelism and periodicity */
  cs_mesh_sync_var_vect((cs_real_t *)wrk);

  {
    int iflmb0 = (cs_glob_ale > CS_ALE_NONE) ? 0 : 1;

    /* VOF algorithm: the pressure step corresponds to the
       correction of the volumetric flux, not the mass flux */
    int itypfl = (vof_parameters->vof_model > 0 || idilat >= 4) ? 0 : 1;

    cs_mass_flux(m,
                 fvq,
                 -1,
                 itypfl,
                 iflmb0,
                 1,  /* init */
                 1,  /* inc */
                 eqp_p->imrgra,
                 1,  /* nswrgu (to save time, no space reconstruction */
                 eqp_u->imligr,
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 wrk,
                 coefav, coefbv,
                 iflux, bflux);

    cs_divergence(m, 1, iflux, bflux, res);
  }

  BFT_FREE(iflux);
  BFT_FREE(bflux);

  /* Weakly compressible algorithm: semi analytic scheme */
  if (idilat >= 4) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      res[c_id] *= crom[c_id];
  }

  /* It is: div(dt/rho*rho grad P) + div(rho u*) - Gamma
     NB: if iphydr=1, div(rho u*) contains div(dt d fext).  */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    res[c_id] += cpro_divu[c_id];
  }

  /* Pressure norm */
  cs_real_t rnormp = sqrt(cs_gdot(n_cells, res, res));

  if (eqp_p->iwarni >= 2)
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" %-16s : normed residuals = %14.6e\n"),
                  f_p->name, rnormp);

  if (iterns <= 1) {
    sinfo->n_it = 0;
  }

  /* Dynamic relaxation initialization
     --------------------------------- */

  cs_real_t nadxk = -1;
  cs_real_t rnorm2 = -1;

  if (eqp_p->iswdyn >= 1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
      adxkm1[c_id] = 0.;
      adxk[c_id] = 0.;
    }

    /* ||A.dx^0||^2 = 0 */
    nadxk = 0.;
    rnorm2 = cs_math_pow2(rnormp);
  }

  /* Initial right hand side */

  {
    int inc = 0;  /* by increment */
    if (   vp_param->iphydr == 1
        || vp_param->iifren == 1
        || vp_param->staggered == 1)
      inc = 1;    /* not by increment */

    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
      cs_diffusion_potential(-1,
                             m,
                             fvq,
                             1,  /* init */
                             inc,
                             eqp_p->imrgra,
                             eqp_p->nswrgr,
                             eqp_p->imligr,
                             vp_param->iphydr,
                             eqp_p->iwgrec,
                             eqp_p->verbosity,
                             eqp_p->epsrgr,
                             eqp_p->climgr,
                             dfrcxt,
                             phi,
                             coefa_dp, coefb_dp,
                             coefaf_dp, coefbf_dp,
                             i_visc, b_visc,
                             c_visc,
                             rhs);

    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
      cs_anisotropic_diffusion_potential(-1,
                                         m,
                                         fvq,
                                         1,  /* init */
                                         inc,
                                         eqp_p->imrgra,
                                         eqp_p->nswrgr,
                                         eqp_p->imligr,
                                         eqp_p->ircflu,
                                         vp_param->iphydr,
                                         eqp_p->iwgrec,
                                         eqp_p->verbosity,
                                         eqp_p->epsrgr,
                                         eqp_p->climgr,
                                         dfrcxt,
                                         phi,
                                         coefa_dp, coefb_dp,
                                         coefaf_dp, coefbf_dp,
                                         i_visc, b_visc,
                                         vitenp,
                                         weighf, weighb,
                                         rhs);
  }

  /* Dynamic relaxation: stores the initial rhs */

  if (eqp_p->iswdyn >= 1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rhs0[c_id] = rhs[c_id];
  }

  /* Finalize the rhs initialization */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    rhs[c_id] = -rhs[c_id] - cpro_divu[c_id] - rovsdt[c_id]*phi[c_id];
  }

  /* Right hand side residual */
  cs_real_t residual = sqrt(cs_gdot(n_cells, rhs, rhs));

  sinfo->rhs_norm = residual;

  /* Pressure variation for the log */
  if (rnormp < cs_math_epzero)
    sinfo->derive = - sinfo->rhs_norm;
  else
    sinfo->derive = sinfo->rhs_norm/rnormp;

  int isweep = 1;

  const char fmt_sweep_residual_relax[]
    = N_(" %-16s : sweep = %d, right hand side norm = %14.6e, relaxp = %g\n");
  const char fmt_sweep_residual_info[]
    = N_(" %-16s : Current reconstruction sweep = %d\n"
         "         sweep residual = %12.5e norm = %12.5e\n"
         "         number of sweeps for solver = %d\n");

  /* Writing */
  if (eqp_p->verbosity >= 2) {
    cs_log_printf(CS_LOG_DEFAULT,
                  _(fmt_sweep_residual_relax),
                  f_p->name, isweep, residual, eqp_p->relaxv);
  }

  /* Reconstruction loop (beginning)
     ------------------------------- */

  int niterf = 0;

  while (isweep <= nswmpr && residual > eqp_p->epsrsm*rnormp) {

    /* Solving on the increment dphi
        ----------------------------- */

    if (eqp_p->iswdyn >= 1) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        dphim1[c_id] = dphi[c_id];
        dphi[c_id] = 0.;
      }
    }
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        dphi[c_id] = 0.;
    }

    cs_real_t ressol = residual;   /* solver residual */

    cs_sles_solve_native(f_p->id, NULL,
                         symmetric, 1, 1,
                         dam, xam,
                         eqp_p->epsilo,
                         rnormp,
                         &niterf,
                         &ressol,
                         rhs,
                         dphi);


    /* Dynamic relaxation of the system
       -------------------------------- */

    cs_real_t alph = -1, beta = -1;

    if (eqp_p->iswdyn >= 1) {

      // Computation of the variable relaxation coefficient

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
        adxkm1[c_id] = adxk[c_id];
        adxk[c_id] = - rhs0[c_id];
      }

      int init = 0;
      int inc = 0;  /* by increment */
      if (   vp_param->iphydr == 1
          || vp_param->iifren == 1)
        inc = 1;    /* not by increment */

      if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
        cs_diffusion_potential(-1,
                               m,
                               fvq,
                               init,
                               inc,
                               eqp_p->imrgra,
                               eqp_p->nswrgr,
                               eqp_p->imligr,
                               vp_param->iphydr,
                               eqp_p->iwgrec,
                               eqp_p->verbosity,
                               eqp_p->epsrgr,
                               eqp_p->climgr,
                               dfrcxt,
                               dphi,
                               coefa_dp, coefb_dp,
                               coefaf_dp, coefbf_dp,
                               i_visc, b_visc,
                               c_visc,
                               adxk);

      else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
        cs_anisotropic_diffusion_potential(-1,
                                           m,
                                           fvq,
                                           init,
                                           inc,
                                           eqp_p->imrgra,
                                           eqp_p->nswrgr,
                                           eqp_p->imligr,
                                           eqp_p->ircflu,
                                           vp_param->iphydr,
                                           eqp_p->iwgrec,
                                           eqp_p->verbosity,
                                           eqp_p->epsrgr,
                                           eqp_p->climgr,
                                           dfrcxt,
                                           dphi,
                                           coefa_dp, coefb_dp,
                                           coefaf_dp, coefbf_dp,
                                           i_visc, b_visc,
                                           vitenp,
                                           weighf, weighb,
                                           adxk);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        adxk[c_id] = - adxk[c_id];
      }

      // ||E.dx^(k-1)-E.0||^2
      cs_real_t nadxkm1 = nadxk;

      // ||E.dx^k-E.0||^2
      nadxk = cs_gdot(n_cells, adxk, adxk);

      // < E.dx^k-E.0; r^k >
      cs_real_t paxkrk = cs_gdot(n_cells, rhs, adxk);

      // Relaxation with respect to dx^k and dx^(k-1)

      cs_real_t paxm1ax = 1., paxm1rk = 0.;
      beta = 0;
      paxm1ax = 0.;

      if (eqp_p->iswdyn >= 2) {

        // < E.dx^(k-1)-E.0; r^k >
        paxm1rk = cs_gdot(n_cells, rhs, adxkm1);

        // < E.dx^(k-1)-E.0; E.dx^k -E.0 >
        paxm1ax = cs_gdot(n_cells, adxk, adxkm1);
        cs_real_t paxm1ax2 = cs_math_pow2(paxm1ax);

        if (   nadxkm1 > 1.e-30*rnorm2
            && (nadxk*nadxkm1-paxm1ax2) > 1.e-30*rnorm2)
          beta = (paxkrk*paxm1ax - nadxk*paxm1rk)/(nadxk*nadxkm1-paxm1ax2);
        else
          beta = 0.;

      }

      // The first sweep is not relaxed
      if (isweep == 1) {
        alph = 1.;
        beta = 0.;
      }
      else if (isweep == 2) {
        beta = 0.;
        if (fmax(nadxk, 1.e-30*rnorm2) <= 1.e-30)
          alph = 1.;
        else
          alph = -paxkrk/fmax(nadxk, 1.e-30*rnorm2);
      }
      else {
        alph = -(paxkrk + beta*paxm1ax)/fmax(nadxk, 1.e-30*rnorm2);
      }

      /* Writing */
      if (eqp_p->verbosity >= 2) {
        cs_log_printf
          (CS_LOG_DEFAULT,
           _(" %-16s : sweep = %d, Dynamic relaxation: alpha = %12.5e "
             "beta = %12.5e\n"
             "    < dI^k  ; R^k > = %12.5e ||dI^k  ||^2 = %12.5e\n"
             "    < dI^k-1; R^k > = %12.5e ||dI^k-1||^2 = %12.5e\n"
             "    < dI^k-1; dI^k > = %12.5e\n"),
           f_p->name, isweep, alph, beta,
           paxkrk, nadxk, paxm1rk, nadxkm1, paxm1ax);
      }

    } /* End if dynamic relaxation ((eqp_p->iswdyn >= 1)) */

    /* Update the  pressure increment
       ------------------------------ */

    if (eqp_p->iswdyn <= 0) {
      if (idtvar >= 0 && isweep <= nswmpr && residual > eqp_p->epsrsm*rnormp) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          phia[c_id] = phi[c_id];
          phi[c_id] = phi[c_id] + eqp_p->relaxv*dphi[c_id];
        }
      }
      /* If it is the last sweep, update with the total increment */
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          phia[c_id] = phi[c_id];
          phi[c_id] = phi[c_id] + dphi[c_id];
        }
      }
    }
    else if (eqp_p->iswdyn == 1) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        phia[c_id] = phi[c_id];
        phi[c_id] = phi[c_id] + alph*dphi[c_id];
      }
    }
    else if (eqp_p->iswdyn >= 2) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        phia[c_id] = phi[c_id];
        phi[c_id] = phi[c_id] + alph*dphi[c_id] + beta*dphim1[c_id];
      }
    }

    /* Update the right hand side and update the residual
     *   rhs^{k+1} = - div(rho u^n) - D(dt, delta delta p^{k+1})
     * --------------------------------------------------------- */

    {
      int init = 1;
      int inc = 0;  /* by increment */
      if (   vp_param->iphydr == 1
          || vp_param->iifren == 1
          || vp_param->staggered == 1)
        inc = 1;    /* not by increment */

      if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
        cs_diffusion_potential(-1,
                               m,
                               fvq,
                               1,  /* init */
                               inc,
                               eqp_p->imrgra,
                               eqp_p->nswrgr,
                               eqp_p->imligr,
                               vp_param->iphydr,
                               eqp_p->iwgrec,
                               eqp_p->verbosity,
                               eqp_p->epsrgr,
                               eqp_p->climgr,
                               dfrcxt,
                               phi,
                               coefa_dp, coefb_dp,
                               coefaf_dp, coefbf_dp,
                               i_visc, b_visc,
                               c_visc,
                               rhs);

      else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
        cs_anisotropic_diffusion_potential(-1,
                                           m,
                                           fvq,
                                           init,
                                           inc,
                                           eqp_p->imrgra,
                                           eqp_p->nswrgr,
                                           eqp_p->imligr,
                                           eqp_p->ircflu,
                                           vp_param->iphydr,
                                           eqp_p->iwgrec,
                                           eqp_p->verbosity,
                                           eqp_p->epsrgr,
                                           eqp_p->climgr,
                                           dfrcxt,
                                           phi,
                                           coefa_dp, coefb_dp,
                                           coefaf_dp, coefbf_dp,
                                           i_visc, b_visc,
                                           vitenp,
                                           weighf, weighb,
                                           rhs);
    }

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      rhs[c_id] = - cpro_divu[c_id] - rhs[c_id] - rovsdt[c_id]*phi[c_id];
    }

    /* Convergence test */
    residual = sqrt(cs_gdot(n_cells, rhs, rhs));

    /*  Writing */
    sinfo->n_it += niterf;

    if (eqp_p->iwarni >= 2) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _(fmt_sweep_residual_relax),
                     f_p->name, isweep, residual, eqp_p->relaxv);

      cs_log_printf(CS_LOG_DEFAULT,
                    _(fmt_sweep_residual_info),
                     f_p->name, isweep, residual, rnormp, niterf);
    }

    isweep += 1;

  }  /* End of reconstruction loop
        -------------------------- */

  /* For logging */
  if (fabs(rnormp) > 0.)
    sinfo->res_norm = residual/rnormp;
  else
    sinfo->res_norm = 0.;

  /*  Writing */
  if (eqp_p->iwarni >= 1) {
    if (residual <= eqp_p->epsrsm*rnormp)
      cs_log_printf(CS_LOG_DEFAULT,
                    _(fmt_sweep_residual_info),
                    f_p->name, isweep-1, residual, rnormp, niterf);

    else if(isweep > nswmpr) {  /* non-convergence */
      cs_log_printf(CS_LOG_DEFAULT,
                    _("@\n"
                      "@ @@ Warning: %s (pressure correction step)\n"
                      "     =======\n"
                      "  Maximum number of iterations (%d) reached\n"),
                    f_p->name, nswmpr);
    }
  }

  /* Compute the indicator, taken the volume into account (L2 norm) or not */

  cs_field_t *f_err_est = cs_field_by_name_try("est_error_der_1");
  if (f_err_est != NULL) {
    cs_real_t *c_estim_der = f_err_est->val;
    const cs_real_t *restrict cell_vol = fvq->cell_vol;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_estim_der[c_id] = fabs(rhs[c_id]) / cell_vol[c_id];
  }
  f_err_est = cs_field_by_name_try("est_error_der_2");
  if (f_err_est != NULL) {
    cs_real_t *c_estim_der = f_err_est->val;
    const cs_real_t *restrict cell_vol = fvq->cell_vol;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_estim_der[c_id] = fabs(rhs[c_id]) / sqrt(cell_vol[c_id]);
  }

  /* Update the mass flux
     -------------------- */

  /* We cancel the face viscosity for coupled faces so as not to modify the
     boundary mass flux in the case of a pressure Dirichlet:
     pressure correction and filter are canceled. */

  if (cs_sat_coupling_n_couplings() > 0) {
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      if (bc_type[f_id] == CS_COUPLED_FD)
        b_visc[f_id] = 0.;
    }
  }

  {
    int inc  = 0;
    /* In case of hydrostatic pressure, inc is set to 1 to take explicit
       boundary conditions on the pressure (coefa) */
    if (   vp_param->iphydr == 1
        || vp_param->iifren == 1
        || vp_param->staggered == 1)
      inc = 1;

    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
      cs_face_diffusion_potential(-1,
                                  m,
                                  fvq,
                                  0,  /* init */
                                  inc,
                                  eqp_p->imrgra,
                                  eqp_p->nswrgr,
                                  eqp_p->imligr,
                                  vp_param->iphydr,
                                  eqp_p->iwgrec,
                                  eqp_p->verbosity,
                                  eqp_p->epsrgr,
                                  eqp_p->climgr,
                                  dfrcxt,
                                  phia,
                                  coefa_dp, coefb_dp,
                                  coefaf_dp, coefbf_dp,
                                  i_visc, b_visc, c_visc,
                                  imasfl, bmasfl);

    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
      cs_face_anisotropic_diffusion_potential(-1,
                                              m,
                                              fvq,
                                              0,  /* init */
                                              inc,
                                              eqp_p->imrgra,
                                              eqp_p->nswrgr,
                                              eqp_p->imligr,
                                              eqp_p->ircflu,
                                              vp_param->iphydr,
                                              eqp_p->iwgrec,
                                              eqp_p->verbosity,
                                              eqp_p->epsrgr,
                                              eqp_p->climgr,
                                              dfrcxt,
                                              phia,
                                              coefa_dp, coefb_dp,
                                              coefaf_dp, coefbf_dp,
                                              i_visc, b_visc,
                                              vitenp,
                                              weighf, weighb,
                                              imasfl, bmasfl);

    /* The last increment is not reconstructed so as to fulfill exactly
       the continuity equation (see theory guide). The value of dfrcxt has
       no importance. */

    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
      cs_face_diffusion_potential(-1,
                                  m,
                                  fvq,
                                  0,  /* init */
                                  0,  /* inc */
                                  eqp_p->imrgra,
                                  0,  /* nswrgr (no reconstruction) */
                                  eqp_p->imligr,
                                  vp_param->iphydr,
                                  eqp_p->iwgrec,
                                  eqp_p->verbosity,
                                  eqp_p->epsrgr,
                                  eqp_p->climgr,
                                  dfrcxt,
                                  dphi,
                                  coefa_dp, coefb_dp,
                                  coefaf_dp, coefbf_dp,
                                  i_visc, b_visc, c_visc,
                                  imasfl, bmasfl);

    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
      cs_face_anisotropic_diffusion_potential(-1,
                                              m,
                                              fvq,
                                              0, /* init */
                                              0, /* inc */
                                              eqp_p->imrgra,
                                              0, /* nswrgr (no reconstruction) */
                                              eqp_p->imligr,
                                              0, /* ircflu */
                                              vp_param->iphydr,
                                              eqp_p->iwgrec,
                                              eqp_p->verbosity,
                                              eqp_p->epsrgr,
                                              eqp_p->climgr,
                                              dfrcxt,
                                              dphi,
                                              coefa_dp, coefb_dp,
                                              coefaf_dp, coefbf_dp,
                                              i_visc, b_visc,
                                              vitenp,
                                              weighf, weighb,
                                              imasfl, bmasfl);

  }

  /* Update density
     -------------- */

  if (compressible_flag == 3) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      crom_eos[c_id] = crom_eos[c_id] + phi[c_id]/c2[c_id];
  }

  BFT_FREE(c2);

  /* Suppression of the grid hierarchy (release solver setup)
     --------------------------------- */

  cs_sles_free_native(f_p->id, NULL);

  BFT_FREE(dam);
  BFT_FREE(xam);

  /* Weakly compressible algorithm: semi analytic scheme
     2nd step solving a convection diffusion equation
     =================================================== */

  if (idilat == 5) {

    cs_real_t *ddphi;
    cs_real_3_t *coefar, *cofafr;
    cs_real_33_t *coefbr, *cofbfr;
    CS_MALLOC_HD(ddphi, n_cells_ext, cs_real_t, cs_alloc_mode);
    BFT_MALLOC(coefar, n_b_faces, cs_real_3_t);
    BFT_MALLOC(cofafr, n_b_faces, cs_real_3_t);
    BFT_MALLOC(coefbr, n_b_faces, cs_real_33_t);
    BFT_MALLOC(cofbfr, n_b_faces, cs_real_33_t);

    /* Convective flux: dt/rho grad(rho) */

    /* Dirichlet Boundary Condition on rho
       ----------------------------------- */

    cs_real_t *coefa_rho, *coefb_rho;
    BFT_MALLOC(coefa_rho, n_b_faces, cs_real_t);
    BFT_MALLOC(coefb_rho, n_b_faces, cs_real_t);

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      coefa_rho[f_id] = brom[f_id];
      coefb_rho[f_id] = 0.;
    }

    cs_halo_type_t halo_type = CS_HALO_STANDARD;
    cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;
    cs_gradient_type_by_imrgra(eqp_u->imrgra,
                               &gradient_type,
                               &halo_type);

    cs_gradient_scalar("Work array",
                       gradient_type,
                       halo_type,
                       1,             /* inc */
                       eqp_u->nswrgr,
                       0,             /* iphydp */
                       1,             /* w_stride */
                       eqp_p->verbosity,
                       eqp_u->imligr,
                       eqp_u->epsrgr,
                       eqp_u->climgr,
                       NULL,          /* f_ext */
                       coefa_rho, coefb_rho,
                       crom,
                       NULL,         /* c_weight */
                       NULL,         /* cpl */
                       gradp);

    BFT_FREE(coefa_rho);
    BFT_FREE(coefb_rho);

    /* dt/rho * grad rho */
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        wrk[c_id][j] = gradp[c_id][j] * dt[c_id] / crom[c_id];
      }
    }

    /* Viscosity */
    cs_face_viscosity(m,
                      fvq,
                      eqp_u->imvisf,
                      dt,
                      i_visc, b_visc);

    /* (dt/rho * grad rho) . S */

    /* Neumann boundary Conditions for the convective flux */

    const cs_real_t qimpv[3] = {0., 0., 0.};

    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION) {

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        cs_lnum_t c_id = b_face_cells[f_id];
        cs_real_t hint = dt[c_id] / b_dist[f_id];

        cs_boundary_conditions_set_neumann_vector
          (coefar[f_id], cofafr[f_id],
           coefbr[f_id], cofbfr[f_id],
           qimpv, hint);
      }

    }
    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) {

      /* Symmetric tensor diffusivity */

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        cs_lnum_t c_id = b_face_cells[f_id];

        cs_real_t visci[3][3], dist[3];
        const cs_real_t *n = b_face_normal[f_id];

        visci[0][0] = vitenp[c_id][0];
        visci[1][1] = vitenp[c_id][1];
        visci[2][2] = vitenp[c_id][2];
        visci[1][0] = vitenp[c_id][3];
        visci[0][1] = vitenp[c_id][3];
        visci[2][1] = vitenp[c_id][4];
        visci[1][2] = vitenp[c_id][4];
        visci[2][0] = vitenp[c_id][5];
        visci[0][2] = vitenp[c_id][5];

        dist[0] = b_face_cog[f_id][0] - cell_cen[c_id][0];
        dist[1] = b_face_cog[f_id][1] - cell_cen[c_id][1];
        dist[2] = b_face_cog[f_id][2] - cell_cen[c_id][2];

        /* ||Ki.S||^2 */
        cs_real_t viscis =   cs_math_pow2(  visci[0][0]*n[0]
                                          + visci[1][0]*n[1]
                                          + visci[2][0]*n[2])
                           + cs_math_pow2(  visci[0][1]*n[0]
                                          + visci[1][1]*n[1]
                                          + visci[2][1]*n[2])
                           + cs_math_pow2(  visci[0][2]*n[0]
                                          + visci[1][2]*n[1]
                                          + visci[2][2]*n[2]);

        /* IF.Ki.S */
        cs_real_t fikis
          = (  cs_math_3_dot_product(dist, visci[0]) * n[0]
             + cs_math_3_dot_product(dist, visci[1]) * n[1]
             + cs_math_3_dot_product(dist, visci[2]) * n[2]);

        /* Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
           NB: eps =1.d-1 must be consistent with vitens.f90 */

        fikis = fmax(fikis, 1.e-1*sqrt(viscis)*b_dist[f_id]);

        cs_real_t hint = viscis/b_face_surf[f_id]/fikis;

        cs_boundary_conditions_set_neumann_vector
          (coefar[f_id], cofafr[f_id],
           coefbr[f_id], cofbfr[f_id],
           qimpv, hint);

      }

    } /* end of boundary condition definitions */

    int iflmb0 = (cs_glob_ale > CS_ALE_NONE) ? 0 : 1;

    cs_mass_flux(m,
                 fvq,
                 -1,
                 0,  /* itypfl */
                 iflmb0,
                 1,  /* init */
                 1,  /* inc */
                 eqp_u->imrgra,
                 eqp_u->nswrgr,
                 eqp_u->imligr,
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 wrk,
                 coefar, coefbr,
                 velflx, velflb);

    /* Boundary condition for the pressure increment
       coefb, coefbf are those of the pressure */

    cs_real_t *coefa_dp2, *coefaf_dp2;
    BFT_MALLOC(coefa_dp2, n_b_faces, cs_real_t);
    BFT_MALLOC(coefaf_dp2, n_b_faces, cs_real_t);

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      coefa_dp2[f_id] = 0.;
      coefaf_dp2[f_id] = 0.;
    }

    /* Convective source term */

    cs_array_real_fill_zero(n_cells, rhs);

    cs_equation_param_t  eqp_loc = *eqp_p;

    /*  All boundary convective flux with upwind */

    eqp_loc.iconv  = 1;
    eqp_loc.istat  = -1;
    eqp_loc.idiff  = 0;
    eqp_loc.idifft = -1;
    eqp_loc.iswdyn = -1;
    eqp_loc.nswrgr = eqp_u->nswrgr;  /* FIXME: based on Fortran version,
                                        but value from pressure would
                                        seem more logical */
    eqp_loc.nswrsm = -1;
    eqp_loc.iwgrec = 0;
    eqp_loc.thetav = 1;
    eqp_loc.blend_st = 0; /* Warning, may be overwritten if a field */
    eqp_loc.epsilo = -1;
    eqp_loc.epsrsm = -1;

    cs_balance_scalar(idtvar,
                      -1,
                      0,     /* imucpp */
                      1,     /* imasac */
                      1,     /* inc */
                      &eqp_loc,
                      phi, phi,
                      coefa_dp2, coefb_p,
                      coefaf_dp2, coefbf_p,
                      velflx, velflb,
                      i_visc, b_visc,
                      NULL,  /* viscel */
                      NULL,  /* xcpp */
                      NULL,  /* weighf */
                      NULL,  /* weighb */
                      0,     /* icvflb; upwind scheme */
                      NULL,
                      rhs);

    /* Initialization of the variable to solve */
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      rovsdt[c_id] = 340.0/dt[c_id] * cell_f_vol[c_id];
      dphi[c_id]   = 0.;
      ddphi[c_id]  = 0.;
      rhs[c_id]    = - rhs[c_id];
    }

    /* Solve the convection diffusion equation */

    const char var_name[] = "Pr compress";

    cs_sles_push(f_p->id, var_name);

    eqp_loc = *eqp_p;

    eqp_loc.iconv  = 1;
    eqp_loc.istat  = -1;
    eqp_loc.icoupl = -1;
    eqp_loc.ndircl = 0;    /* to reinforce the diagonal */
    eqp_loc.idiff  = 1;
    eqp_loc.idifft = -1;
    eqp_loc.iwgrec = 0;    /* Warning, may be overwritten if a field */
    eqp_loc.blend_st = 0;  /* Warning, may be overwritten if a field */

    cs_equation_iterative_solve_scalar(idtvar,
                                       iterns,
                                       f_p->id,
                                       var_name,
                                       0,      /* iescap */
                                       0,      /* imucpp */
                                       -1,     /* normp */
                                       &eqp_loc,
                                       dphi, dphi,
                                       coefa_dp2, coefb_p,
                                       coefaf_dp2, coefbf_p,
                                       velflx, velflb,
                                       i_visc, b_visc,
                                       i_visc, b_visc,
                                       NULL,   /* viscel */
                                       weighf, weighb,
                                       0,      /* icvflb (upwind conv. flux) */
                                       NULL,   /* icvfli */
                                       rovsdt,
                                       rhs,
                                       dphi, ddphi,
                                       NULL,   /* xcpp */
                                       NULL);  /* eswork */

    cs_sles_pop(f_p->id);

    /* Update the pressure increment */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      phi[c_id] += dphi[c_id];
      /* Remove the last increment */
      dphi[c_id] -= ddphi[c_id];
    }

    /* Update the mass flux */

    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
      cs_face_diffusion_potential(-1,
                                  m,
                                  fvq,
                                  0,  /* init */
                                  1,  /* inc */
                                  eqp_p->imrgra,
                                  eqp_p->nswrgr,
                                  eqp_p->imligr,
                                  vp_param->iphydr,
                                  eqp_p->iwgrec,
                                  eqp_p->verbosity,
                                  eqp_p->epsrgr,
                                  eqp_p->climgr,
                                  dfrcxt,
                                  dphi,
                                  coefa_dp2, coefb_p,
                                  coefaf_dp2, coefbf_p,
                                  i_visc, b_visc,
                                  dt,
                                  imasfl, bmasfl);

    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
      cs_face_anisotropic_diffusion_potential(-1,
                                              m,
                                              fvq,
                                              0,  /* init */
                                              1,  /* inc */
                                              eqp_p->imrgra,
                                              eqp_p->nswrgr,
                                              eqp_p->imligr,
                                              eqp_p->ircflu,
                                              vp_param->iphydr,
                                              eqp_p->iwgrec,
                                              eqp_p->verbosity,
                                              eqp_p->epsrgr,
                                              eqp_p->climgr,
                                              dfrcxt,
                                              dphi,
                                              coefa_dp2, coefb_p,
                                              coefaf_dp2, coefbf_p,
                                              i_visc, b_visc,
                                              vitenp,
                                              weighf, weighb,
                                              imasfl, bmasfl);

    /* The last increment is not reconstructed so as to fulfill exactly
       the continuity equation (see theory guide). The value of dfrcxt has
       no importance. */

    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
      cs_face_diffusion_potential(-1,
                                  m,
                                  fvq,
                                  0,  /* init */
                                  0,  /* inc */
                                  eqp_p->imrgra,
                                  0,  /* nswrgr (no reconstruction) */
                                  eqp_p->imligr,
                                  vp_param->iphydr,
                                  eqp_p->iwgrec,
                                  eqp_p->verbosity,
                                  eqp_p->epsrgr,
                                  eqp_p->climgr,
                                  dfrcxt,
                                  ddphi,
                                  coefa_dp2, coefb_p,
                                  coefaf_dp2,coefbf_p,
                                  i_visc, b_visc,
                                  dt,
                                  imasfl, bmasfl);

    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
      cs_face_anisotropic_diffusion_potential(-1,
                                              m,
                                              fvq,
                                              0, /* init */
                                              0, /* inc */
                                              eqp_p->imrgra,
                                              0, /* nswrgr (no reconstruction) */
                                              eqp_p->imligr,
                                              0, /* ircflu */
                                              vp_param->iphydr,
                                              eqp_p->iwgrec,
                                              eqp_p->verbosity,
                                              eqp_p->epsrgr,
                                              eqp_p->climgr,
                                              dfrcxt,
                                              ddphi,
                                              coefa_dp2, coefb_p,
                                              coefaf_dp2, coefbf_p,
                                              i_visc, b_visc,
                                              vitenp,
                                              weighf, weighb,
                                              imasfl, bmasfl);

    /* Free memory */
    CS_FREE_HD(ddphi);
    BFT_FREE(coefa_dp2);
    BFT_FREE(coefaf_dp2);
    BFT_FREE(coefar);
    BFT_FREE(coefbr);
    BFT_FREE(cofafr);
    BFT_FREE(cofbfr);

  } /* End if weaky compressible algorithm (idilat = 5) */

  BFT_FREE(velflx);
  BFT_FREE(velflb);

  /* Update the pressure field
     ========================= */

  // Pressure at the last sub-iteration

  cs_real_t *pk1;
  BFT_MALLOC(pk1, n_cells_ext, cs_real_t);

  if (idtvar < 0) {
    const cs_real_t relaxv = eqp_p->relaxv;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      pk1[c_id] = cvar_pr[c_id];
      cvar_pr[c_id] += relaxv*phi[c_id];
    }
  }
  else {
    if (vp_param->staggered == 0) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        pk1[c_id] = cvar_pr[c_id];
        cvar_pr[c_id] += phi[c_id];
      }
    }
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        pk1[c_id] = cvar_pr[c_id];
        cvar_pr[c_id] = phi[c_id];
      }
    }
  }

  /* Transformation of volume fluxes into mass fluxes */

  /* Update the density when solving the Helmholtz equation */
  if (ieos != CS_EOS_NONE) {
    cs_real_t *cpro_rho_mass = cs_field_by_name("density_mass")->val;
    cs_real_t drop;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      drop = (_coef * cvar_pr[c_id] - (_coef - 1.) * cvara_pr[c_id]
           - pk1[c_id]) * dc2[c_id];
      crom_eos[c_id] += drop;
      cpro_rho_mass[c_id] = crom_eos[c_id];
    }
  }

  BFT_FREE(dc2);

  /* Compute the isobaric heat capacity if needed */
  cs_field_t  *f_cv = cs_field_by_name_try("isobaric_heat_capacity");
  if (f_cv != NULL) {
     cs_thermal_model_cv(f_cv->val);
  }

  /* Correction of the temperature and yv after the pressure */
  if (   idilat == 2
      && thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY
      && ieos == CS_EOS_MOIST_AIR) {
    /* Last argument is the method used, 1 for the newton, 2 for the
     * pressure increment (explicit correction)*/
    cs_thermal_model_newton_t(2,
                              pk1,
                              cvar_th,
                              cvar_pr,
                              cvara_pr,
                              yw,
                              yv,
                              tempk);
  }

  BFT_FREE(pk1);

  /* Save some information */
  if (idilat == 2 && ieos != CS_EOS_NONE) {
    /* CFL conditions related to the pressure equation */
    cs_real_t *cflp = NULL;
    cs_field_t *f_cflp = cs_field_by_name_try("cfl_p");
    if (f_cflp != NULL) {
      cflp = f_cflp->val;
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
        cflp[c_id] = 0.;
      }

      /* Parallelism and periodicity */
      cs_mesh_sync_var_vect((cs_real_t *)wrk2);

      /* Call the function to compute the cfl condition related to the pressure
       * equation */
      cs_thermal_model_cflp(croma,
                            wrk2,
                            cvara_pr,
                            imasfl,
                            cflp);
    }
    if (kinetic_st == 1) {
      /* Save rho k-1 , est ce le bon endroit ?*/
      cs_real_t *rho_k_prev = cs_field_by_name("rho_k_prev")->val;
      // Test compute kinetic energy st
      cs_thermal_model_kinetic_st_finalize(rho_k_prev, crom_eos);
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
        rho_k_prev[c_id] =  crom_eos[c_id];
      }
    }
  }

  if (idilat == 4) {

    const cs_real_t *restrict fw = fvq->weight;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      bmasfl[f_id] *= brom[f_id];
    }

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id_0 = i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = i_face_cells[f_id][1];

      // FIXME: should be coherent with the convective scheme of the species...
      imasfl[f_id] *=   fw[f_id]*crom[c_id_0]
                      + (1.-fw[f_id])*crom[c_id_1];
    }

  }

  /*  Free memory */
  BFT_FREE(taui);
  BFT_FREE(taub);
  BFT_FREE(wrk2);
  BFT_FREE(wrk);
  BFT_FREE(res);
  BFT_FREE(phia);
  BFT_FREE(dphi);
  BFT_FREE(_cpro_divu);
  BFT_FREE(gradp);
  CS_FREE_HD(rhs);
  BFT_FREE(rovsdt);
  BFT_FREE(weighf);
  BFT_FREE(weighb);
  BFT_FREE(adxk);
  BFT_FREE(adxkm1);
  BFT_FREE(dphim1);
  BFT_FREE(rhs0);
  BFT_FREE(xdtsro);
  BFT_FREE(tpusro);

  BFT_FREE(cpro_rho_tc);
  BFT_FREE(bpro_rho_tc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the pressure correction step of the Navier-Stokes equations
 *        for incompressible or slightly compressible flows.
 *
 * \remark:
 * - CDO face-based scheme is used to solve the Poisson equation.
 *  \f[
 *     - DIV \cdot Hdg \cdot GRAD \left(\delta p \right) =
 * - \divs \left( \rho \vect{\widetilde{u}}\right)
 * \f]
 * The mass flux is then updated as follows:
 * \f[
 *  \dot{m}^{n+1} = \dot{m}^{n}
 *                 - Hodge \cdot GRAD \left(\delta p \right)
 * \f]
 *
 * \param[in] vel     velocity
 * \param[in] coefav  boundary condition array for the variable (explicit part)
 * \param[in] coefbv  boundary condition array for the variable (implicit part)
 */
/*----------------------------------------------------------------------------*/

static void
_pressure_correction_cdo(cs_real_t  vel[restrict][3],
                         cs_real_t  coefav[restrict][3],
                         cs_real_t  coefbv[restrict][3][3])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  const cs_cdo_quantities_t  *quant = cs_glob_domain->cdo_quantities;
  const cs_cdo_connect_t  *connect = cs_glob_domain->connect;

  const cs_lnum_t  n_cells = m->n_cells;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;

  cs_field_t  *f_p = CS_F_(p);
  const cs_field_t  *f_vel = CS_F_(vel);
  assert((cs_real_t *)vel == f_vel->val);

  const cs_equation_param_t  *eqp_u = cs_field_get_equation_param_const(f_vel);
  const cs_equation_param_t  *eqp_p = cs_field_get_equation_param_const(f_p);

  const cs_velocity_pressure_model_t  *vp_model
    = cs_glob_velocity_pressure_model;
  const cs_velocity_pressure_param_t  *vp_param
    = cs_glob_velocity_pressure_param;

  const int idilat = vp_model->idilat;
  const int iphydr = vp_param->iphydr;
  const int compressible_flag
    = cs_glob_physical_model_flag[CS_COMPRESSIBLE];

  cs_pressure_correction_cdo_t *prcdo = cs_pressure_correction_cdo;
  if (prcdo == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_prcdo));

  cs_equation_t  *eq_dp = prcdo->pressure_incr;

  /* Allocate temporary arrays */

  cs_real_3_t  *wrk = NULL;
  BFT_MALLOC(wrk, n_cells_ext, cs_real_3_t);

  cs_field_t  *f_dp = cs_field_by_id(eq_dp->field_id);
  cs_real_t  *phi = f_dp->val;
  cs_real_t  *divu = prcdo->div_st;

  /* Associate pointers to pressure diffusion coefficients */

  cs_real_t  *restrict dt = CS_F_(dt)->val;
  cs_real_t  *c_visc = dt;

  /* Index of the field */

  cs_solving_info_t  *sinfo
    = cs_field_get_key_struct_ptr(f_p,
                                  cs_field_key_id("solving_info"));

  /* Physical quantities */

  cs_real_t  *crom_eos = CS_F_(rho)->val;
  cs_real_t  *brom_eos = CS_F_(rho_b)->val;
  cs_real_t  *crom = crom_eos;
  cs_real_t  *brom = brom_eos;

  cs_real_t  *cvar_pr = f_p->vals[0];

  const int  kimasf = cs_field_key_id("inner_mass_flux_id");
  const int  kbmasf = cs_field_key_id("boundary_mass_flux_id");

  cs_real_t  *imasfl = cs_field_by_id(cs_field_get_key_int(f_p, kimasf))->val;
  cs_real_t  *bmasfl = cs_field_by_id(cs_field_get_key_int(f_p, kbmasf))->val;

  cs_real_t  *ipotfl = prcdo->inner_potential_flux;
  cs_real_t  *bpotfl = prcdo->bdy_potential_flux;

  /* Solving options */

  const cs_time_step_t  *ts = cs_glob_time_step;

  cs_real_3_t  *gradp = (cs_real_3_t*) prcdo->pressure_gradient->val;

  if (   ts->nt_cur <= ts->nt_ini
      && iphydr == 0 && idilat <= 1
      && compressible_flag < 0
      && cs_restart_present() == false) {

    cs_array_real_fill_zero(3*n_cells, prcdo->pressure_gradient->val);
    cs_array_real_fill_zero(n_i_faces, ipotfl);
    cs_array_real_fill_zero(n_b_faces, bpotfl);

  }

  /* Building the linear system to solve.
   * ==================================== */

  const cs_real_t arak = vp_param->arak;

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (cs_lnum_t j = 0; j < 3; j++)
      wrk[c_id][j] = vel[c_id][j] + arak*gradp[c_id][j]*dt[c_id]/crom[c_id];
  }

  /* Sync for parallelism and periodicity */

  cs_mesh_sync_var_vect((cs_real_t *)wrk);

  {
    int inc = 1;
    int iflmb0 = (cs_glob_ale > CS_ALE_NONE) ? 0 : 1;
    int itypfl = 1;

    cs_mass_flux(m,
                 fvq,
                 -1,            /* f_id */
                 itypfl,        /* =1 --> multiply by rho */
                 iflmb0,
                 1,             /* init */
                 inc,
                 eqp_u->imrgra,
                 eqp_u->nswrgr,
                 eqp_u->imligr,
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 wrk,
                 coefav, coefbv,
                 imasfl, bmasfl);
  }

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    imasfl[f_id] += arak*ipotfl[f_id];

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    bmasfl[f_id] += arak*bpotfl[f_id];

  /*
   * Solving
   * =========================================== */

  /* Variables are set to 0
   *   phi  is the increment of the pressure
   *   divu is the initial divergence of the predicted mass flux */

  cs_array_real_fill_zero(n_cells, phi);

  /* Initial divergence (= rhs of the system to solve) */

  cs_divergence(m,
                1,  /* init */
                imasfl,
                bmasfl,
                divu);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    divu[c_id] *= -1;

  /* Compute and set the right-hand side residual */

  cs_real_t residual = sqrt(cs_gdot(n_cells, divu, divu));

  sinfo->rhs_norm = residual;

  /* Solve the pressure increment with CDO
     ------------------------------------- */

  cs_equation_solve_steady_state(m, eq_dp);

  /* Update the mass flux and potential flux
     --------------------------------------- */

  cs_real_t *diff_flux = NULL;
  BFT_MALLOC(diff_flux, quant->n_faces, cs_real_t);

  cs_equation_compute_diffusive_flux(eq_dp,
                                     NULL, /* eqp --> default*/
                                     NULL, /* diff_pty --> default */
                                     NULL, /* dof_values --> default */
                                     NULL, /* cell_values --> default */
                                     cs_flag_primal_face,
                                     cs_glob_time_step->t_cur,
                                     diff_flux);

  for (cs_lnum_t f_id = 0; f_id < quant->n_i_faces; f_id++) {
    imasfl[f_id] += diff_flux[f_id];
    ipotfl[f_id] += diff_flux[f_id];
  }

  for (cs_lnum_t f_id = 0; f_id < quant->n_b_faces; f_id++) {
    bmasfl[f_id] += diff_flux[f_id + quant->n_i_faces];
    bpotfl[f_id] += diff_flux[f_id + quant->n_i_faces];
  }

  /* Reconstruct the cell gradient
     ----------------------------- */

  cs_real_3_t  *grddp_c = (cs_real_3_t*) prcdo->pressure_incr_gradient->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_lnum_t  s = connect->c2f->idx[c_id], e = connect->c2f->idx[c_id+1];
    const cs_lnum_t  *c2f_ids = connect->c2f->ids + s;
    const short int  *c2f_sgn = connect->c2f->sgn + s;

    const cs_lnum_t  n_fc = e - s;
    const cs_real_t  *cell_center = quant->cell_centers + 3*c_id;
    cs_real_3_t  grddp_reco = {0.0, 0.0, 0.0};

    for (short int j = 0; j < n_fc; j++) {

      cs_lnum_t  f_id = c2f_ids[j];

      const cs_real_t *face_center  = (f_id < quant->n_i_faces) ?
        quant->i_face_center + 3*f_id :
        quant->b_face_center + 3*(f_id - quant->n_i_faces);

      for (short int k = 0; k < 3; k++)
        grddp_reco[k] +=
          - c2f_sgn[j] * diff_flux[f_id] * (face_center[k] - cell_center[k]);

    } /* Loop on cell faces */

    const cs_real_t  ccoef = 1./(c_visc[c_id]*quant->cell_vol[c_id]);
    for (short int k = 0; k < 3; k++) {
      const cs_real_t  incr_k = ccoef*grddp_reco[k];
      grddp_c[c_id][k] = incr_k;
      gradp[c_id][k] += incr_k;
    }

  } /* Loop on cells */

  /*  Free memory */

  BFT_FREE(wrk);
  BFT_FREE(diff_flux);

#if defined(DEBUG) && !defined(NDEBUG) && CS_PRESSURE_CORRECTION_CDO_DBG > 0
  cs_real_t *res = NULL;
  BFT_MALLOC(res, n_cells_ext, cs_real_t);

  cs_divergence(m, 1, imasfl, bmasfl, res);

  cs_real_t rnormp = sqrt(cs_gdot(n_cells, res, res));

  cs_log_printf(CS_LOG_DEFAULT,
                ">> Divergence of mass flux after correction step: %15.7e \n",
                rnormp);

  BFT_FREE(res);
#endif

  /* Update the pressure field
     ------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    cvar_pr[c_id] += phi[c_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a cs_pressure_correction_cdo_t structure
 *
 * \return a pointer to a newly allocated
 *         \ref cs_pressure_correction_cdo_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_pressure_correction_cdo_t *
_pressure_correction_cdo_create(void)
{
  cs_pressure_correction_cdo_t *prcdo = NULL;

  BFT_MALLOC(prcdo, 1, cs_pressure_correction_cdo_t);

  /* Equation
     -------- */

  prcdo->pressure_incr = NULL;

  /* Fields
     ------ */

  prcdo->pressure_incr_gradient = NULL;
  prcdo->pressure_gradient = NULL;

  /* Other arrays
     ------------ */

  prcdo->div_st = NULL;
  prcdo->inner_potential_flux = NULL;
  prcdo->bdy_potential_flux = NULL;
  prcdo->bdy_pressure_incr = NULL;

  return prcdo;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate the pressure increment solving with Legacy FV
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_fv_activate(void)
{
  if (cs_pressure_correction_cdo_active)
    bft_error (__FILE__, __LINE__, 0,
               _("\n The pressure correction step is treated by CDO,"
                 "\n  Check the pressure correction model"));

  if (CS_F_(p) != NULL) {
    cs_field_t *f = cs_field_create("pressure_increment",
                                    CS_FIELD_INTENSIVE,
                                    CS_MESH_LOCATION_CELLS,
                                    1,
                                    false);
    cs_field_set_key_int(f,
                         cs_field_key_id("parent_field_id"),
                         CS_F_(p)->id);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate the pressure increment solving with CDO
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_cdo_activate(void)
{
  if (cs_pressure_correction_cdo_active)
    return;

  cs_pressure_correction_cdo_active = true;

  cs_pressure_correction_cdo_t *prcdo = NULL;

  if (cs_pressure_correction_cdo == NULL)
    cs_pressure_correction_cdo = _pressure_correction_cdo_create();

  prcdo = cs_pressure_correction_cdo;
  assert(prcdo != NULL);

  /* Activate the CDO module along with the FV framework */

  cs_param_cdo_mode_set(CS_PARAM_CDO_MODE_WITH_FV);

  /* Add a new equation related to the pressure correction with a variable
   * field named "pressure_increment" */

  cs_equation_t  *eq =
    cs_equation_add("pressure_increment", /* equation name */
                    "pressure_increment", /* associated variable field name */
                    CS_EQUATION_TYPE_PREDEFINED,
                    1,                        /* dimension of the unknown */
                    CS_PARAM_BC_HMG_NEUMANN); /* default boundary */

  cs_equation_param_t  *eqp = cs_equation_get_param(eq);

  /* Predefined settings */

  cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");

  /* BC settings */

  cs_equation_param_set(eqp, CS_EQKEY_BC_ENFORCEMENT, "algebraic");

  /* System to solve is SPD by construction */

  cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "fcg");
  cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
  cs_equation_param_set(eqp, CS_EQKEY_AMG_TYPE, "k_cycle");

  cs_equation_param_set(eqp, CS_EQKEY_ITSOL_MAX_ITER, "2500");
  cs_equation_param_set(eqp, CS_EQKEY_ITSOL_RTOL, "1e-5");
  cs_equation_param_set(eqp, CS_EQKEY_ITSOL_RESNORM_TYPE, "filtered");

  prcdo->pressure_incr = eq;    /* Keep a link the equation pointer */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the pressure increment, either FV or CDO
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_model_activate(void)
{
  const cs_velocity_pressure_model_t
    *vp_model = cs_glob_velocity_pressure_model;

  if (vp_model->iprcdo > 0)
    cs_pressure_correction_cdo_activate();
  else
    cs_pressure_correction_fv_activate();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if pressure solving with CDO is activated
 *
 * \return true if solving with CDO is requested, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_pressure_correction_cdo_is_activated(void)
{
  if (cs_pressure_correction_cdo_active)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the pressure increment equation
 *         At this stage, numerical settings should be completely determined
 *         but connectivity and geometrical information is not yet available.
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_cdo_init_setup(void)
{
  cs_pressure_correction_cdo_t *prcdo = cs_pressure_correction_cdo;

  if (prcdo == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_prcdo));

  cs_equation_t *eq = prcdo->pressure_incr;
  cs_equation_param_t* eqp = cs_equation_get_param(eq);

  /* Add the variable field */

  cs_equation_predefined_create_field(1, eq); /* 1 = keep two states */

  /* Associate the pressure increment variable field to the pressure field */

  cs_field_set_key_int(cs_field_by_id(eq->field_id),
                       cs_field_key_id("parent_field_id"),
                       CS_F_(p)->id);

  /* Additional fields */

  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");
  const int  field_post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

  prcdo->pressure_gradient =
    cs_field_find_or_create("pressure_gradient",
                            CS_FIELD_INTENSIVE,
                            CS_MESH_LOCATION_CELLS,
                            3,
                            false);

  cs_field_set_key_int(prcdo->pressure_gradient, post_key, 1);
  cs_field_set_key_int(prcdo->pressure_gradient, log_key, field_post_flag);

  prcdo->pressure_incr_gradient =
    cs_field_find_or_create("pressure_increment_gradient",
                            CS_FIELD_INTENSIVE,
                            CS_MESH_LOCATION_CELLS,
                            3,
                            false);

  cs_field_set_key_int(prcdo->pressure_incr_gradient, post_key, 1);
  cs_field_set_key_int(prcdo->pressure_incr_gradient, log_key, field_post_flag);

 /* Activate the diffusion term
  * ---------------------------
  * Diffusivity coefficent for pressure correction equation is represented by
  * time step.Thus the diffusion is added after the field dt is created.
  */

  cs_property_t *pty = cs_property_by_name("time_step");
  cs_property_def_by_field(pty, cs_field_by_name("dt"));
  cs_equation_add_diffusion(eqp, pty);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize setting-up the pressure increment equation
 *         At this stage, numerical settings should be completely determined
 *
 * \param[in] domain     pointer to a cs_domaint_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_cdo_finalize_setup(const cs_domain_t   *domain)
{
  cs_pressure_correction_cdo_t *prcdo = cs_pressure_correction_cdo;

  if (prcdo == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_prcdo));

  cs_equation_t *eq = prcdo->pressure_incr;
  cs_equation_param_t* eqp = cs_equation_get_param(eq);

  /* Allocate useful arrays
     --------------------- */

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  BFT_MALLOC(prcdo->div_st, n_cells_ext, cs_real_t);
  BFT_MALLOC(prcdo->inner_potential_flux, n_i_faces, cs_real_t);
  BFT_MALLOC(prcdo->bdy_potential_flux, n_b_faces, cs_real_t);
  BFT_MALLOC(prcdo->bdy_pressure_incr, n_b_faces, cs_real_t);

  /* Affect source term for the equation
     ----------------------------------- */

  cs_array_real_fill_zero(n_cells_ext, prcdo->div_st);

  cs_equation_add_source_term_by_array(eqp,
                                       NULL,  /* all cells */
                                       cs_flag_primal_cell,
                                       prcdo->div_st,
                                       false, /* is owner */
                                       true); /* full length */

  /* Handle the boundary conditions for the correction step
     ------------------------------------------------------ */

  for (int i = 0; i < domain->boundaries->n_boundaries; i++) {

    const int z_id = domain->boundaries->zone_ids[i];
    const cs_zone_t *z = cs_boundary_zone_by_id(z_id);

    cs_boundary_type_t b_type = domain->boundaries->types[i];

    // TODO handle more types of pressure boundary conditions
    if (b_type & CS_BOUNDARY_OUTLET || b_type & CS_BOUNDARY_IMPOSED_P) {

      for (cs_lnum_t iel = 0; iel < z->n_elts; iel++) {
        cs_lnum_t f_id = z->elt_ids[iel];
        prcdo->bdy_pressure_incr[f_id] = 0.0;
      }

      cs_equation_add_bc_by_array(eqp,
                                  CS_PARAM_BC_DIRICHLET,
                                  z->name,
                                  cs_flag_primal_face,
                                  prcdo->bdy_pressure_incr,
                                  false, /* Do not transfer ownership */
                                  true); /* full length */

    }

  } /* Loop on pressure definitions */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the pressure correction
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_cdo_destroy_all(void)
{
  cs_pressure_correction_cdo_t  *prcdo = cs_pressure_correction_cdo;
  if (prcdo == NULL)
    return;

  BFT_FREE(prcdo->div_st);
  BFT_FREE(prcdo->inner_potential_flux);
  BFT_FREE(prcdo->bdy_potential_flux);
  BFT_FREE(prcdo->bdy_pressure_incr);

  BFT_FREE(prcdo);

  cs_pressure_correction_cdo = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the pressure correction step of the Navier-Stokes equations
 *        for incompressible or slightly compressible flows.
 *
 * This function solves the following Poisson equation on the pressure:
 * \f[
 *     D \left( \Delta t, \delta p \right) =
 * \divs \left( \rho \vect{\widetilde{u}}\right)
 *     - \Gamma^n
 *     + \dfrac{\rho^n - \rho^{n-1}}{\Delta t}
 * \f]
 *
 * Either the legacy FV method or a CDO face-based scheme is used.
 *
 * For the legacy case, the mass flux is  updated as follows:
 * \f[
 *  \dot{m}^{n+1}_\ij = \dot{m}^{n}_\ij
 *                    - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
 * \f]
 *
 * \remark:
 * - an iterative process is used to solve the Poisson equation.
 * - if the arak coefficient is set to 1, the the Rhie & Chow filter is
 *   activated.
 *
 * Please refer to the
 * <a href="../../theory.pdf#resopv"><b>resopv</b></a>
 * section of the theory guide for more information.
 *
 * \param[in]       iterns    Navier-Stokes iteration number
 * \param[in]       nfbpcd    number of faces with condensation source term
 * \param[in]       ncmast    number of cells with condensation source terms
 * \param[in]       ifbpcd    index of faces with condensation source term
 * \param[in]       ltmast    list of cells with condensation source terms
 *                            (1 to n numbering)
 * \param[in]       isostd    indicator of standard outlet and index
 *                            of the reference outlet face
 * \param[in]       vel       velocity
 * \param[in, out]  da_uu     velocity matrix
 * \param[in]       coefav    boundary condition array for the variable
 *                            (explicit part)
 * \param[in]       coefbv    boundary condition array for the variable
 *                            (implicit part)
 * \param[in]       coefa_dp  boundary conditions for the pressure increment
 * \param[in]       coefb_dp  boundary conditions for the pressure increment
 * \param[in]       spcond    variable value associated to the condensation
 *                            source term (for ivar=ipr, spcond is the
 *                            flow rate
 *                            \f$ \Gamma_{s,cond}^n \f$)
 * \param[in]       svcond    variable value associated to the condensation
 *                            source term (for ivar=ipr, svcond is the flow rate
 *                            \f$ \Gamma_{v, cond}^n \f$)
 * \param[in]       frcxt     external forces making hydrostatic pressure
 * \param[in]       dfrcxt    variation of the external forces
 *                            composing the hydrostatic pressure
 * \param[in]       i_visc    visc*surface/dist aux faces internes
 * \param[in]       b_visc    visc*surface/dist aux faces de bord
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction(int        iterns,
                       cs_lnum_t  nfbpcd,
                       cs_lnum_t  ncmast,
                       cs_lnum_t  ifbpcd[],
                       cs_lnum_t  ltmast[],
                       const int  isostd[],
                       cs_real_t  vel[restrict][3],
                       cs_real_t  da_uu[restrict][6],
                       cs_real_t  coefav[restrict][3],
                       cs_real_t  coefbv[restrict][3][3],
                       cs_real_t  coefa_dp[restrict],
                       cs_real_t  coefb_dp[restrict],
                       cs_real_t  spcond[restrict],
                       cs_real_t  svcond[restrict],
                       cs_real_t  frcxt[restrict][3],
                       cs_real_t  dfrcxt[restrict][3],
                       cs_real_t  i_visc[restrict],
                       cs_real_t  b_visc[restrict])
{
 const cs_velocity_pressure_model_t  *vp_model
    = cs_glob_velocity_pressure_model;

 if (vp_model->iprcdo == 0)
   _pressure_correction_fv(iterns,
                           nfbpcd,
                           ncmast,
                           ifbpcd,
                           ltmast,
                           isostd,
                           vel,
                           da_uu,
                           coefav,
                           coefbv,
                           coefa_dp,
                           coefb_dp,
                           spcond,
                           svcond,
                           frcxt,
                           dfrcxt,
                           i_visc,
                           b_visc);
 else
   _pressure_correction_cdo(vel,
                            coefav,
                            coefbv);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
