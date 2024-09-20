/*============================================================================
 * Pressure correction.
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
#include "cs_base_accel.h"
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
#include "cs_wall_condensation.h"

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

static cs_pressure_correction_cdo_t *cs_pressure_correction_cdo = nullptr;

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

  int *c_disable_flag = mq->c_disable_flag;
  cs_lnum_t has_dc = mq->has_disable_flag;

  const int ksinfo = cs_field_key_id("solving_info");
  cs_field_t *f = cs_field_by_name_try("hydrostatic_pressure");
  cs_solving_info_t *sinfo =
    static_cast<cs_solving_info_t*>(cs_field_get_key_struct_ptr(f, ksinfo));
  const cs_equation_param_t *eqp_pr = cs_field_get_equation_param_const(f);

  cs_dispatch_context ctx, ctx_c;
#if defined(HAVE_CUDA)
  ctx_c.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  if (iterns < 2)
    sinfo->n_it = 0;

  cs_field_bc_coeffs_t *bc_coeffs_hp = f->bc_coeffs;
  cs_real_t *cofaf_hp = bc_coeffs_hp->af;
  cs_real_t *cofbf_hp = bc_coeffs_hp->bf;
  cs_real_t *coefa_hp = bc_coeffs_hp->a;
  cs_real_t *coefb_hp = bc_coeffs_hp->b;

  /* Check for variation of the hydrostatic pressure at outlet
   *
   * We check if the source term has changed. We exit directly
   * if we do not have standard outlets.
   * The precisiton for tests if more or less arbitrary. */

  int *ical;
  CS_MALLOC_HD(ical, 1, int, cs_alloc_mode); // allocation on GPU
  *ical = 0.;

  const cs_real_t precab = 1.e2*cs_math_epzero;
  const cs_real_t precre = sqrt(cs_math_epzero);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    cs_real_t rnrmf = cs_math_3_square_norm(frcxt[c_id]);
    cs_real_t rnrmdf = cs_math_3_square_norm(dfrcxt[c_id]);

    if ((rnrmdf >= precre*rnrmf) && (rnrmdf >= precab)) {
      *ical = 1;
    }
  });

  ctx.wait();

  /* copy device to host for the next parall sum on ical */
  int _ical;
  #if defined(HAVE_ACCEL)
    cs_copy_d2h(&_ical, ical, sizeof(int));
  #else
    _ical = *ical;
  #endif

  cs_parall_sum(1, CS_INT_TYPE, &_ical);
  if ((_ical == 0) && (cs_glob_atmo_option->open_bcs_treatment == 0)) {
    *indhyd = 0;
    CS_FREE_HD(ical);
    return;
  }

  CS_FREE_HD(ical);

  if (cs_log_default_is_active() || eqp_pr->verbosity > 0)
    bft_printf("  Hydrostatic pressure computation:\n"
               "    updating the Dirichlets at the end"
               " (_hydrostatic_pressure_compute)\n");

  *indhyd = 1;

  cs_real_3_t *next_fext;
  CS_MALLOC_HD(next_fext, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  cs_real_t *rovsdt = nullptr, *viscce = nullptr;
  CS_MALLOC_HD(rovsdt, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(viscce, n_cells_ext, cs_real_t, cs_alloc_mode);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    const int c_act = 1 - (has_dc * c_disable_flag[has_dc * c_id]);

    rovsdt[c_id] = 0.0;
    dphi[c_id] = 0.;
    viscce[c_id] = 1.0;

    for (cs_lnum_t ii = 0; ii < 3; ii++)
      next_fext[c_id][ii] = frcxt[c_id][ii]*c_act + dfrcxt[c_id][ii];
  });

  ctx.wait();
  cs_mesh_sync_var_vect((cs_real_t*)next_fext);

  /* Prepare matrix and boundary conditions
     -------------------------------------- */

  ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

    /* Neumann BC (qimp=0 --> a=af=bf=0, b=1) */

    // Gradient BCs
    coefa_hp[face_id] = 0.;
    coefb_hp[face_id] = 1.;

    // Flux BCs
    cofaf_hp[face_id] = 0.;
    cofbf_hp[face_id] = 0.;

  });

  ctx_c.wait();

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
                           eqp_pr->theta,
                           0,
                           bc_coeffs_hp,
                           &rovsdt[0],
                           iflux,
                           bflux,
                           i_visc,
                           b_visc,
                           nullptr,
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
                    cofbf_hp,
                    iflux,
                    bflux,
                    i_visc,
                    b_visc,
                    viscce,
                    viscce,
                    viscce);

  cs_real_t *div_fext = nullptr;
  CS_MALLOC_HD(div_fext, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_divergence(m, 1, iflux, bflux, div_fext);

  cs_real_t residual = 0;
  const cs_real_t rnorm = sqrt(cs_gdot(n_cells, div_fext, div_fext));

  /* Loops on non-orthogonalities (resolution)
     ----------------------------------------- */

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
                           bc_coeffs_hp,
                           i_visc,
                           b_visc,
                           viscce,
                           rhs);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rhs[c_id] = -div_fext[c_id] - rhs[c_id];
      dphi[c_id] = 0.;
    });

    ctx.wait();

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

    cs_sles_solve_native(f->id,
                         nullptr,
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

    //cs_axpy(n_cells, 1.0, dphi, cvar_hydro_pres);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cvar_hydro_pres[c_id] += dphi[c_id];
    });
    ctx.wait();

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

  cs_sles_free_native(f->id, nullptr);

  CS_FREE_HD(viscce);
  CS_FREE_HD(rovsdt);
  CS_FREE_HD(div_fext);
  CS_FREE_HD(next_fext);
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
 * \param[in]       iterns        Navier-Stokes iteration number
 * \param[in]       nfbpcd        number of faces with condensation source term
 * \param[in]       ncmast        number of cells with condensation source terms
 * \param[in]       ifbpcd        index of faces with condensation source term
 * \param[in]       ltmast        list of cells with condensation source terms
 *                                (1 to n numbering)
 * \param[in]       isostd        indicator of standard outlet and index
 *                                of the reference outlet face
 * \param[in]       vel           velocity
 * \param[in, out]  da_uu         velocity matrix
 * \param[in]       bc_coeffs_v   boundary condition structure for the variable
 * \param[in]       bc_coeffs_dp  boundary conditions structure for the
 *                                pressure increment
 * \param[in]       spcond        variable value associated to the condensation
 *                                source term (for ivar=ipr, spcond is the
 *                                flow rate
 *                                \f$ \Gamma_{s,cond}^n \f$)
 * \param[in]       svcond        variable value associated to the condensation
 *                                source term (for ivar=ipr, svcond is the flow rate
 *                                \f$ \Gamma_{v, cond}^n \f$)
 * \param[in]       frcxt         external forces making hydrostatic pressure
 * \param[in]       dfrcxt        variation of the external forces
 *                                composing the hydrostatic pressure
 * \param[in]       i_visc        visc*surface/dist aux faces internes
 * \param[in]       b_visc        visc*surface/dist aux faces de bord
 */
/*----------------------------------------------------------------------------*/

static void
_pressure_correction_fv(int                   iterns,
                        cs_lnum_t             nfbpcd,
                        cs_lnum_t             ncmast,
                        cs_lnum_t             ifbpcd[],
                        cs_lnum_t             ltmast[],
                        const int             isostd[],
                        cs_real_t             vel[][3],
                        cs_real_t             da_uu[][6],
                        cs_field_bc_coeffs_t *bc_coeffs_v,
                        cs_field_bc_coeffs_t *bc_coeffs_dp,
                        cs_real_t             spcond[],
                        cs_real_t             svcond[],
                        cs_real_t             frcxt[][3],
                        cs_real_t             dfrcxt[][3],
                        cs_real_t             i_visc[],
                        cs_real_t             b_visc[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *)fvq->cell_cen;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *)fvq->b_face_cog;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_real_t *restrict b_dist = (const cs_real_t *)fvq->b_dist;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_t *restrict i_f_face_surf = fvq->i_f_face_surf;
  const cs_real_t *restrict b_f_face_surf = fvq->b_f_face_surf;
  const cs_real_t *restrict cell_f_vol = fvq->cell_f_vol;

  const int *bc_type = cs_glob_bc_type;

  const int meteo_profile = cs_glob_atmo_option->meteo_profile;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx, ctx_c;
#if defined(HAVE_CUDA)
  ctx_c.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  cs_field_t *f_iddp = cs_field_by_name("pressure_increment");
  cs_real_t *phi = f_iddp->val;

  /* Boundary conditions for delta P */
  cs_real_t *coefa_dp = bc_coeffs_dp->a;
  cs_real_t *coefb_dp = bc_coeffs_dp->b;
  cs_real_t *coefaf_dp = bc_coeffs_dp->af;
  cs_real_t *coefbf_dp = bc_coeffs_dp->bf;

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

  cs_real_t *cvar_hydro_pres = nullptr, *cvar_hydro_pres_prev = nullptr;
  cs_real_t *c_visc = nullptr;
  cs_real_6_t *vitenp = nullptr;
  cs_real_t *taui = nullptr, *taub = nullptr;

  cs_field_t  *f_hp = cs_field_by_name_try("hydrostatic_pressure");
  if (f_hp != nullptr) {
    cvar_hydro_pres = f_hp->vals[0];
    cvar_hydro_pres_prev = f_hp->vals[1];
  }

  /* Allocate temporary arrays */

  cs_real_t *dam, *xam, *rhs, *res;
  CS_MALLOC_HD(dam, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(xam, m->n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(res, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(rhs, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_t *phia, *iflux, *bflux, *dphi;
  CS_MALLOC_HD(phia, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(iflux, m->n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(bflux, m->n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(dphi, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_3_t *wrk;
  CS_MALLOC_HD(wrk, n_cells_ext, cs_real_3_t, cs_alloc_mode);
  cs_real_3_t *wrk2;
  CS_MALLOC_HD(wrk2, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  cs_real_t *adxk = nullptr, *adxkm1 = nullptr, *dphim1 = nullptr, *rhs0 = nullptr;
  if (eqp_p->iswdyn > 0) {
    CS_MALLOC_HD(adxk, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(adxkm1, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(dphim1, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(rhs0, n_cells_ext, cs_real_t, cs_alloc_mode);
  }

  /* Associate pointers to pressure diffusion coefficients */
  c_visc = dt;
  if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
    vitenp = (cs_real_6_t *)(cs_field_by_name("dttens")->val);

  /* Index of the field */
  const int ksinfo = cs_field_key_id("solving_info");
  cs_solving_info_t *sinfo
    = static_cast<cs_solving_info_t*>(cs_field_get_key_struct_ptr(f_p, ksinfo));

  cs_field_t *f_weight = nullptr;
  if (eqp_p->iwgrec == 1) {
    /* Weighting field for gradient */
    int kwgrec = cs_field_key_id("gradient_weighting_id");
    f_weight = cs_field_by_id(cs_field_get_key_int(f_p, kwgrec));
  }

  cs_real_t *cpro_divu = nullptr, *_cpro_divu = nullptr;
  cs_field_t *f_divu
    = cs_field_by_name_try("algo:predicted_velocity_divergence");
  if (f_divu != nullptr)
    cpro_divu = f_divu->val;
  else {
    CS_MALLOC_HD(_cpro_divu, n_cells_ext, cs_real_t, cs_alloc_mode);
    cpro_divu = _cpro_divu;
  }

  /* Boundary conditions */

  const cs_field_bc_coeffs_t *bc_coeffs_p = f_p->bc_coeffs;
  cs_real_t *coefa_p = bc_coeffs_p->a;
  cs_real_t *coefb_p = bc_coeffs_p->b;
  cs_real_t *coefaf_p = bc_coeffs_p->af;
  cs_real_t *coefbf_p = bc_coeffs_p->bf;

  /* Physical quantities */

  cs_real_t *crom = nullptr;
  const cs_real_t *brom = nullptr;

  cs_real_t *crom_eos = CS_F_(rho)->val;
  const cs_real_t *croma = nullptr;
  if (vp_param->icalhy == 1 || idilat > 1 || fluid_props->irovar) {
    croma = CS_F_(rho)->val_pre;
  }
  const cs_real_t *brom_eos = CS_F_(rho_b)->val;
  const cs_real_t *broma = nullptr;
  if (fluid_props->irovar)
    broma = CS_F_(rho_b)->val_pre;

  /* Time-interpolated density */

  cs_real_t *cpro_rho_tc = nullptr, *bpro_rho_tc = nullptr;

  if (fluid_props->irovar && (   idilat > 1
                              || vof_parameters->vof_model > 0
                              || compressible_flag == 3)) {

    cs_real_t *cpro_rho_mass = cs_field_by_name("density_mass")->val;
    cs_real_t *bpro_rho_mass = cs_field_by_name("boundary_density_mass")->val;

    /* Staggered in time velocity and pressure */
    if (eqp_u->theta < 1 && iterns > 1 && vp_param->itpcol == 0) {
      CS_MALLOC_HD(cpro_rho_tc, n_cells_ext, cs_real_t, cs_alloc_mode);
      CS_MALLOC_HD(bpro_rho_tc, m->n_b_faces, cs_real_t, cs_alloc_mode);

      const cs_real_t theta =  eqp_u->theta;
      const cs_real_t one_m_theta = 1.0 - theta;
      crom = cpro_rho_tc;

      ctx_c.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cpro_rho_tc[c_id] =  theta * cpro_rho_mass[c_id]
                            + one_m_theta * croma[c_id];
      });

      ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
        bpro_rho_tc[f_id] =   theta * bpro_rho_mass[f_id]
                            + one_m_theta * broma[f_id];
      });

      ctx.wait();

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
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t j = 0; j < 6; j++)
        da_uu[c_id][j] = vitenp[c_id][j];
    });
    ctx.wait();

    cs_mesh_sync_var_sym_tens(da_uu);
  }

  /* Calculation of dt/rho */

  cs_real_t *xdtsro = nullptr;
  cs_real_6_t *tpusro = nullptr;

  if (vof_parameters->vof_model > 0 || idilat == 4) {

    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION) {
      CS_MALLOC_HD(xdtsro, n_cells_ext, cs_real_t, cs_alloc_mode);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        xdtsro[c_id] = dt[c_id]/crom[c_id];
      });
      ctx.wait();

      cs_halo_sync_var(m->halo, CS_HALO_STANDARD, xdtsro);
    }
    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) {
      CS_MALLOC_HD(tpusro, n_cells_ext, cs_real_6_t, cs_alloc_mode);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        const cs_real_t drom = 1. / crom[c_id];
        for (cs_lnum_t j = 0; j < 6; j++)
          tpusro[c_id][j] = vitenp[c_id][j] * drom;
      });
      ctx.wait();

      cs_mesh_sync_var_sym_tens(tpusro);

      vitenp = tpusro;
    }

    /* Associate pointers to pressure diffusion coefficient */

    c_visc = xdtsro;

  }

  if (vp_param->staggered == 1) {

    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, dt);

    const cs_real_t *hli = cs_field_by_name("inner_face_head_loss")->val;
    const cs_real_t *hlb = cs_field_by_name("boundary_face_head_loss")->val;

    CS_MALLOC_HD(taui, n_i_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(taub, n_b_faces, cs_real_t, cs_alloc_mode);

    ctx_c.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      cs_lnum_t c_id_0 = i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = i_face_cells[f_id][1];
      cs_real_t dt_f = 0.5 * (dt[c_id_0] + dt[c_id_1]);

      taui[f_id] = dt_f / (1. + hli[f_id] * dt_f);
    });

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      cs_lnum_t c_id = b_face_cells[f_id];
      taub[f_id] = dt[c_id] / (1. + hlb[f_id] * dt[c_id]);
    });

  }

  /* Compute an approximated pressure increment if needed,
   * that is when there are buoyancy terms (gravity and variable density)
   * with a free outlet.
   * ==================================================================== */

  /* Standard initialization */

  cs_arrays_set_value<cs_real_t, 1>(n_i_faces, 0., iflux);

  if (vp_param->staggered == 0) {
    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      coefa_dp[f_id] = 0.;
      coefaf_dp[f_id] = 0.;
      coefb_dp[f_id] = coefb_p[f_id];
      coefbf_dp[f_id] = coefbf_p[f_id];
      bflux[f_id] = 0.;
    });
  }
  else {
    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      coefa_dp[f_id] = coefa_p[f_id];
      coefaf_dp[f_id] = coefaf_p[f_id];
      coefb_dp[f_id] = coefb_p[f_id];
      coefbf_dp[f_id] = coefbf_p[f_id];
      bflux[f_id] = 0.;
    });
  }

  ctx.wait();
  ctx_c.wait();

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

    if (f_hp != nullptr && indhyd == 1) {

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
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cvar_hydro_pres[c_id] -= phydr0;
      });

    }

    /* If hydrostatic pressure increment or free entrance Inlet. */

    if (indhyd == 1 || vp_param->iifren == 1) {

      const int *auto_flag = cs_glob_bc_pm_info->iautom;
      const cs_real_t *b_head_loss
        = cs_boundary_conditions_get_b_head_loss(false);

      ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {

        int iautof = 0;

        /* automatic inlet/outlet face for atmospheric flow */
        if (meteo_profile > 0)
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
               NB: eps =1.d-1 must be consistent
               with `cs_face_anisotropic_viscosity_scalar`. */

            fikis = cs_math_fmax(fikis, 1.e-1*sqrt(viscis)*b_dist[f_id]);

            hint = viscis/b_face_surf[f_id]/fikis;

          }

          if (indhyd == 1) {
            if (f_hp != nullptr) {
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

              /* Convective_outlet BC */

              // Gradient BCs
              coefb_dp[f_id] = cfl / (1.0 + cfl);
              coefa_dp[f_id] = (1.0 - coefb_dp[f_id]) * pimp;

              // Flux BCs
              coefaf_dp[f_id] = - hint * coefa_dp[f_id];
              coefbf_dp[f_id] =   hint * (1.0 - coefb_dp[f_id]);
            }

            else
              coefaf_dp[f_id] = - hint*coefa_dp[f_id];

          }

          /* Outher boundary face types */

          else
            coefaf_dp[f_id] = - hint*coefa_dp[f_id];

        } /* if (isostd[f_id] == 1 || (open_bcs_flag >= 1 && iautof >= 1)) */

      }); /* End of loop on boundary faces */

    }  /* if (indhyd == 1 || vp_param->iifren == 1) */

  } /* if (vp_param->iphydr == 1 && vp_param->iifren == 1) */

  /* Building the linear system to solve.
   * ==================================== */

  /* Implicit term */

  cs_real_t *rovsdt;
  CS_MALLOC_HD(rovsdt, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_arrays_set_value<cs_real_t, 1>(n_cells_ext, 0., rovsdt);

  /* Compressible scheme implicit part;
     Getting the thermal parameters */

  int ieos = cs_glob_cf_model->ieos;
  int thermal_variable = cs_glob_thermal_model->thermal_variable;
  int kinetic_st = cs_glob_thermal_model->has_kinetic_st;

  const cs_real_t *temp = nullptr;
  cs_real_t *xcpp = nullptr, *dc2 = nullptr;
  cs_real_t *cvar_th = nullptr, *tempk = nullptr;
  cs_real_t _coef = 0;
  cs_real_t *yw = nullptr, *yv = nullptr;

  if (idilat == 2) {

    /* Get the temperature */
    if (thermal_variable == CS_THERMAL_MODEL_TEMPERATURE)
      temp = CS_F_(t)->val;
    else {
      const cs_field_t *f_t = cs_field_by_name_try("temperature");
      if (f_t != nullptr) {
        tempk = f_t->val;
        temp = f_t->val;
      }
    }

    if (temp != nullptr) {

      cvar_th = CS_F_(t)->val;

      /* Allocation */
      CS_MALLOC_HD(dc2, n_cells_ext, cs_real_t, cs_alloc_mode);
      CS_MALLOC_HD(xcpp, n_cells_ext, cs_real_t, cs_alloc_mode);

      /* Theta scheme related term */
      _coef = 1. + 2. * (1. - eqp_u->theta);

      /* Get cp */
      if (fluid_props->icp > 0) {
        cs_real_t *cpro_cp = CS_F_(cp)->val;
        ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          xcpp[c_id] = cpro_cp[c_id];
        });
      }
      else {
        ctx_c.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          xcpp[c_id] = 1.;
        });
      }

      ctx_c.wait();

      /* Get mass fractions if needed */
      if (ieos == CS_EOS_MOIST_AIR) {
        yw = cs_field_by_name("yw")->val;
        yv = cs_field_by_name("yv")->val;
      }
      cs_real_t *cvar_fracm = nullptr;

      /* Compute dc2 */
      cs_thermal_model_c_square(xcpp,
                                temp,
                                cvar_pr,
                                yv,
                                cvar_fracm,
                                yw,
                                dc2);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        rovsdt[c_id] += cell_f_vol[c_id] * _coef * dc2[c_id] / dt[c_id];
      });

      BFT_FREE(xcpp);

    }
  }

  /* Implicit part of the cavitation source */
  if (i_vof_mass_transfer != 0 && cavitation_parameters->itscvi == 1) {
    cs_real_t *dgdpca = cs_field_by_name("model:cavitation_st_dgdp")->val;
    cs_real_t rho1 = vof_parameters->rho1;
    cs_real_t rho2 = vof_parameters->rho2;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rovsdt[c_id] -= cell_f_vol[c_id] * dgdpca[c_id] * (1./rho2 - 1./rho1);
    });
  }

  /* Strengthen the diagonal for Low Mach Algorithm */
  if (idilat == 3) {
    const cs_real_t epsdp = vp_param->epsdp;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rovsdt[c_id] += epsdp*cell_f_vol[c_id]/dt[c_id];
    });
  }

  cs_real_t *c2 = nullptr;

  if (compressible_flag == 3) {

    cs_real_t *cvar_fracv = nullptr;
    cs_real_t *cvar_fracm = nullptr;
    cs_real_t *cvar_frace = nullptr;

    cs_real_t *cpro_cp = nullptr;
    cs_real_t *cpro_cv = nullptr;

    if (fluid_props->icp >= 0)
      cpro_cp = cs_field_by_id(fluid_props->icp)->val;

    if (fluid_props->icv >= 0)
      cpro_cv = cs_field_by_id(fluid_props->icv)->val;

    CS_MALLOC_HD(c2, n_cells_ext, cs_real_t, cs_alloc_mode);

    cs_cf_thermo_c_square(cpro_cp, cpro_cv, cvar_pr, crom,
                          cvar_fracv, cvar_fracm, cvar_frace, c2, n_cells);

    int istat = eqp_p->istat;
    if (istat != 0) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        rovsdt[c_id] += istat*(cell_f_vol[c_id]/(dt[c_id]*c2[c_id]));
      });
    }

  }

  /* Face diffusivity */

  cs_real_t *weighb = nullptr;
  cs_real_2_t *weighf = nullptr;

  if (eqp_p->idiff >= 1) {

    /* Scalar diffusivity */
    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION) {

      cs_face_viscosity(m,
                        fvq,
                        eqp_p->imvisf,
                        c_visc,
                        i_visc,
                        b_visc);

      if (f_weight != nullptr) {  /* Weighting for gradient */
        cs_real_t *weight = f_weight->val;

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          weight[c_id] = c_visc[c_id];
        });
        ctx.wait();

        cs_halo_sync_var(m->halo, CS_HALO_STANDARD, weight);
      }

    }

    /* Tensor diffusivity */
    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) {

      /* Allocate temporary arrays */
      CS_MALLOC_HD(weighf, n_i_faces, cs_real_2_t, cs_alloc_mode);
      CS_MALLOC_HD(weighb, n_b_faces, cs_real_t, cs_alloc_mode);

      cs_face_anisotropic_viscosity_scalar(m,
                                           fvq,
                                           vitenp,
                                           eqp_p->verbosity,
                                           weighf,
                                           weighb,
                                           i_visc,
                                           b_visc);

      if (f_weight != nullptr) { /* Weighting for gradient */
        cs_real_6_t *weight = (cs_real_6_t *)f_weight->val;

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          for (cs_lnum_t j = 0; j < 6; j++)
            weight[c_id][j] = vitenp[c_id][j];
        });
        ctx.wait();

        cs_mesh_sync_var_sym_tens((cs_real_6_t *)(f_weight->val));

      }

    }

  }
  else {
    ctx_c.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      i_visc[f_id] = 0.;
    });

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      b_visc[f_id] = 0.;
    });
  }

  if (vp_param->staggered == 1) {
    ctx_c.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      cs_lnum_t c_id_0 = i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = i_face_cells[f_id][1];
      cs_real_t dtm = 0.5 * (dt[c_id_0] + dt[c_id_1]);
      i_visc[f_id] = taui[f_id] / dtm * i_visc[f_id];
    });

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      cs_lnum_t c_id = b_face_cells[f_id];
      b_visc[f_id] = taub[f_id] / dt[c_id] * b_visc[f_id];
    });
  }

  ctx.wait(); // needed for rovsdt and b_visc
  ctx_c.wait(); // needed for i_visc

  cs_matrix_wrapper_scalar(eqp_p->iconv, eqp_p->idiff, eqp_p->ndircl,
                           isym,
                           1.0, /* thetap */
                           0,   /* imucpp */
                           bc_coeffs_dp, rovsdt,
                           imasfl, bmasfl, i_visc, b_visc,
                           nullptr, dam, xam);

  /* Mass flux initialization
     ------------------------ */

  /* Predicted mass flux and first Rhie and Chow component */

  cs_real_3_t  *gradp;
  CS_MALLOC_HD(gradp, n_cells_ext, cs_real_3_t, cs_alloc_mode);

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

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          for (cs_lnum_t j = 0; j < 3; j++) {
            wrk[c_id][j] =   gradp[c_id][j] - frcxt[c_id][j]
                           - cpro_poro_div_duq[c_id][j];
          }
        });
      }

      else {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          for (cs_lnum_t j = 0; j < 3; j++) {
            wrk[c_id][j] = gradp[c_id][j] - frcxt[c_id][j];
          }
        });
      }

    }
    else {
      cs_array_copy<cs_real_t>(3*n_cells,
                               (const cs_real_t *)gradp,
                               (cs_real_t *)wrk);
    }

    if (   (eqp_p->idften & CS_ISOTROPIC_DIFFUSION)
        && vp_param->rcfact == 0) {

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t arsrdt = (arak/crom[c_id]) * dt[c_id];

        for (cs_lnum_t j = 0; j < 3; j++)
          wrk[c_id][j] *= arsrdt;
      });

    }

    else if (   (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
             ||  vp_param->rcfact == 1) {

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
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
      });

    }

  }

  else {
    cs_arrays_set_value<cs_real_t, 1>(3*n_cells, 0., (cs_real_t *)wrk);
  }

  if (idilat < 4) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        wrk[c_id][j] += vel[c_id][j];
        wrk2[c_id][j] = wrk[c_id][j];
      }
    });
  }

  ctx.wait();

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
                 static_cast<cs_gradient_limit_t>(eqp_u->imligr),
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 wrk,
                 bc_coeffs_v,
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
      ctx_c.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
        imasfla[f_id] = imasfl[f_id];
      });

      ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
        bmasfla[f_id] = bmasfl[f_id];
      });

    }

    if (cs_glob_porous_model >= 1) {

      const cs_real_t *c_porosity = cs_field_by_name("porosity")->val;

      ctx_c.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
        cs_lnum_t c_id_0 = i_face_cells[f_id][0];
        cs_lnum_t c_id_1 = i_face_cells[f_id][1];
        cs_real_t dtm = 0.5 * (dt[c_id_0]+dt[c_id_1]);
        cs_real_t porosf = cs_math_fmin(c_porosity[c_id_0], c_porosity[c_id_1]);
        imasfl[f_id] =   taui[f_id] / dtm
                       * imasfla[f_id]+porosf*taui[f_id]*sti[f_id]
                       * i_f_face_surf[f_id];
      });

    }

    else { /* cs_glob_porous_model == 0) */

      ctx_c.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
        cs_lnum_t c_id_0 = i_face_cells[f_id][0];
        cs_lnum_t c_id_1 = i_face_cells[f_id][1];
        cs_real_t dtm = 0.5 * (dt[c_id_0] + dt[c_id_1]);
        imasfl[f_id] =   taui[f_id] / dtm
                       * imasfla[f_id]+taui[f_id]*sti[f_id]
                       * i_f_face_surf[f_id];
      });

    } /* end of test on cs_glob_porous_model */

    ctx_c.wait();

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      cs_lnum_t c_id = b_face_cells[f_id];

      if (bc_type[f_id] == CS_INLET) {

        bmasfl[f_id] = taub[f_id] / dt[c_id] * bmasfl[f_id];

        cs_real_t dimp =   -(1. - dt[c_id]/taub[f_id])
                         * bmasfl[f_id]/b_f_face_surf[f_id];
        cs_real_t hint = taub[f_id] / b_dist[f_id];

        /* Neumann_scalar BC */

        // Gradient BCs
        coefa_dp[f_id] = -dimp/cs_math_fmax(hint, 1.e-300);
        coefb_dp[f_id] = 1.;

        // Flux BCs
        coefaf_dp[f_id] = dimp;
        coefbf_dp[f_id] = 0.;
      }
      else
        bmasfl[f_id] = taub[f_id] / dt[c_id] * bmasfla[f_id];
    });

    ctx.wait();

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

    cs_real_t *ipro_visc = nullptr, *bpro_visc = nullptr;

    CS_MALLOC_HD(ipro_visc, n_i_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(bpro_visc, n_b_faces, cs_real_t, cs_alloc_mode);

    /* Scalar diffusivity */
    if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION && vp_param->rcfact == 0) {

      cs_real_t *cpro_visc;
      CS_MALLOC_HD(cpro_visc, n_cells_ext, cs_real_t, cs_alloc_mode);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cpro_visc[c_id] = arak * c_visc[c_id];
      });
      ctx.wait();

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
        ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
          if (bc_type[f_id] == CS_COUPLED_FD)
            bpro_visc[f_id] = 0.;
        });
        ctx.wait();
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
                                  bc_coeffs_p,
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

      CS_FREE_HD(cpro_visc);

    }

    /* Tensor diffusivity */
    else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION || vp_param->rcfact == 1) {

      cs_real_6_t *cpro_vitenp;
      CS_MALLOC_HD(cpro_vitenp, n_cells_ext, cs_real_6_t, cs_alloc_mode);

      if (idilat == 4 || vof_parameters->vof_model > 0) {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          cs_real_t arsr = arak / crom[c_id];
          for (cs_lnum_t j = 0; j < 6; j++)
            cpro_vitenp[c_id][j] = arsr * da_uu[c_id][j];
        });
      }
      else {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          for (cs_lnum_t j = 0; j < 6; j++)
            cpro_vitenp[c_id][j] = arak * da_uu[c_id][j];
        });
      }

      ctx.wait();

      cs_real_2_t *weighftp = nullptr;
      cs_real_t *weighbtp = nullptr;
      CS_MALLOC_HD(weighftp, n_i_faces, cs_real_2_t, cs_alloc_mode);
      CS_MALLOC_HD(weighbtp, n_b_faces, cs_real_t, cs_alloc_mode);

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
        ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
          if (bc_type[f_id] == CS_COUPLED_FD)
            bpro_visc[f_id] = 0.;
        });
        ctx.wait();
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
                                              bc_coeffs_p,
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

      CS_FREE_HD(cpro_vitenp);
      CS_FREE_HD(weighftp);
      CS_FREE_HD(weighbtp);
    }

    CS_FREE_HD(ipro_visc);
    CS_FREE_HD(bpro_visc);
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
    if (f_hp != nullptr) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        phi[c_id]  = cvar_hydro_pres[c_id] - cvar_hydro_pres_prev[c_id];
        dphi[c_id] = 0.;
        phia[c_id] = phi[c_id];
      });
    }
    else {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        phi[c_id] = 0.;
        dphi[c_id] = 0.;
        phia[c_id] = phi[c_id];
      });
    }
  }
  else {
    if (f_hp != nullptr) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        phi[c_id]  = cvara_pr[c_id] + cvar_hydro_pres[c_id]
                                    - cvar_hydro_pres_prev[c_id];
        dphi[c_id] = 0.;
        phia[c_id] = phi[c_id];
      });
    }
    else {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        phi[c_id]  = cvara_pr[c_id];
        dphi[c_id] = 0.;
        phia[c_id] = phi[c_id];
      });
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

  cs_real_t *velflx = nullptr, *velflb = nullptr;

  if (idilat >= 4) {

    cs_real_t *cpro_tsrho = cs_field_by_name("dila_st")->val;

    CS_MALLOC_HD(velflx, n_i_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(velflb, n_b_faces, cs_real_t, cs_alloc_mode);

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
                 static_cast<cs_gradient_limit_t>(eqp_u->imligr),
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 vel,
                 bc_coeffs_v,
                 velflx, velflb);

    cs_divergence(m, init, velflx, velflb, res);

    if (idilat == 4) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cpro_divu[c_id] +=   res[c_id]
        /* Add the dilatation source term D(rho)/Dt */
                           + cpro_tsrho[c_id] / crom[c_id];
      });
    }
    else {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cpro_divu[c_id] +=   res[c_id]*crom[c_id]
        /* Add the dilatation source term D(rho)/Dt */
                           + cpro_tsrho[c_id];
      });
    }

    /* The mass flux is completed by u*.S (idilat=4)
     *                               rho u* . S (idilat=5)
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
                 static_cast<cs_gradient_limit_t>(eqp_u->imligr),
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 vel,
                 bc_coeffs_v,
                 imasfl, bmasfl);

  }

  /* Mass source terms adding for volumic flow rate */

  cs_lnum_t n_elts = 0;
  const cs_lnum_t *elt_ids = nullptr;
  cs_real_t *mst_val_p = nullptr;

  cs_volume_mass_injection_get_arrays(f_p, &n_elts, &elt_ids, nullptr,
                                      &mst_val_p, nullptr);

  if (n_elts > 0) {
    ctx.parallel_for(n_elts, [=] CS_F_HOST_DEVICE (cs_lnum_t c_idx) {
      cs_lnum_t c_id = elt_ids[c_idx];
      cpro_divu[c_id] -= cell_f_vol[c_id] * mst_val_p[c_idx];
    });
  }

  /* Source term adding for condensation modelling */

  if (nfbpcd > 0) {
    const int var_id_key = cs_field_key_id("variable_id");
    const int ipr = cs_field_get_key_int(f_p, var_id_key);

    cs_real_t *_spcond = spcond + (ipr-1)*nfbpcd;

    ctx.parallel_for(nfbpcd, [=] CS_F_HOST_DEVICE (cs_lnum_t f_idx) {
      cs_lnum_t f_id = ifbpcd[f_idx];
      cs_lnum_t c_id = b_face_cells[f_id];
      cpro_divu[c_id] -= b_face_surf[f_id] * _spcond[f_idx];
    });
  }

  /* volume Gamma source for metal mass structures
     condensation modelling */

  if (ncmast > 0) {
    const int var_id_key = cs_field_key_id("variable_id");
    const int ipr = cs_field_get_key_int(f_p, var_id_key);

    cs_real_t *_svcond = svcond + (ipr-1)*ncmast;
    cs_real_t *surfbm = nullptr;
    CS_MALLOC_HD(surfbm, ncmast, cs_real_t, cs_alloc_mode);

    cs_wall_condensation_volume_exchange_surf_at_cells(surfbm);

    ctx.parallel_for(ncmast, [=] CS_F_HOST_DEVICE (cs_lnum_t c_idx) {
      cs_lnum_t c_id = ltmast[c_idx];
      cpro_divu[c_id] -= surfbm[c_idx] * _svcond[c_idx];
    });

    ctx.wait();

    CS_FREE_HD(surfbm);
  }

  /* Source term associated to the mass aggregation */

  if ((idilat == 2 || idilat == 3) && compressible_flag != 3) {
    if (ieos == CS_EOS_NONE) { // If no particular EOS is set
      if (vp_param->itpcol == 1 && eqp_u->theta < 1.) {
        cs_real_t *imasfla
          = cs_field_by_id(cs_field_get_key_int(f_p, kimasf))->val_pre;
        cs_real_t *bmasfla
          = cs_field_by_id(cs_field_get_key_int(f_p, kbmasf))->val_pre;
        cs_real_t *divu_prev;
        CS_MALLOC_HD(divu_prev, n_cells_ext, cs_real_t, cs_alloc_mode);

        cs_divergence(m, 1, imasfla, bmasfla, divu_prev);

        cs_real_t theta = eqp_u->theta;

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          cs_real_t drom = crom_eos[c_id] - croma[c_id];
          cpro_divu[c_id] +=  (1. + theta) *drom *cell_f_vol[c_id] /dt[c_id]
                             + theta *divu_prev[c_id];
        });

        ctx.wait();

        CS_FREE_HD(divu_prev);
      }
      else {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          cs_real_t drom = crom_eos[c_id] - croma[c_id];
          cpro_divu[c_id] += drom *cell_f_vol[c_id] /dt[c_id];
        });
      }
    }
    else { /* compressible scheme explicit part */
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t drom = crom_eos[c_id] - croma[c_id];
        cs_real_t drop =   (-1 + _coef ) * (cvar_pr[c_id] - cvara_pr[c_id])
                         * dc2[c_id];
        cpro_divu[c_id] += (drom + drop) * cell_f_vol[c_id] / dt[c_id];
      });
    }

  }

  /* Lagrangian source terms */

  if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
      && cs_glob_lagr_source_terms->ltsmas == 1) {

    cs_real_t *lag_st_m = cs_field_by_name("lagr_st_pressure")->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cpro_divu[c_id] -= lag_st_m[c_id];
    });
  }

  /* Cavitation source term */

  if (i_vof_mass_transfer != 0) {
    cs_real_t *gamcav = cs_field_by_name("model:cavitation_gamma")->val;
    cs_real_t rho1 = vof_parameters->rho1;
    cs_real_t rho2 = vof_parameters->rho2;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cpro_divu[c_id] -= cell_f_vol[c_id]*gamcav[c_id]*(1./rho2 - 1./rho1);
    });
  }

  /* Norm residual
   * Historical norm for the pressure step:
   *   div(rho u* + dt gradP^(n))-Gamma
   *   i.e.  RHS of the pressure + div(dt gradP^n) (otherwise there is a risk
   *   a 0 norm at steady states...). Represents terms that pressure has to
   *   balance. */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    const cs_real_t dt_d_rho = dt[c_id] / crom[c_id];

    for (cs_lnum_t j = 0; j < 3; j++)
      wrk[c_id][j] = dt_d_rho * gradp[c_id][j];

  });

  ctx.wait();

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
                 static_cast<cs_gradient_limit_t>(eqp_u->imligr),
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 wrk,
                 bc_coeffs_v,
                 iflux, bflux);

    cs_divergence(m, 1, iflux, bflux, res);
  }

  CS_FREE_HD(iflux);
  CS_FREE_HD(bflux);

  if (idilat >= 4) {
    /* Weakly compressible algorithm: semi analytic scheme */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      res[c_id] =   res[c_id]*crom[c_id]
      /* It is: div(dt/rho*rho grad P) + div(rho u*) - Gamma
         NB: if iphydr=1, div(rho u*) contains div(dt d fext).  */
                  + cpro_divu[c_id];
    });
  }
  else {
    /* It is: div(dt/rho*rho grad P) + div(rho u*) - Gamma
       NB: if iphydr=1, div(rho u*) contains div(dt d fext).  */
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      res[c_id] += cpro_divu[c_id];
    });
  }

  ctx.wait();

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
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      adxkm1[c_id] = 0.;
      adxk[c_id] = 0.;
    });

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
                             bc_coeffs_dp,
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
                                         bc_coeffs_dp,
                                         i_visc, b_visc,
                                         vitenp,
                                         weighf, weighb,
                                         rhs);
  }

  if (eqp_p->iswdyn >= 1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* Dynamic relaxation: stores the initial rhs */
      rhs0[c_id] = rhs[c_id];

      /* Finalize the rhs initialization */
      rhs[c_id] = -rhs[c_id] - cpro_divu[c_id] - rovsdt[c_id]*phi[c_id];
    });
  }
  else {
    /* Finalize the rhs initialization */
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rhs[c_id] = -rhs[c_id] - cpro_divu[c_id] - rovsdt[c_id]*phi[c_id];
    });
  }

  ctx.wait();

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
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        dphim1[c_id] = dphi[c_id];
        dphi[c_id] = 0.;
      });
    }
    else {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        dphi[c_id] = 0.;
      });
    }

    ctx.wait();

    cs_real_t ressol = residual;   /* solver residual */

    cs_sles_solve_native(f_p->id, nullptr,
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

      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        adxkm1[c_id] = adxk[c_id];
        adxk[c_id] = - rhs0[c_id];
      });

      ctx.wait();

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
                               bc_coeffs_dp,
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
                                           bc_coeffs_dp,
                                           i_visc, b_visc,
                                           vitenp,
                                           weighf, weighb,
                                           adxk);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        adxk[c_id] = - adxk[c_id];
      });

      ctx.wait();

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
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          phia[c_id] = phi[c_id];
          phi[c_id] = phi[c_id] + eqp_p->relaxv*dphi[c_id];
        });
      }
      /* If it is the last sweep, update with the total increment */
      else {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          phia[c_id] = phi[c_id];
          phi[c_id] = phi[c_id] + dphi[c_id];
        });
      }
    }
    else if (eqp_p->iswdyn == 1) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        phia[c_id] = phi[c_id];
        phi[c_id] = phi[c_id] + alph*dphi[c_id];
      });
    }
    else if (eqp_p->iswdyn >= 2) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        phia[c_id] = phi[c_id];
        phi[c_id] = phi[c_id] + alph*dphi[c_id] + beta*dphim1[c_id];
      });
    }

    ctx.wait();

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
                               bc_coeffs_dp,
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
                                           bc_coeffs_dp,
                                           i_visc, b_visc,
                                           vitenp,
                                           weighf, weighb,
                                           rhs);
    }

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rhs[c_id] = - cpro_divu[c_id] - rhs[c_id] - rovsdt[c_id]*phi[c_id];
    });

    ctx.wait();

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
  if (cs_math_fabs(rnormp) > 0.)
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
  if (f_err_est != nullptr) {
    cs_real_t *c_estim_der = f_err_est->val;
    const cs_real_t *restrict cell_vol = fvq->cell_vol;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      c_estim_der[c_id] = cs_math_fabs(rhs[c_id]) / cell_vol[c_id];
    });
  }
  f_err_est = cs_field_by_name_try("est_error_der_2");
  if (f_err_est != nullptr) {
    cs_real_t *c_estim_der = f_err_est->val;
    const cs_real_t *restrict cell_vol = fvq->cell_vol;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      c_estim_der[c_id] = cs_math_fabs(rhs[c_id]) / sqrt(cell_vol[c_id]);
    });
  }

  /* Update the mass flux
     -------------------- */

  /* We cancel the face viscosity for coupled faces so as not to modify the
     boundary mass flux in the case of a pressure Dirichlet:
     pressure correction and filter are canceled. */

  if (cs_sat_coupling_n_couplings() > 0) {
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      if (bc_type[f_id] == CS_COUPLED_FD)
        b_visc[f_id] = 0.;
    });
    ctx.wait();
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
                                  bc_coeffs_dp,
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
                                              bc_coeffs_dp,
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
                                  bc_coeffs_dp,
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
                                              bc_coeffs_dp,
                                              i_visc, b_visc,
                                              vitenp,
                                              weighf, weighb,
                                              imasfl, bmasfl);

  }

  /* Update density
     -------------- */

  if (compressible_flag == 3) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      crom_eos[c_id] += phi[c_id]/c2[c_id];
    });
    ctx.wait();
  }

  /* Suppression of the grid hierarchy (release solver setup)
     --------------------------------- */

  cs_sles_free_native(f_p->id, nullptr);

  CS_FREE_HD(dam);
  CS_FREE_HD(xam);
  CS_FREE_HD(c2);

  /* Weakly compressible algorithm: semi analytic scheme
     2nd step solving a convection diffusion equation
     =================================================== */

  if (idilat == 5) {

    cs_real_t *ddphi;
    CS_MALLOC_HD(ddphi, n_cells_ext, cs_real_t, cs_alloc_mode);

    cs_field_bc_coeffs_t bc_coeffs_loc;
    CS_MALLOC_HD(bc_coeffs_loc.a,  3*n_b_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(bc_coeffs_loc.af, 3*n_b_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(bc_coeffs_loc.b,  9*n_b_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(bc_coeffs_loc.bf, 9*n_b_faces, cs_real_t, cs_alloc_mode);

    cs_real_3_t  *coefar = (cs_real_3_t  *)bc_coeffs_loc.a;
    cs_real_3_t  *cofafr = (cs_real_3_t  *)bc_coeffs_loc.af;
    cs_real_33_t *coefbr = (cs_real_33_t *)bc_coeffs_loc.b;
    cs_real_33_t *cofbfr = (cs_real_33_t *)bc_coeffs_loc.bf;

    /* Boundary condition for the pressure increment
       coefb, coefbf are those of the pressure */

    cs_field_bc_coeffs_t bc_coeffs_dp2;
    cs_field_bc_coeffs_shallow_copy(bc_coeffs_p, &bc_coeffs_dp2);
    CS_MALLOC_HD(bc_coeffs_dp2.a,  n_b_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(bc_coeffs_dp2.af, n_b_faces, cs_real_t, cs_alloc_mode);

    cs_real_t *coefa_dp2  = bc_coeffs_dp2.a;
    cs_real_t *coefaf_dp2 = bc_coeffs_dp2.af;

    /* Convective flux: dt/rho grad(rho) */

    cs_field_bc_coeffs_t bc_coeffs_rho_loc;
    cs_field_bc_coeffs_init(&bc_coeffs_rho_loc);
    CS_MALLOC_HD(bc_coeffs_rho_loc.a, n_b_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(bc_coeffs_rho_loc.b, n_b_faces, cs_real_t, cs_alloc_mode);

    cs_real_t *coefa_rho = bc_coeffs_rho_loc.a;
    cs_real_t *coefb_rho = bc_coeffs_rho_loc.b;

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {

      /* Dirichlet Boundary Condition on rho */
      coefa_rho[f_id] = brom[f_id];
      coefb_rho[f_id] = 0.;

      coefa_dp2[f_id] = 0.;
      coefaf_dp2[f_id] = 0.;

      if (   eqp_p->idften & CS_ISOTROPIC_DIFFUSION
          || eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) { // test inutile ?

        /* Neumann boundary Conditions for the convective flux (qimpv = 0) */

        //cs_lnum_t c_id = b_face_cells[f_id];
        //cs_real_t hint = dt[c_id] / b_dist[f_id];

        // Gradient and Flux BCs (qimpv=0 --> a=af=bf=0, b=Id)

        for (cs_lnum_t i = 0; i < 3; i++) {
          coefar[f_id][i] = 0.;
          cofafr[f_id][i] = 0.;

          for (cs_lnum_t j = 0; j < 3; j++) {
            cofbfr[f_id][i][j] = 0.;
            coefbr[f_id][i][j] = 0.;
          }
          coefbr[f_id][i][i] = 1.0;
        }
      }
    });

    ctx.wait();

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
                       static_cast<cs_gradient_limit_t>(eqp_u->imligr),
                       eqp_u->epsrgr,
                       eqp_u->climgr,
                       nullptr,          /* f_ext */
                       &bc_coeffs_rho_loc,
                       crom,
                       nullptr,         /* c_weight */
                       nullptr,         /* cpl */
                       gradp);

    CS_FREE_HD(bc_coeffs_rho_loc.a);
    CS_FREE_HD(bc_coeffs_rho_loc.b);

    /* dt/rho * grad rho */
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rhs[c_id] = 0.;

      for (cs_lnum_t j = 0; j < 3; j++) {
        wrk[c_id][j] = gradp[c_id][j] * dt[c_id] / crom[c_id];
      }
    });

    ctx.wait();

    /* Viscosity */
    cs_face_viscosity(m,
                      fvq,
                      eqp_u->imvisf,
                      dt,
                      i_visc, b_visc);

    /* (dt/rho * grad rho) . S */

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
                 static_cast<cs_gradient_limit_t>(eqp_u->imligr),
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 wrk,
                 &bc_coeffs_loc,
                 velflx, velflb);

    /* Convective source term */

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
    eqp_loc.theta = 1;
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
                      &bc_coeffs_dp2,
                      velflx, velflb,
                      i_visc, b_visc,
                      nullptr,  /* viscel */
                      nullptr,  /* xcpp */
                      nullptr,  /* weighf */
                      nullptr,  /* weighb */
                      0,     /* icvflb; upwind scheme */
                      nullptr,
                      rhs);

    /* Initialization of the variable to solve */
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rovsdt[c_id] = 340.0/dt[c_id] * cell_f_vol[c_id];
      dphi[c_id]   = 0.;
      ddphi[c_id]  = 0.;
      rhs[c_id]    = - rhs[c_id];
    });

    ctx.wait();

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
                                       &bc_coeffs_dp2,
                                       velflx, velflb,
                                       i_visc, b_visc,
                                       i_visc, b_visc,
                                       nullptr,   /* viscel */
                                       weighf, weighb,
                                       0,      /* icvflb (upwind conv. flux) */
                                       nullptr,   /* icvfli */
                                       rovsdt,
                                       rhs,
                                       dphi, ddphi,
                                       nullptr,   /* xcpp */
                                       nullptr);  /* eswork */

    cs_sles_pop(f_p->id);

    /* Update the pressure increment */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      phi[c_id] += dphi[c_id];
      /* Remove the last increment */
      dphi[c_id] -= ddphi[c_id];
    });

    ctx.wait();

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
                                  &bc_coeffs_dp2,
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
                                              &bc_coeffs_dp2,
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
                                  &bc_coeffs_dp2,
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
                                              &bc_coeffs_dp2,
                                              i_visc, b_visc,
                                              vitenp,
                                              weighf, weighb,
                                              imasfl, bmasfl);

    /* Free memory */
    CS_FREE_HD(ddphi);
    CS_FREE_HD(coefar);
    CS_FREE_HD(coefbr);
    CS_FREE_HD(cofafr);
    CS_FREE_HD(cofbfr);

    cs_field_bc_coeffs_free_copy(bc_coeffs_p, &bc_coeffs_dp2);

  } /* End if weaky compressible algorithm (idilat = 5) */

  CS_FREE_HD(velflx);
  CS_FREE_HD(velflb);

  /* Update the pressure field
     ========================= */

  // Pressure at the last sub-iteration

  cs_real_t *pk1;
  CS_MALLOC_HD(pk1, n_cells_ext, cs_real_t, cs_alloc_mode);

  if (idtvar < 0) {
    const cs_real_t relaxv = eqp_p->relaxv;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      pk1[c_id] = cvar_pr[c_id];
      cvar_pr[c_id] += relaxv*phi[c_id];
    });
  }
  else {
    if (vp_param->staggered == 0) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        pk1[c_id] = cvar_pr[c_id];
        cvar_pr[c_id] += phi[c_id];
      });
    }
    else {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        pk1[c_id] = cvar_pr[c_id];
        cvar_pr[c_id] = phi[c_id];
      });
    }
  }

  ctx.wait();

  /* Transformation of volume fluxes into mass fluxes */

  /* Update the density when solving the Helmholtz equation */
  if (ieos != CS_EOS_NONE) {
    cs_real_t *cpro_rho_mass = cs_field_by_name("density_mass")->val;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t drop = (  _coef * cvar_pr[c_id] - (_coef - 1.) * cvara_pr[c_id]
                        - pk1[c_id]) * dc2[c_id];
      crom_eos[c_id] += drop;
      cpro_rho_mass[c_id] = crom_eos[c_id];
    });
    ctx.wait();
  }

  CS_FREE_HD(dc2);

  /* Compute the isobaric heat capacity if needed */
  cs_field_t  *f_cv = cs_field_by_name_try("isobaric_heat_capacity");
  if (f_cv != nullptr) {
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

  CS_FREE_HD(pk1);

  /* Save some information */
  if (idilat == 2 && ieos != CS_EOS_NONE) {
    /* CFL conditions related to the pressure equation */
    cs_real_t *cflp = nullptr;
    cs_field_t *f_cflp = cs_field_by_name_try("algo:cfl_p");
    if (f_cflp != nullptr) {
      cflp = f_cflp->val;

      cs_arrays_set_value<cs_real_t, 1>(n_cells_ext, 0., cflp);

      ctx.wait();

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

      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        rho_k_prev[c_id] = crom_eos[c_id];
      });
    }

  }

  if (idilat == 4) {

    const cs_real_t *restrict fw = fvq->weight;

    ctx_c.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      bmasfl[f_id] *= brom[f_id];
    });

    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
      cs_lnum_t c_id_0 = i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = i_face_cells[f_id][1];

      // FIXME: should be coherent with the convective scheme of the species...
      imasfl[f_id] *=   fw[f_id]*crom[c_id_0]
                      + (1.-fw[f_id])*crom[c_id_1];
    });

  }

  /* Finalize optional post_processing algo field */
  if (f_divu != nullptr) {
    int *c_disable_flag = fvq->c_disable_flag;
    cs_lnum_t has_dc = fvq->has_disable_flag;
     ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
       cs_real_t dvol = 0;
       const int ind = has_dc * c_id;
       const int c_act = (1 - (has_dc * c_disable_flag[ind]));
       if (c_act == 1)
         dvol = 1.0/cell_f_vol[c_id];
       cpro_divu[c_id] *= dvol;
    });
  }

  ctx.wait();
  ctx_c.wait();

  /*  Free memory */
  CS_FREE_HD(taui);
  CS_FREE_HD(taub);
  CS_FREE_HD(wrk2);
  CS_FREE_HD(wrk);
  CS_FREE_HD(res);
  CS_FREE_HD(phia);
  CS_FREE_HD(_cpro_divu);
  CS_FREE_HD(gradp);
  CS_FREE_HD(rovsdt);
  CS_FREE_HD(weighf);
  CS_FREE_HD(weighb);
  CS_FREE_HD(adxk);
  CS_FREE_HD(adxkm1);
  CS_FREE_HD(dphim1);
  CS_FREE_HD(rhs0);
  CS_FREE_HD(xdtsro);
  CS_FREE_HD(tpusro);
  CS_FREE_HD(rhs);
  CS_FREE_HD(dphi);

  CS_FREE_HD(cpro_rho_tc);
  CS_FREE_HD(bpro_rho_tc);
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
 * \param[in] vel          velocity values
 * \param[in] bc_coeffs_v  boundary condition coefficients
 */
/*----------------------------------------------------------------------------*/

static void
_pressure_correction_cdo(cs_real_t              vel[][3],
                         cs_field_bc_coeffs_t  *bc_coeffs_v)
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
  if (prcdo == nullptr) bft_error(__FILE__, __LINE__, 0, _(_err_empty_prcdo));

  cs_equation_t  *eq_dp = prcdo->pressure_incr;

  /* Allocate temporary arrays */

  cs_real_3_t  *wrk = nullptr;
  BFT_MALLOC(wrk, n_cells_ext, cs_real_3_t);

  cs_field_t  *f_dp = cs_field_by_id(eq_dp->field_id);
  cs_real_t  *phi = f_dp->val;
  cs_real_t  *divu = prcdo->div_st;

  /* Associate pointers to pressure diffusion coefficients */

  cs_real_t  *restrict dt = CS_F_(dt)->val;
  cs_real_t  *c_visc = dt;

  /* Index of the field */
  const int ksinfo = cs_field_key_id("solving_info");
  cs_solving_info_t  *sinfo
    = static_cast<cs_solving_info_t*>(cs_field_get_key_struct_ptr(f_p, ksinfo));

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
                 static_cast<cs_gradient_limit_t>(eqp_u->imligr),
                 eqp_p->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 wrk,
                 bc_coeffs_v,
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

  cs_real_t *diff_flux = nullptr;
  BFT_MALLOC(diff_flux, quant->n_faces, cs_real_t);

  cs_equation_compute_diffusive_flux(eq_dp,
                                     nullptr, /* eqp --> default*/
                                     nullptr, /* diff_pty --> default */
                                     nullptr, /* dof_values --> default */
                                     nullptr, /* cell_values --> default */
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
  cs_real_t *res = nullptr;
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
  cs_pressure_correction_cdo_t *prcdo = nullptr;

  BFT_MALLOC(prcdo, 1, cs_pressure_correction_cdo_t);

  /* Equation
     -------- */

  prcdo->pressure_incr = nullptr;

  /* Fields
     ------ */

  prcdo->pressure_incr_gradient = nullptr;
  prcdo->pressure_gradient = nullptr;

  /* Other arrays
     ------------ */

  prcdo->div_st = nullptr;
  prcdo->inner_potential_flux = nullptr;
  prcdo->bdy_potential_flux = nullptr;
  prcdo->bdy_pressure_incr = nullptr;

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

  if (CS_F_(p) != nullptr) {
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

  cs_pressure_correction_cdo_t *prcdo = nullptr;

  if (cs_pressure_correction_cdo == nullptr)
    cs_pressure_correction_cdo = _pressure_correction_cdo_create();

  prcdo = cs_pressure_correction_cdo;
  assert(prcdo != nullptr);

  /* Activate the CDO module along with the FV framework */

  cs_param_cdo_mode_set(CS_PARAM_CDO_MODE_WITH_FV);

  /* Add a new equation related to the pressure correction with a variable
   * field named "pressure_increment" */

  cs_equation_t  *eq =
    cs_equation_add("pressure_increment", /* equation name */
                    "pressure_increment", /* associated variable field name */
                    CS_EQUATION_TYPE_PREDEFINED,
                    1,                        /* dimension of the unknown */
                    CS_BC_SYMMETRY); /* default boundary */

  cs_equation_param_t  *eqp = cs_equation_get_param(eq);

  /* Predefined settings */

  cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");

  /* BC settings */

  cs_equation_param_set(eqp, CS_EQKEY_BC_ENFORCEMENT, "algebraic");

  /* System to solve is SPD by construction */

  cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "fcg");
  cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
  cs_equation_param_set(eqp, CS_EQKEY_AMG_TYPE, "k_cycle");

  cs_equation_param_set(eqp, CS_EQKEY_SOLVER_MAX_ITER, "2500");
  cs_equation_param_set(eqp, CS_EQKEY_SOLVER_RTOL, "1e-5");
  cs_equation_param_set(eqp, CS_EQKEY_SOLVER_RESNORM_TYPE, "filtered");

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

  if (prcdo == nullptr) bft_error(__FILE__, __LINE__, 0, _(_err_empty_prcdo));

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
    cs_field_find_or_create("algo:pressure_gradient",
                            CS_FIELD_INTENSIVE,
                            CS_MESH_LOCATION_CELLS,
                            3,
                            false);

  cs_field_set_key_int(prcdo->pressure_gradient, post_key, 1);
  cs_field_set_key_int(prcdo->pressure_gradient, log_key, field_post_flag);

  prcdo->pressure_incr_gradient =
    cs_field_find_or_create("algo:pressure_increment_gradient",
                            CS_FIELD_INTENSIVE,
                            CS_MESH_LOCATION_CELLS,
                            3,
                            false);

  cs_field_set_key_int(prcdo->pressure_incr_gradient, post_key, 1);
  cs_field_set_key_int(prcdo->pressure_incr_gradient, log_key, field_post_flag);

 /* Activate the diffusion term
  * ---------------------------
  * Diffusivity coefficient for pressure correction equation is represented by
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

  if (prcdo == nullptr) bft_error(__FILE__, __LINE__, 0, _(_err_empty_prcdo));

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
                                       nullptr,  /* all cells */
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
                                  CS_BC_DIRICHLET,
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
  if (prcdo == nullptr)
    return;

  BFT_FREE(prcdo->div_st);
  BFT_FREE(prcdo->inner_potential_flux);
  BFT_FREE(prcdo->bdy_potential_flux);
  BFT_FREE(prcdo->bdy_pressure_incr);

  BFT_FREE(prcdo);

  cs_pressure_correction_cdo = nullptr;
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
 * \param[in]       iterns        Navier-Stokes iteration number
 * \param[in]       nfbpcd        number of faces with condensation source term
 * \param[in]       ncmast        number of cells with condensation source terms
 * \param[in]       ifbpcd        index of faces with condensation source term
 * \param[in]       ltmast        list of cells with condensation source terms
 *                                (1 to n numbering)
 * \param[in]       isostd        indicator of standard outlet and index
 *                                of the reference outlet face
 * \param[in]       vel           velocity
 * \param[in, out]  da_uu         velocity matrix
 * \param[in]       bc_coeffs_v   boundary condition structure for the variable
 * \param[in]       bc_coeffs_dp  boundary conditions structure for the
 *                                pressure increment
 * \param[in]       spcond        variable value associated to the condensation
 *                                source term (for ivar=ipr, spcond is the
 *                                flow rate
 *                                \f$ \Gamma_{s,cond}^n \f$)
 * \param[in]       svcond        variable value associated to the condensation
 *                                source term (for ivar=ipr, svcond is the flow rate
 *                                \f$ \Gamma_{v, cond}^n \f$)
 * \param[in]       frcxt         external forces making hydrostatic pressure
 * \param[in]       dfrcxt        variation of the external forces
 *                                composing the hydrostatic pressure
 * \param[in]       i_visc        visc*surface/dist aux faces internes
 * \param[in]       b_visc        visc*surface/dist aux faces de bord
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction(int                   iterns,
                       cs_lnum_t             nfbpcd,
                       cs_lnum_t             ncmast,
                       cs_lnum_t             ifbpcd[],
                       cs_lnum_t             ltmast[],
                       const int             isostd[],
                       cs_real_t             vel[][3],
                       cs_real_t             da_uu[][6],
                       cs_field_bc_coeffs_t *bc_coeffs_v,
                       cs_field_bc_coeffs_t *bc_coeffs_dp,
                       cs_real_t             spcond[],
                       cs_real_t             svcond[],
                       cs_real_t             frcxt[][3],
                       cs_real_t             dfrcxt[][3],
                       cs_real_t             i_visc[],
                       cs_real_t             b_visc[])
{
  /* Pointers to BC coefficients */

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
                            bc_coeffs_v,
                            bc_coeffs_dp,
                            spcond,
                            svcond,
                            frcxt,
                            dfrcxt,
                            i_visc,
                            b_visc);
  else
    _pressure_correction_cdo(vel,
                             bc_coeffs_v);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
