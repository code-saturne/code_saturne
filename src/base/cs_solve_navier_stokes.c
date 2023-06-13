/*============================================================================
 * Solve the Navier-Stokes equations.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_ale.h"
#include "cs_array.h"
#include "cs_assert.h"
#include "cs_atmo.h"
#include "cs_bad_cells_regularisation.h"
#include "cs_blas.h"
#include "cs_boundary_conditions.h"
#include "cs_divergence.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_lagr.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_porous_model.h"
#include "cs_pressure_correction.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_sat_coupling.h"
//#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"
#include "cs_volume_mass_injection.h"
#include "cs_wall_condensation.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_solve_navier_stokes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for Fortran functions and variables.
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Fortran function prototypes for subroutines from field.f90.
 *============================================================================*/

void
cs_f_navier_stokes_total_pressure(void);

void
cs_f_solve_navier_stokes(const int      iterns,
                         int            *icvrge,
                         const int      itrale,
                         int            impale[],
                         int            isostd[],
                         int            ale_bc_type[]);

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update flux mass for turbomachinery.
 *
 * \param[in]      m       pointer to associated mesh structure
 * \param[in]      mq      pointer to associated mesh quantities structure
 * \param[in]      crom    density at cells
 * \param[in]      brom    density at boundary faces
 * \param[in]      impale  indicator of imposed displacement
 * \param[in, out] imasfl  interior face mass flux
 * \param[in, out] bmasfl  boundary face mass flux
 */
/*----------------------------------------------------------------------------*/

static void
_turbomachinery_mass_flux(const cs_mesh_t             *m,
                          const cs_mesh_quantities_t  *mq,
                          const cs_real_t              crom[],
                          const cs_real_t              brom[],
                          cs_real_t                    imasfl[],
                          cs_real_t                    bmasfl[])
{
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_real_3_t  *restrict surfbo
    = (const cs_real_3_t  *restrict) mq->b_face_normal;
  const cs_real_3_t *restrict surfac
    = (const cs_real_3_t *restrict) mq->i_face_normal;

  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)mq->b_face_cog;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)mq->i_face_cog;

  const int *irotce = cs_turbomachinery_get_cell_rotor_num();

# pragma omp parallel for if (n_i_faces > CS_THR_MIN)
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    const cs_lnum_t c_id1 = i_face_cells[face_id][0];
    const cs_lnum_t c_id2 = i_face_cells[face_id][1];
    if ((irotce[c_id1] != 0) || (irotce[c_id2] != 0)) {
      const cs_real_t rhofac = 0.5*(crom[c_id1] + crom[c_id2]);
      cs_real_t vr1[3], vr2[3];
      const cs_rotation_t *r_num1 = cs_glob_rotation + irotce[c_id1];
      const cs_rotation_t *r_num2 = cs_glob_rotation + irotce[c_id2];
      cs_rotation_velocity(r_num1, i_face_cog[face_id], vr1);
      cs_rotation_velocity(r_num2, i_face_cog[face_id], vr2);

      imasfl[face_id] -= 0.5*rhofac*(  surfac[face_id][0]*(vr1[0] + vr2[0])
                                     + surfac[face_id][1]*(vr1[1] + vr2[1])
                                     + surfac[face_id][2]*(vr1[2] + vr2[2]));
    }
  }

# pragma omp parallel if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    cs_real_t vr[3];
    const cs_lnum_t c_id = b_face_cells[face_id];
    if (irotce[c_id] != 0) {
      const cs_real_t rhofac = brom[face_id];
      const cs_rotation_t *r_num = cs_glob_rotation + irotce[c_id];
      cs_rotation_velocity(r_num, b_face_cog[face_id], vr);

      bmasfl[face_id] -= rhofac*(  surfbo[face_id][0]*vr[0]
                                 + surfbo[face_id][1]*vr[1]
                                 + surfbo[face_id][2]*vr[2]);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief In the ALE framework, update mass flux by adding mesh velocity.
 *
 * \param[in]      m       pointer to associated mesh structure
 * \param[in]      mq      pointer to associated mesh quantities structure
 * \param[in]      dt      time step at cells
 * \param[in]      crom    density at cells
 * \param[in]      brom    density at boundary faces
 * \param[in]      impale  indicator of imposed displacement
 * \param[in, out] imasfl  interior face mass flux
 * \param[in, out] bmasfl  boundary face mass flux
 */
/*----------------------------------------------------------------------------*/

static void
_mesh_velocity_mass_flux(const cs_mesh_t             *m,
                         const cs_mesh_quantities_t  *mq,
                         const cs_real_t              dt[],
                         const cs_real_t              crom[],
                         const cs_real_t              brom[],
                         cs_real_t                    imasfl[],
                         cs_real_t                    bmasfl[])
{
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_lnum_t *i_face_vtx_idx = m->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx_lst = m->i_face_vtx_lst;
  const cs_lnum_t *b_face_vtx_idx = m->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx_lst = m->b_face_vtx_lst;

  const cs_real_3_t *vtx_coord = (const cs_real_3_t *)(m->vtx_coord);
  const cs_real_3_t *surfbo = (const cs_real_3_t *) mq->b_face_normal;
  const cs_real_3_t *surfac = (const cs_real_3_t *) mq->i_face_normal;

  const cs_real_3_t *mshvel = (const cs_real_3_t *)CS_F_(mesh_u)->val;

  const cs_real_3_t *xyzno0
    = (const cs_real_3_t *)cs_field_by_name("vtx_coord0")->val;

  const cs_real_3_t *disale
    = (const cs_real_3_t *)cs_field_by_name("mesh_displacement")->val;

  if (cs_glob_space_disc->iflxmw > 0) {

    /* One temporary array needed for internal faces,
     * in case some internal vertices are moved directly by the user */

    cs_real_t *intflx = NULL, *bouflx = NULL;
    BFT_MALLOC(intflx, n_i_faces, cs_real_t);
    BFT_MALLOC(bouflx, n_b_faces, cs_real_t);

    const cs_real_3_t *claale
      = (const cs_real_3_t *)CS_F_(mesh_u)->bc_coeffs->a;
    const cs_real_33_t *clbale
      = (const cs_real_33_t *)CS_F_(mesh_u)->bc_coeffs->b;

    const cs_equation_param_t *eqp_mesh
      = cs_field_get_equation_param_const(CS_F_(mesh_u));

    cs_mass_flux(m,
                 mq,
                 CS_F_(mesh_u)->id,
                 1,  /* itypfl */
                 1,  /* iflmb0 */
                 1,  /* init */
                 1,  /* inc */
                 eqp_mesh->imrgra,
                 eqp_mesh->nswrgr,
                 eqp_mesh->imligr,
                 eqp_mesh->iwarni,
                 eqp_mesh->epsrgr,
                 eqp_mesh->climgr,
                 crom, brom,
                 mshvel,
                 claale, clbale,
                 intflx, bouflx);

    cs_axpy(n_b_faces, -1, bouflx, bmasfl);
    cs_axpy(n_i_faces, -1, intflx, imasfl);

    BFT_FREE(intflx);
    BFT_FREE(bouflx);
  }

  /* Here we need of the opposite of the mesh velocity. */

  else { /* if (cs_glob_space_disc->iflxmw == 0) */

    /* Compute the mass flux using the nodes displacement */

#   pragma omp parallel if (n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      cs_real_t disp_fac[3] = {0, 0, 0};
      const cs_lnum_t s_id = b_face_vtx_idx[face_id];
      const cs_lnum_t e_id = b_face_vtx_idx[face_id+1];
      const cs_lnum_t icpt = e_id - s_id;
      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t inod = b_face_vtx_lst[ii];
        for (cs_lnum_t jj = 0; jj < 3; jj++)
          disp_fac[jj] +=   disale[inod][jj]
                          - (vtx_coord[inod][jj] - xyzno0[inod][jj]);
      }
      const cs_lnum_t c_id = b_face_cells[face_id];
      bmasfl[face_id] -= brom[face_id] * (  disp_fac[0]*surfbo[face_id][0]
                                          + disp_fac[1]*surfbo[face_id][1]
                                          + disp_fac[2]*surfbo[face_id][2])
                         / dt[c_id]/icpt;
    }

#   pragma omp parallel if (n_i_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
      cs_real_t disp_fac[3] = {0, 0, 0};
      const cs_lnum_t s_id = i_face_vtx_idx[face_id];
      const cs_lnum_t e_id = i_face_vtx_idx[face_id+1];
      const cs_lnum_t icpt = e_id - s_id;
      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t inod = i_face_vtx_lst[ii];
        for (cs_lnum_t jj = 0; jj < 3; jj++)
          disp_fac[jj] +=   disale[inod][jj]
                          - (vtx_coord[inod][jj] - xyzno0[inod][jj]);
      }

      /* For inner vertices, the mass flux due to the mesh displacement is
       * recomputed from the nodes displacement */
      const cs_lnum_t c_id1 = i_face_cells[face_id][0];
      const cs_lnum_t c_id2 = i_face_cells[face_id][1];
      const cs_real_t dtfac = 0.5*(dt[c_id1] + dt[c_id2]);
      const cs_real_t rhofac = 0.5*(crom[c_id1] + crom[c_id2]);
      imasfl[face_id] -= rhofac * (  disp_fac[0]*surfac[face_id][0]
                                   + disp_fac[1]*surfac[face_id][1]
                                   + disp_fac[2]*surfac[face_id][2])
                         / dtfac / icpt;
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update of the fluid velocity field.
 *
 * \param[in]      m       pointer to associated mesh structure
 * \param[in]      mq      pointer to associated mesh quantities structure
 * \param[in]      dt      time step at cells
 * \param[in]      crom    density at cells
 * \param[in]      cromk1  density at cells
 * \param[in, out] imasfl  interior face mass flux
 * \param[in, out] bmasfl  boundary face mass flux
 */
/*----------------------------------------------------------------------------*/

static void
_update_fluid_vel(const cs_mesh_t             *m,
                  const cs_mesh_quantities_t  *mq,
                  const cs_equation_param_t   *eqp_p,
                  const cs_vof_parameters_t   *vof_param,
                  const cs_real_t              dt[],
                  const cs_real_t              crom[],
                  const cs_real_t              cromk1[],
                  cs_real_t                    imasfl[],
                  cs_real_t                    bmasfl[],
                  cs_real_t                    coefa_dp[],
                  cs_real_3_t                  vel[],
                  cs_real_3_t                  dfrcxt[],
                  cs_real_3_t                  frcxt[],
                  cs_real_6_t                  dttens[],
                  const int                    isostd[])

{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)mq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)mq->b_face_cog;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;

  const cs_velocity_pressure_param_t *vp_param = cs_glob_velocity_pressure_param;

  /* irevmc = 0: Update the velocity with the pressure gradient. */

  if (vp_param->irevmc == 0) {

    /* The predicted velocity is corrected by the cell gradient of the
     * pressure increment. */
    int inc = 0;

    cs_gradient_porosity_balance(inc);

    if (vp_param->iphydr == 1 || vp_param->iifren == 1)
      inc = 1;

    /* Pressure increment gradient */

    cs_real_3_t *cpro_gradp = NULL, *gradp = NULL;
    cs_field_t *f_inc = cs_field_by_name_try("pressure_increment_gradient");
    if (f_inc != NULL)
      cpro_gradp = (cs_real_3_t *)f_inc->val;
    else {
      BFT_MALLOC(gradp, n_cells_ext, cs_real_3_t);
      cpro_gradp = gradp;
    }

    /* Scalar diffusivity */

    cs_real_t *cpro_wgrec_s = NULL;
    cs_real_6_t *cpro_wgrec_v = NULL;

    if (vof_param->vof_model != 0) {
      const int kwgrec = cs_field_key_id_try("gradient_weighting_id");
      const int iflwgr = cs_field_get_key_int(CS_F_(p), kwgrec);
      cs_field_t *f_g = cs_field_by_id(iflwgr);
      if (f_g->dim == 1) {
        cpro_wgrec_s = f_g->val;
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          cpro_wgrec_s[c_id] = dt[c_id] / crom[c_id];
        cs_mesh_sync_var_scal(cpro_wgrec_s);
      }
      else if (f_g->dim == 6) {
        cpro_wgrec_v = (cs_real_6_t *)f_g->val;
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          for (cs_lnum_t ii = 0; ii < 6; ii++)
            cpro_wgrec_v[c_id][ii] = dttens[c_id][ii] / crom[c_id];
        }
        cs_mesh_sync_var_sym_tens(cpro_wgrec_v);
      }
    }

    if (cs_glob_velocity_pressure_model->iprcdo == 0) {
      const cs_field_t *f_ddp = cs_field_by_name("pressure_increment");
      cs_field_gradient_potential(f_ddp,
                                  false,
                                  inc,
                                  vp_param->iphydr,
                                  dfrcxt,
                                  cpro_gradp);
    }

    /*  Update the velocity field */

    const cs_real_t thetap = eqp_p->thetav;

    /* Specific handling of hydrostatic pressure */

    if (vp_param->iphydr == 1) {

      /* Scalar diffusion for the pressure */
      if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION) {
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          const cs_real_t dtsrom = thetap*dt[c_id] / crom[c_id];
          const cs_real_t rhok1drhok = cromk1[c_id] / crom[c_id];
          for (cs_lnum_t isou = 0; isou < 3; isou++)
            vel[c_id][isou] =   vel[c_id][isou] * rhok1drhok
                              + dtsrom*(  dfrcxt[c_id][isou]
                                        - cpro_gradp[c_id][isou]);
        }
      }

      /* Tensorial diffusion for the pressure */
      else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) {

#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          const cs_real_t unsrom = thetap / crom[c_id];
          const cs_real_t rhok1drhok = cromk1[c_id] / crom[c_id];

          vel[c_id][0]
            =   vel[c_id][0] * rhok1drhok
              + unsrom*(  dttens[c_id][0]*(dfrcxt[c_id][0]-cpro_gradp[c_id][0])
                        + dttens[c_id][3]*(dfrcxt[c_id][1]-cpro_gradp[c_id][1])
                        + dttens[c_id][5]*(dfrcxt[c_id][2]-cpro_gradp[c_id][2]));
          vel[c_id][1]
            =   vel[c_id][1] * rhok1drhok
              + unsrom*(  dttens[c_id][3]*(dfrcxt[c_id][0]-cpro_gradp[c_id][0])
                        + dttens[c_id][1]*(dfrcxt[c_id][1]-cpro_gradp[c_id][1])
                        + dttens[c_id][4]*(dfrcxt[c_id][2]-cpro_gradp[c_id][2]));

          vel[c_id][2]
            =   vel[c_id][2] * rhok1drhok
              + unsrom*(  dttens[c_id][5]*(dfrcxt[c_id][0]-cpro_gradp[c_id][0])
                        + dttens[c_id][4]*(dfrcxt[c_id][1]-cpro_gradp[c_id][1])
                        + dttens[c_id][2]*(dfrcxt[c_id][2]-cpro_gradp[c_id][2]));
        }

      }

      /* Update of the Dirichlet boundary conditions on the
       * pressure for the outlet */

      const int *iautom = NULL;
      if (   cs_glob_atmo_option->open_bcs_treatment > 0
          && cs_glob_atmo_option->meteo_profile > 0) {
        iautom = cs_atmo_get_auto_flag();
      }

      cs_real_t *coefa_p = CS_F_(p)->bc_coeffs->a;

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        /*  automatic inlet/outlet face for atmospheric flow */
        int iautof = 0;
        if (iautom != NULL)
          iautof = iautom[face_id];

        if (isostd[face_id] == 1 || iautof > 0)
          coefa_p[face_id] += coefa_dp[face_id];
      }

    }

    /* Standard handling of hydrostatic pressure */

    else {  /* if (vp_param->iphydr == 0) */

      /* Scalar diffusion for the pressure */
      if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION) {
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          const cs_real_t dtsrom = thetap*dt[c_id] / crom[c_id];
          const cs_real_t rhok1drhok = cromk1[c_id] / crom[c_id];
          for (cs_lnum_t isou = 0; isou < 3; isou++) {
            vel[c_id][isou] =   vel[c_id][isou] * rhok1drhok
                              - dtsrom * cpro_gradp[c_id][isou];

          }
        }
      }

      /* Tensorial diffusion for the pressure */
      else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) {
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          const cs_real_t unsrom = thetap / crom[c_id];
          const cs_real_t rhok1drhok = cromk1[c_id] / crom[c_id];

          vel[c_id][0] =   vel[c_id][0] * rhok1drhok
                         - unsrom*(  dttens[c_id][0]*(cpro_gradp[c_id][0])
                                   + dttens[c_id][3]*(cpro_gradp[c_id][1])
                                   + dttens[c_id][5]*(cpro_gradp[c_id][2]));
          vel[c_id][1] =   vel[c_id][1] * rhok1drhok
                         - unsrom*(  dttens[c_id][3]*(cpro_gradp[c_id][0])
                                   + dttens[c_id][1]*(cpro_gradp[c_id][1])
                                   + dttens[c_id][4]*(cpro_gradp[c_id][2]));
          vel[c_id][2] =   vel[c_id][2] * rhok1drhok
                         - unsrom*(  dttens[c_id][5]*(cpro_gradp[c_id][0])
                                   + dttens[c_id][4]*(cpro_gradp[c_id][1])
                                   + dttens[c_id][2]*(cpro_gradp[c_id][2]));
        }
      }

    } /* vp_param->iphydr */

    if (gradp != NULL)
      BFT_FREE(gradp);
  }

  /* RT0 update from the mass fluxes */
  else { /* vp_param->irevmc != 0) */

    cs_array_real_fill_zero(3*n_cells_ext, (cs_real_t *)vel);

    /* vel = 1 / (rho Vol) SUM mass_flux (X_f - X_i) */
    if (vof_param->vof_model == 0) {

      const cs_real_t *cell_f_vol = mq->cell_f_vol;

      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        const cs_lnum_t c_id1 = i_face_cells[face_id][0];
        const cs_lnum_t c_id2 = i_face_cells[face_id][1];

        cs_real_t vol_fl_drhovol1 = 0,  vol_fl_drhovol2 = 0;

        /* If it is not a solid cell */
        if (cs_mesh_quantities_cell_is_active(mq, c_id1) == 1)
          vol_fl_drhovol1 = imasfl[face_id] / (crom[c_id1]*cell_f_vol[c_id1]);

        /* If it is not a solid cell */
        if (cs_mesh_quantities_cell_is_active(mq, c_id2) == 1)
          vol_fl_drhovol2 = imasfl[face_id] / (crom[c_id2]*cell_f_vol[c_id2]);

        vel[c_id1][0] += vol_fl_drhovol1 * (  i_face_cog[face_id][0]
                                            - cell_cen[c_id1][0]);
        vel[c_id1][1] += vol_fl_drhovol1 * (  i_face_cog[face_id][1]
                                            - cell_cen[c_id1][1]);
        vel[c_id1][2] += vol_fl_drhovol1 * (  i_face_cog[face_id][2]
                                            - cell_cen[c_id1][2]);

        vel[c_id2][0] -= vol_fl_drhovol2 * (  i_face_cog[face_id][0]
                                            - cell_cen[c_id2][0]);
        vel[c_id2][1] -= vol_fl_drhovol2 * (  i_face_cog[face_id][1]
                                            - cell_cen[c_id2][1]);
        vel[c_id2][2] -= vol_fl_drhovol2 * (  i_face_cog[face_id][2]
                                            - cell_cen[c_id2][2]);
      }

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        const cs_lnum_t c_id1 = b_face_cells[face_id];

        cs_real_t vol_fl_drhovol1 = 0;
        /* If it is not a solid cell */
        if (cs_mesh_quantities_cell_is_active(mq, c_id1) == 1)
          vol_fl_drhovol1 = bmasfl[face_id]/(crom[c_id1]*cell_f_vol[c_id1]);

        vel[c_id1][0] += vol_fl_drhovol1 * (  b_face_cog[face_id][0]
                                            - cell_cen[c_id1][0]);
        vel[c_id1][1] += vol_fl_drhovol1 * (  b_face_cog[face_id][1]
                                            - cell_cen[c_id1][1]);
        vel[c_id1][2] += vol_fl_drhovol1 * (  b_face_cog[face_id][2]
                                            - cell_cen[c_id1][2]);

      }
    }

    else { /* if (vof_param->vof_model > 1) */

      const cs_real_t *cell_f_vol = mq->cell_f_vol;

      /* Id of the volume flux */
      const int kimasf = cs_field_key_id("inner_mass_flux_id");
      const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
      const int ivolfl_id
        = cs_field_get_key_int(cs_field_by_name("void_fraction"), kimasf);
      const int bvolfl_id
        = cs_field_get_key_int(cs_field_by_name("void_fraction"), kbmasf);
      const cs_real_t *ivolfl = cs_field_by_id(ivolfl_id)->val;
      const cs_real_t *bvolfl = cs_field_by_id(bvolfl_id)->val;

      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        const cs_lnum_t c_id1 = i_face_cells[face_id][0];
        const cs_lnum_t c_id2 = i_face_cells[face_id][1];

        cs_real_t vol_fl_drhovol1 = 0, vol_fl_drhovol2 = 0;

        /* If it is not a solid cell */
        if (cs_mesh_quantities_cell_is_active(mq, c_id1) == 1)
          vol_fl_drhovol1 = ivolfl[face_id] / cell_f_vol[c_id1];

        /* If it is not a solid cell */
        if (cs_mesh_quantities_cell_is_active(mq, c_id2) == 1)
          vol_fl_drhovol2 = ivolfl[face_id] / cell_f_vol[c_id2];

        vel[c_id1][0] += vol_fl_drhovol1 * (  i_face_cog[face_id][0]
                                            - cell_cen[c_id1][0]);
        vel[c_id1][1] += vol_fl_drhovol1 * (  i_face_cog[face_id][1]
                                            - cell_cen[c_id1][1]);
        vel[c_id1][2] += vol_fl_drhovol1 * (  i_face_cog[face_id][2]
                                            - cell_cen[c_id1][2]);

        vel[c_id2][0] -= vol_fl_drhovol2 * (  i_face_cog[face_id][0]
                                            - cell_cen[c_id2][0]);
        vel[c_id2][1] -= vol_fl_drhovol2 * (  i_face_cog[face_id][1]
                                            - cell_cen[c_id2][1]);
        vel[c_id2][2] -= vol_fl_drhovol2 * (  i_face_cog[face_id][2]
                                            - cell_cen[c_id2][2]);
      }

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        const cs_lnum_t c_id1 = b_face_cells[face_id];

        cs_real_t vol_fl_drhovol1 = 0;
        /* If it is not a solid cell */
        if (cs_mesh_quantities_cell_is_active(mq, c_id1) == 1)
          vol_fl_drhovol1 = bvolfl[face_id] / cell_f_vol[c_id1];

        vel[c_id1][0] += vol_fl_drhovol1 * (  b_face_cog[face_id][0]
                                            - cell_cen[c_id1][0]);
        vel[c_id1][1] += vol_fl_drhovol1 * (  b_face_cog[face_id][1]
                                            - cell_cen[c_id1][1]);
        vel[c_id1][2] += vol_fl_drhovol1 * (  b_face_cog[face_id][2]
                                            - cell_cen[c_id1][2]);
      }

    }
  } /* vp_param->irevmc */

  cs_mesh_sync_var_vect((cs_real_t *)vel);

  if (vp_param->iphydr == 1) {
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      int i_active = cs_mesh_quantities_cell_is_active(mq, c_id);
      for (cs_lnum_t ii = 0; ii < 3; ii++)
        frcxt[c_id][ii] =   frcxt[c_id][ii]*i_active
                          + dfrcxt[c_id][ii];
    }
    cs_mesh_sync_var_vect((cs_real_t *)frcxt);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print norms of density, velocity and pressure in listing.
 *
 * \param[in]  m         pointer to associated mesh structure
 * \param[in]  mq        pointer to associated mesh quantities structure
 * \param[in]  iterns    sub-iteration count
 * \param[in]  icvrge    convergence indicator
 * \param[in]  crom      density at cells
 * \param[in]  brom      density at boundary faces
 * \param[in]  imasfl    interior face mass flux
 * \param[in]  bmasfl    boundary face mass flux
 * \param[in]  cvar_pr   pressure
 * \param[in]  cvar_vel  velocity
 */
/*----------------------------------------------------------------------------*/

static void
_log_norm(const cs_mesh_t                *m,
          const cs_mesh_quantities_t     *mq,
          int                             iterns,
          int                             icvrge,
          const cs_real_t                 crom[],
          const cs_real_t                 brom[],
          const cs_real_t                 imasfl[],
          const cs_real_t                 bmasfl[],
          const cs_real_t                 cvar_pr[],
          const cs_real_3_t               cvar_vel[])
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;
  const cs_real_t *i_face_surf = mq->i_face_surf;
  const cs_real_t *i_f_face_surf = mq->i_f_face_surf;
  const cs_real_t *b_face_surf = mq->b_face_surf;
  const cs_real_t *b_f_face_surf = mq->b_f_face_surf;

  bft_printf(" AFTER CONTINUITY PRESSURE\n"
             " -------------------------\n");
  cs_real_t rnorm = -1.0, rnormt = -1.0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    rnorm = fmax(rnorm, fabs(cvar_pr[c_id]));
  cs_parall_max(1, CS_REAL_TYPE, &rnorm);

  bft_printf("Max. pressure, %e12.4, (max. absolute value)\n", rnorm);

  rnorm = -1.0;
  cs_lnum_t imax = 1, imaxt = -1;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t vitnor = cs_math_3_norm(cvar_vel[c_id]);
   if (vitnor >= rnormt) {
     imaxt  = c_id;
     rnormt = vitnor;
   }
  }
  if (rnormt > rnorm) {
    imax = imaxt;
    rnorm = rnormt;
  }

  cs_real_t xyzmax[3] = {cell_cen[imax][0],
                         cell_cen[imax][1],
                         cell_cen[imax][2]};

  cs_parall_max_loc_vals(3, &rnorm, xyzmax);

  bft_printf("Max. velocity, %e12.4, in, %e11.3, %e11.3, %e11.3\n",
             rnorm, xyzmax[0], xyzmax[1], xyzmax[2]);

  cs_lnum_t imin = 1, imint = 1;
  rnorm = cs_math_3_norm(cvar_vel[0]);
  rnormt = rnorm;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t vitnor = cs_math_3_norm(cvar_vel[c_id]);
     if (vitnor <= rnormt) {
       imint  = c_id;
       rnormt = vitnor;
     }
  }
  if (rnormt < rnorm) {
    imin = imint;
    rnorm = rnormt;
  }

  cs_real_t xyzmin[3] = {cell_cen[imin][0],
                         cell_cen[imin][1],
                         cell_cen[imin][2]};

  cs_parall_min_loc_vals(3, &rnorm, xyzmin);

  bft_printf("Min. velocity,%e12.4, in, %e11.3, %e11.3, %e11.3\n",
             rnorm, xyzmin[0], xyzmin[1], xyzmin[2]);

  const cs_real_t *ivolfl = NULL, *bvolfl = NULL;

  const int iporos = cs_glob_porous_model;
  cs_real_t *porosi = NULL;

  /* With porosity */
  if (iporos > 0) {
    porosi = CS_F_(poro)->val;
    cs_mesh_sync_var_scal(porosi);
  }

  if (cs_glob_vof_parameters->vof_model > 0) {
    const int kimasf = cs_field_key_id("inner_mass_flux_id");
    const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
    const int ivolfl_id
      = cs_field_get_key_int(cs_field_by_name("void_fraction"), kimasf);
    const int bvolfl_id
      = cs_field_get_key_int(cs_field_by_name("void_fraction"), kbmasf);

    ivolfl = cs_field_by_id(ivolfl_id)->val;
    bvolfl = cs_field_by_id(bvolfl_id)->val;
  }

  cs_real_t rnormi = cs_math_big_r;
  cs_real_t rnorma = -cs_math_big_r;

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    const cs_lnum_t c_id1 = i_face_cells[face_id][0];
    const cs_lnum_t c_id2 = i_face_cells[face_id][1];
    cs_real_t rhom;
    if (iporos == 1 || iporos == 2)
      rhom = (porosi[c_id1]*crom[c_id1] + porosi[c_id2]*crom[c_id2])*0.5;
    else
      rhom = (crom[c_id1] + crom[c_id2])*0.5;
    /* Deal with null fluid section */
    rnorm = 0.;
    if (i_f_face_surf[face_id] / i_face_surf[face_id] > cs_math_epzero) {
      rnorm = fabs(imasfl[face_id]) / (i_f_face_surf[face_id]*rhom);
      if (cs_glob_vof_parameters->vof_model > 0)
        rnorm = fabs(ivolfl[face_id]) / i_f_face_surf[face_id];
    }
    rnorma = fmax(rnorma, rnorm);
    rnormi = fmin(rnormi, rnorm);
  }
  cs_parall_min(1, CS_REAL_TYPE, &rnormi);
  cs_parall_max(1, CS_REAL_TYPE, &rnorma);

  bft_printf(" Max. velocity at interior faces %e12.4; min. %e12.4\n",
             rnorma, rnormi);

  rnormi = cs_math_big_r;
  rnorma = -cs_math_big_r;

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    if (bvolfl != NULL) {
      /*  Deal with null fluid section */
      rnorm = 0;
      if (b_f_face_surf[face_id] / b_face_surf[face_id] > cs_math_epzero)
        rnorm = bvolfl[face_id] / (b_f_face_surf[face_id]);
    }
    else {
      const cs_lnum_t c_id = b_face_cells[face_id];
      if ((iporos == 1) || (iporos == 2))
        rnorm = bmasfl[face_id]
               / (b_face_surf[face_id]*brom[face_id]*porosi[c_id]);
      else {
      /* Deal with null fluid section */
        rnorm = 0;
        if (mq->b_f_face_surf[face_id]/mq->b_face_surf[face_id] > cs_math_epzero)
          rnorm = bmasfl[face_id]/(mq->b_f_face_surf[face_id]*brom[face_id]);
      }
    }
    rnorma = fmax(rnorma, rnorm);
    rnormi = fmin(rnormi, rnorm);
  }
  cs_parall_min(1, CS_REAL_TYPE, &rnormi);
  cs_parall_max(1, CS_REAL_TYPE, &rnorma);

  bft_printf(" Max. velocity at boundary faces %e12.4; min. %e12.4\n",
             rnorma, rnormi);

  rnorm = cs_sum(n_b_faces, bmasfl);
  cs_parall_sum(1, CS_REAL_TYPE, &rnorm);

  bft_printf(" Mass balance  at boundary: %e14.6\n", rnorm);
  bft_printf(" ------------------------------------------------------\n");

  const cs_velocity_pressure_param_t *vp_param = cs_glob_velocity_pressure_param;

  if (vp_param->nterup > 1) {
    if (icvrge == 0) {
      bft_printf(" Fixed point for velocity-pressure coupling at iteration: "
                 "%d\n", iterns);
      bft_printf("   norm = %e12.4, norm 0 = %e12.4, toler = %e12.4\n",
                 vp_param->xnrmu, vp_param->xnrmu0, vp_param->epsup);
      bft_printf(" ------------------------------------------------------\n");
      if (iterns == vp_param->nterup) {
        bft_printf(" Non convergence of fixed point for velocity-pressure "
                   "coupling"
                   " ------------------------------------------------------\n");
      }
    }
    else {
      bft_printf(" Fixed point convergence at iteration %d", iterns);
      bft_printf("   norm = %e12.4, norm 0 = %e12.4, toler = %e12.4\n",
                 vp_param->xnrmu, vp_param->xnrmu0, vp_param->epsup);
      bft_printf(" ------------------------------------------------------\n");
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print norms of density, velocity and pressure in listing.
 *
 * \param[in]  m         pointer to associated mesh structure
 * \param[in]  mq        pointer to associated mesh quantities structure
 * \param[in]  iterns    sub-iteration count
 * \param[in]  icvrge    convergence indicator
 * \param[in]  crom      density at cells
 * \param[in]  brom      density at boundary faces
 * \param[in]  imasfl    interior face mass flux
 * \param[in]  bmasfl    boundary face mass flux
 * \param[in]  cvar_pr   pressure
 * \param[in]  cvar_vel  velocity
 */
/*----------------------------------------------------------------------------*/

static void
_resize_non_interleaved_cell_arrays(const cs_mesh_t    *m,
                                    cs_lnum_t           n_sub,
                                    cs_real_t         **array)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  cs_real_t *buffer = NULL;

  BFT_MALLOC(buffer, n_sub*n_cells, cs_real_t);
  for (cs_lnum_t i = 0; i < n_sub; i++) {
    cs_array_real_copy(n_cells, *array + i*n_cells_ext, buffer + i*n_cells);
  }

  BFT_REALLOC(*array, n_sub*n_cells_ext, cs_real_t);

  for (cs_lnum_t i = 0; i < n_sub; i++) {
    cs_real_t *src = buffer + i*n_cells;
    cs_real_t *dst = *array + i*n_cells_ext;
    cs_array_real_copy(n_cells, src, dst);
    cs_mesh_sync_var_scal(dst);
  }

  BFT_FREE(buffer);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

void
cs_f_navier_stokes_total_pressure(void)
{
  cs_solve_navier_stokes_update_total_pressure(cs_glob_mesh,
                                               cs_glob_mesh_quantities,
                                               cs_glob_fluid_properties);
}

void
cs_f_solve_navier_stokes(const int   iterns,
                         int        *icvrge,
                         const int   itrale,
                         int         impale[],
                         int         isostd[],
                         int         ale_bc_type[])
{
  cs_solve_navier_stokes(iterns,
                         icvrge,
                         itrale,
                         impale,
                         isostd,
                         ale_bc_type);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update total pressure (defined as a post-processed property).
 *
 * For the compressible module, the solved pressure is already
 * the total pressure.
 *
 * Note: for Eddy Viscosity Models, the TKE may be included in the
 * solved pressure.
 *
 * \param[in]     m   pointer to mesh structure
 * \param[in]     mq  pointer to mesh quantities structure
 * \param[in]     fp  pointer to fluid properties structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_navier_stokes_update_total_pressure(const cs_mesh_t              *m,
                                             const cs_mesh_quantities_t   *mq,
                                             const cs_fluid_properties_t  *fp)
{
  /* TODO: use a function pointer here to adapt to different cases */

  cs_field_t *f = cs_field_by_name_try("total_pressure");

  if ((CS_F_(p) == NULL) || (f == NULL))
    return;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;
  const cs_real_t *gxyz = cs_get_glob_physical_constants()->gravity;
  const cs_real_t *xyzp0 = fp->xyzp0;
  const cs_real_t p0 = fp->p0, pred0 = fp->pred0, ro0 = fp->ro0;

  cs_real_t *cpro_prtot = f->val;
  const cs_real_t *cvar_pr = CS_F_(p)->val;

  const cs_real_3_t *cpro_momst = NULL;

  if (cs_glob_atmo_option->open_bcs_treatment != 0)
    cpro_momst
      = (const cs_real_3_t *)cs_field_by_name("momentum_source_terms")->val;

  /* Update cell values */

# pragma omp parallel if (n_cells > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n_cells, sizeof(cs_real_t), &s_id, &e_id);

    if (cpro_momst == NULL) {
      for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++) {
        cpro_prtot[c_id] =   cvar_pr[c_id]
                           + ro0 * cs_math_3_distance_dot_product(xyzp0,
                                                                  cell_cen[c_id],
                                                                  gxyz)
                           + p0 - pred0;
      }
    }
    else {
      for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++)
        cpro_prtot[c_id] =   cvar_pr[c_id]
                           + ro0 * cs_math_3_distance_dot_product(xyzp0,
                                                                  cell_cen[c_id],
                                                                  gxyz)
                           + p0 - pred0
                           - cs_math_3_distance_dot_product(xyzp0,
                                                            cell_cen[c_id],
                                                            cpro_momst[c_id]);
    }

    /* For Eddy Viscosity Models, "2/3 rho k"
       is included in the solved pressure */

    if (  (   cs_glob_turb_model->itytur == 2
           || cs_glob_turb_model->itytur == 5
           || cs_glob_turb_model->iturb == CS_TURB_K_OMEGA)
        && cs_glob_turb_rans_model->igrhok != 1) {

      const cs_real_t *cvar_k = CS_F_(k)->val;
      const cs_real_t *cpro_rho = CS_F_(rho)->val;

      for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++)
        cpro_prtot[c_id] -= 2.0/3 * cpro_rho[c_id]*cvar_k[c_id];
    }

  } /* cell values update */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve Navier-Stokes equations for incompressible or slightly
 *        compressible flows for one time step. Both convection-diffusion
 *        and continuity steps are performed.
 *
 * \param[in]     iterns        index of the iteration on Navier-Stokes
 * \param[in]     icvrge        convergence indicator
 * \param[in]     itrale        number of the current ALE iteration
 * \param[in]     impale        indicator of imposed displacement
 * \param[in]     isostd        indicator of standard outlet
 *                              + index of the reference face
 * \param[in]     ale_bc_type   Type of boundary for ALE
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_navier_stokes(const int   iterns,
                       int        *icvrge,
                       const int   itrale,
                       const int  *impale,
                       const int   isostd[],
                       int         ale_bc_type[])
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_lnum_t n_cells = m->n_cells;
  cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  cs_lnum_t n_i_faces = m->n_i_faces;
  cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_time_step_t *ts = cs_glob_time_step;
  const cs_wall_condensation_t *w_condensation = cs_glob_wall_condensation;
  const cs_vof_parameters_t *vof_param = cs_glob_vof_parameters;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  const cs_velocity_pressure_model_t *vp_model = cs_glob_velocity_pressure_model;
  cs_velocity_pressure_param_t *vp_param = cs_get_glob_velocity_pressure_param();

  const cs_equation_param_t *eqp_p
    = cs_field_get_equation_param_const(CS_F_(p));

  const cs_equation_param_t *eqp_u
    = cs_field_get_equation_param_const(CS_F_(vel));

  int nbrcpl = cs_sat_coupling_n_couplings();

  /* Initialization
   * -------------- */

  cs_real_t *dt = CS_F_(dt)->val;
  cs_real_t *cvar_pr = CS_F_(p)->val;
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_3_t *vela = (cs_real_3_t *)CS_F_(vel)->val_pre;

  /* Map some specific field arrays */
  cs_field_t *f_dttens = cs_field_by_name_try("dttens");
  cs_real_6_t *dttens = NULL;
  if (f_dttens != NULL)
    dttens = (cs_real_6_t *)f_dttens->val;

  /* Pointer to velocity at sub iteration k for velocity-pressure
     inner iterations */
  cs_real_3_t *uvwk = NULL, *velk = NULL;

  if (vp_param->nterup > 1) {

    const cs_real_t *cell_f_vol = mq->cell_f_vol;

    BFT_MALLOC(uvwk, n_cells_ext, cs_real_3_t);
    cs_array_real_copy(n_cells*3, (const cs_real_t *)vel, (cs_real_t *)uvwk);

    /* Compute the L2 velocity norm
       (it is zero at the first time step, so we recompute it) */
    if (iterns == 1 || fabs(vp_param->xnrmu0) <= 0) {
      cs_real_t xnrtmp = 0.0;
#     pragma omp parallel for reduction(+:xnrtmp) if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        xnrtmp += cs_math_3_dot_product(vel[c_id],
                                        vel[c_id])*cell_f_vol[c_id];
      }
      cs_parall_sum(1, CS_REAL_TYPE, &xnrtmp);
      vp_param->xnrmu0 = xnrtmp;

      /* When coupling between multiple instances of code_saturne,
         we compute the total velocity norm.
         This is required so that one instance does not stop earlier than
         the others (the numerical options should still be checked) */
      cs_real_t xnrdis[1], xnr_mu[1] = {vp_param->xnrmu0};
      for (int numcpl = 1; numcpl < nbrcpl+1; numcpl++) {
        cs_sat_coupling_array_exchange(numcpl,
                                       1, /* nbrdis */
                                       1, /* nbrloc */
                                       xnr_mu,
                                       xnrdis);
         xnr_mu[0] += xnrdis[1];
      }
      vp_param->xnrmu0 = sqrt(xnr_mu[0]);
    }

    /* Handle parallelism or periodicity of uvwk and pressure */
    cs_mesh_sync_var_scal(cvar_pr);
    cs_mesh_sync_var_vect((cs_real_t *)uvwk);
    velk = uvwk;

  }
  else
    velk = vela;

  /* Physical quantities */
  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;

  /* Pointers to properties */
  cs_real_t *crom_eos = CS_F_(rho)->val;
  const cs_real_t *brom_eos = CS_F_(rho_b)->val;
  const cs_real_t *croma = NULL, *broma = NULL;

  const cs_real_t *brom = NULL;
  cs_real_t *crom, *cpro_rho_mass = NULL, *bpro_rho_mass = NULL;

  const cs_real_t *cromk1 = NULL;
  cs_real_t *cpro_rho_k1 = NULL;
  cs_real_t *cpro_rho_tc = NULL, *bpro_rho_tc = NULL;

  if (   fluid_props->irovar == 1
      && (   vp_model->idilat > 1
          || vof_param->vof_model > 0
          || cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 3)) {

    /* If iterns = 1: this is density at time n */
    cpro_rho_mass = cs_field_by_name("density_mass")->val;
    bpro_rho_mass = cs_field_by_name("boundary_density_mass")->val;

    /* Time interpolated density */
    if (eqp_u->thetav < 1.0 && vp_param->itpcol == 0) {
      croma = CS_F_(rho)->val_pre;
      broma = CS_F_(rho_b)->val_pre;
      BFT_MALLOC(bpro_rho_tc, n_b_faces, cs_real_t);
      BFT_MALLOC(cpro_rho_tc, n_cells_ext, cs_real_t);

      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
        cpro_rho_tc[c_id] =   eqp_u->thetav * cpro_rho_mass[c_id]
                            + (1.0 - eqp_u->thetav) * croma[c_id];

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
         bpro_rho_tc[face_id] =   eqp_u->thetav * bpro_rho_mass[face_id]
                                + (1.0 - eqp_u->thetav) * broma[face_id];

      crom = cpro_rho_tc;
      cromk1 = cpro_rho_tc;  /* rho at time n+1/2,k-1 */
      brom = bpro_rho_tc;
    }
    else {
      BFT_MALLOC(cpro_rho_k1, n_cells_ext, cs_real_t);
      cs_array_real_copy(n_cells_ext, cpro_rho_mass, cpro_rho_k1);
      crom = cpro_rho_mass;
      cromk1 = cpro_rho_k1;  /* rho at time n+1/2,k-1 */
      brom = bpro_rho_mass;
    }
  }

  /* Weakly variable density algo. (idilat <=1) or constant density */
  else {
    crom = crom_eos;
    cromk1 = crom_eos;   /* rho at time n+1/2,k-1 */
    brom = brom_eos;
  }

  /* Prediction of the mass flux in case of Low Mach compressible algorithm
     ---------------------------------------------------------------------- */

  if (   (vp_model->idilat == 2 || vp_model->idilat == 3)
      && ts->nt_cur > 1
      && vp_param->ipredfl != 0)
    cs_mass_flux_prediction(dt);

  /* Hydrostatic pressure prediction in case of Low Mach compressible algorithm
     ---------------------------------------------------------------------------*/

  cs_real_3_t *grdphd = NULL;
  if (vp_param->iphydr == 2) {
    BFT_MALLOC(grdphd, n_cells_ext, cs_real_3_t);
    cs_hydrostatic_pressure_prediction(grdphd, iterns);
  }

  /* Pressure resolution and computation of mass flux for compressible flow
     ---------------------------------------------------------------------- */

  /* Note, for the compressible algorithm written in pressure increment,
   * this step is merged with the pressure correction step of the incompressible
   * algorithm */

  if (   cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1
      && cs_glob_physical_model_flag[CS_COMPRESSIBLE] != 3) {
    if (eqp_p->verbosity >= 1)
      bft_printf("** SOLVING MASS BALANCE EQUATION\n");

    cs_compressible_convective_mass_flux(iterns, dt, vela);
  }

  /* VoF: compute liquid-vapor mass transfer term (cavitating flows)
     --------------------------------------------------------------- */

  if (vof_param->vof_model & CS_VOF_MERKLE_MASS_TRANSFER) {
    const cs_real_t *cpro_prtot = cs_field_by_name("total_pressure")->val;
    const cs_real_t *cvara_voidf = cs_field_by_name("void_fraction")->val_pre;
    cs_cavitation_compute_source_term(cpro_prtot, cvara_voidf);
  }

  /* Velocity prediction step
     ------------------------ */

  bool irijnu_1 = false;
  if (   cs_glob_turb_model->itytur == 3
      && cs_glob_turb_rans_model->irijnu == 1)
    irijnu_1 = true;

  if (eqp_u->verbosity > 0)
    bft_printf("** SOLVING VELOCITY\n");

  cs_real_t *viscf = NULL, *viscb = NULL;
  cs_real_t *secvib = NULL, *secvif = NULL;
  cs_real_t *viscfi = NULL, *viscbi = NULL;
  cs_real_t *wvisbi = NULL, *wvisfi = NULL;
  cs_real_3_t *frcxt = NULL;

  static cs_real_3_t *trava = NULL;  /* TODO: pass this as argument to calling
                                        function when that is moved to C,
                                        so as to avoid requiring a static
                                        variable. */

  if (vp_param->nterup > 1 && trava == NULL)
    BFT_MALLOC(trava, n_cells_ext, cs_real_3_t);

  if (vp_model->ivisse == 1) {
    BFT_MALLOC(secvif, n_i_faces, cs_real_t);
    BFT_MALLOC(secvib, n_b_faces, cs_real_t);
  }

  if (eqp_u->idften & CS_ISOTROPIC_DIFFUSION) {
    BFT_MALLOC(viscf, n_i_faces, cs_real_t);
    BFT_MALLOC(viscb, n_b_faces, cs_real_t);
    if (irijnu_1) {
      BFT_MALLOC(wvisfi, n_i_faces, cs_real_t);
      BFT_MALLOC(wvisbi, n_b_faces, cs_real_t);
      viscfi = wvisfi;
      viscbi = wvisbi;
    }
    else {
      viscfi = viscf;
      viscbi = viscb;
    }
  }
  else if (eqp_u->idften & CS_ANISOTROPIC_LEFT_DIFFUSION) {
    BFT_MALLOC(viscb, n_b_faces, cs_real_t);
    BFT_MALLOC(viscf, 9*n_i_faces, cs_real_t);
    if (irijnu_1) {
      BFT_MALLOC(wvisbi, n_b_faces, cs_real_t);
      BFT_MALLOC(wvisfi, 9*n_i_faces, cs_real_t);
      viscfi = wvisfi;
      viscbi = wvisbi;
    }
    else {
      viscfi = viscf;
      viscbi = viscb;
    }
  }

  cs_real_3_t *trav = NULL, *dfrcxt = NULL;
  cs_real_6_t *da_uu = NULL;

  BFT_MALLOC(trav, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(da_uu, n_cells_ext, cs_real_6_t);
  BFT_MALLOC(dfrcxt, n_cells_ext, cs_real_3_t);

  if (vp_param->iphydr == 1)
    frcxt = (cs_real_3_t *)cs_field_by_name("volume_forces")->val;

  /* Pointers to BC coefficients */
  cs_real_3_t *coefau = (cs_real_3_t *)CS_F_(vel)->bc_coeffs->a;
  cs_real_3_t *cofafu = (cs_real_3_t *)CS_F_(vel)->bc_coeffs->af;
  cs_real_33_t *coefbu = (cs_real_33_t *)CS_F_(vel)->bc_coeffs->b;
  cs_real_33_t *cofbfu  = (cs_real_33_t *)CS_F_(vel)->bc_coeffs->bf;

  if (vp_param->staggered == 0)
    cs_velocity_prediction(1,
                           iterns,
                           dt,
                           vel,
                           vela,
                           velk,
                           da_uu,
                           coefau,
                           coefbu,
                           cofafu,
                           cofbfu,
                           frcxt,
                           grdphd,
                           trava,
                           dfrcxt,
                           dttens,
                           trav,
                           viscf,
                           viscb,
                           viscfi,
                           viscbi,
                           secvif,
                           secvib);
  else {
    /* Account for external forces partially balanced by the pressure gradient
       (only for the first call; the second one is for error estimators) */
    if (vp_param->iphydr == 1) {
      const cs_real_t ro0 = fluid_props->ro0;
      const cs_real_t *gxyz = cs_get_glob_physical_constants()->gravity;
#     pragma omp parallel for if (m->n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {
        const int is_active = cs_mesh_quantities_cell_is_active(mq, c_id);
        const cs_real_t drom =  (crom[c_id] - ro0) * is_active;
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          dfrcxt[c_id][ii] = drom * gxyz[ii] - frcxt[c_id][ii]*is_active;
      }
      cs_mesh_sync_var_vect((cs_real_t *)dfrcxt);
    }
  }

  /* Bad cells regularisation */
  cs_bad_cells_regularisation_vector(vel, 1);

  /* Exit if no pressure-continuity:
   * update mass fluxes and return */

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const int iflmas = cs_field_get_key_int(CS_F_(vel), kimasf);
  const int iflmab = cs_field_get_key_int(CS_F_(vel), kbmasf);

  cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  if (vp_param->iprco < 1) {
    int iflmb0 = 1;
    if (cs_glob_ale > CS_ALE_NONE)
      iflmb0 = 0;

    cs_mass_flux(m,
                 mq,
                 CS_F_(vel)->id,
                 1,  /* itypfl */
                 iflmb0,
                 1,  /* init */
                 1,  /* inc */
                 eqp_u->imrgra,
                 eqp_u->nswrgr,
                 eqp_u->imligr,
                 eqp_u->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 vel,
                 coefau, coefbu,
                 imasfl, bmasfl);

    /* In the ALE framework, we add the mesh velocity */

    if (cs_glob_ale > CS_ALE_NONE)
      _mesh_velocity_mass_flux(m, mq,
                               dt,
                               crom, brom,
                               imasfl, bmasfl);

    /* Ajout de la vitesse du solide dans le flux convectif,
     * si le maillage est mobile (solide rigide)
     * En turbomachine, on connait exactement la vitesse de maillage a ajouter */

    if (cs_turbomachinery_get_model() > CS_TURBOMACHINERY_NONE)
      _turbomachinery_mass_flux(m,
                                mq,
                                crom, brom,
                                imasfl, bmasfl);

    BFT_FREE(trav);
    BFT_FREE(da_uu);
    BFT_FREE(dfrcxt);

    BFT_FREE(viscb);
    BFT_FREE(viscf);

    BFT_FREE(secvib);
    BFT_FREE(secvif);

    BFT_FREE(grdphd);

    BFT_FREE(cpro_rho_tc);
    BFT_FREE(bpro_rho_tc);

    BFT_FREE(wvisfi);
    BFT_FREE(wvisbi);

    BFT_FREE(uvwk);

    return;
  }

  /* Update mesh for unsteady turbomachinery computations */

  cs_real_t rs_ell[2] = {0, 0};

  if (   iterns == 1
      && cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT) {

    cs_turbomachinery_update_mesh(rs_ell);

    const cs_real_t t1 = cs_timer_wtime();

    m = cs_glob_mesh;
    mq = cs_glob_mesh_quantities;
    ts = cs_glob_time_step;

    n_cells = m->n_cells;
    n_cells_ext = m->n_cells_with_ghosts;
    n_i_faces = m->n_i_faces;
    n_b_faces = m->n_b_faces;

    b_face_cells = (const cs_lnum_t *restrict)m->b_face_cells;

    if (cs_turbomachinery_get_n_couplings() < 1) {

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        /* Cancel the mass flux for symmetry BC */
        if (cs_glob_bc_type[face_id] == CS_SYMMETRY)
          mq->b_sym_flag[face_id] = 0;
        else
          mq->b_sym_flag[face_id] = 1;
      }

      /* Resize temporary internal faces arrays */

      BFT_FREE(viscf);
      if (eqp_u->idften & CS_ISOTROPIC_DIFFUSION)
        BFT_MALLOC(viscf, n_i_faces, cs_real_t);
      else if (eqp_u->idften & CS_ANISOTROPIC_LEFT_DIFFUSION)
        BFT_MALLOC(viscf, 9*n_i_faces, cs_real_t);

      if (wvisfi != NULL) {
        BFT_FREE(viscfi);
        if (eqp_u->idften == 1) {
          if (irijnu_1) {
            BFT_MALLOC(wvisfi, n_i_faces, cs_real_t);
            viscfi = wvisfi;
          }
          else
            viscfi = viscf;
        }
        else if (eqp_u->idften == 6) {
          if (irijnu_1) {
            BFT_MALLOC(wvisfi, 9*n_i_faces, cs_real_t);
            viscfi = wvisfi;
          }
          else
            viscfi = viscf;
        }
      }

      if (secvif != NULL) {
        BFT_FREE(secvif);
        BFT_MALLOC(secvif, n_i_faces, cs_real_t);
      }

      /* Resize and reinitialize main internal faces properties array */
      cs_turbomachinery_reinit_i_face_fields();

      /* Update local pointers on "internal faces" fields */
      imasfl = cs_field_by_id(iflmas)->val;

      if (cs_glob_mesh->halo != NULL) {

        cs_turbomachinery_resize_cell_fields();

        /* Update field mappings
           ("owner" fields handled by cs_turbomachinery_update);
           Remark: most of what is done in this call is redundant with the
           original initialization, and this call could probably be removed. */

        cs_field_map_and_init_bcs();
        dt = cs_field_by_name("dt")->val;

        /* Resize auxiliary arrays (pointe module) */
        cs_fortran_resize_aux_arrays();

        /* Update turbomachinery module pointers */
        cs_turbomachinery_update();

        /* Resize other arrays related to the velocity-pressure resolution */
        BFT_REALLOC(da_uu, n_cells_ext, cs_real_6_t);
        cs_mesh_sync_var_sym_tens(da_uu);

        BFT_REALLOC(trav, n_cells_ext, cs_real_3_t);
        cs_mesh_sync_var_vect((cs_real_t *)trav);

        BFT_REALLOC(dfrcxt, n_cells_ext, cs_real_3_t);
        cs_mesh_sync_var_vect((cs_real_t *)dfrcxt);

        /* Resize other arrays, depending on user options */

        if (   cs_glob_lagr_time_scheme->iilagr != CS_LAGR_OFF
            && cs_glob_lagr_dim->ntersl > 0) {
          _resize_non_interleaved_cell_arrays
            (m,
             cs_glob_lagr_dim->ntersl,
             &(cs_glob_lagr_source_terms->st_val));
        }

        if (vp_param->iphydr == 1)
          frcxt = (cs_real_3_t *)cs_field_by_name("volume_forces")->val;
        else if (vp_param->iphydr == 2) {
          BFT_REALLOC(grdphd, n_cells_ext, cs_real_3_t);
          cs_mesh_sync_var_vect((cs_real_t *)grdphd);
        }

        /* Update local pointers on "cells" fields */

        crom = CS_F_(rho)->val;
        crom_eos = CS_F_(rho)->val;

        if (   fluid_props->irovar == 1
            && (   vp_model->idilat > 1
                || vof_param->vof_model > 0
                || cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 3)) {

          cpro_rho_mass = cs_field_by_name("density_mass")->val;

          /* Time interpolated density */
          if (eqp_u->thetav < 1.0 && vp_param->itpcol == 0) {
            croma = CS_F_(rho)->val_pre;
            BFT_REALLOC(cpro_rho_tc, n_cells_ext, cs_real_t);

#           pragma omp parallel for if (n_cells_ext > CS_THR_MIN)
            for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
              cpro_rho_tc[c_id] =   eqp_u->thetav * cpro_rho_mass[c_id]
                                  + (1.0 - eqp_u->thetav) * croma[c_id];

            crom = cpro_rho_tc;
            cromk1 = cpro_rho_tc;
          }
          else {
            crom = cpro_rho_mass;
            /* rho at time n+1,k-1 */
            BFT_REALLOC(cpro_rho_k1, n_cells_ext, cs_real_t);
            cs_array_real_copy(n_cells_ext, cpro_rho_mass, cpro_rho_k1);
            cromk1 = cpro_rho_k1;
          }

        }
        else {
          crom = crom_eos;
          cromk1 = crom_eos; /* rho at time n+1,k-1 */
        }

        viscl = CS_F_(mu)->val;
        visct = CS_F_(mu_t)->val;

        vel = (cs_real_3_t *)CS_F_(vel)->val;
        vela = (cs_real_3_t *)CS_F_(vel)->val_pre;

        cvar_pr = CS_F_(p)->val;

        if (f_dttens != NULL)
          dttens = (cs_real_6_t *)f_dttens->val;

        if (vp_param->nterup > 1) {
          BFT_REALLOC(velk, n_cells_ext, cs_real_3_t);
          cs_mesh_sync_var_vect((cs_real_t *)velk);
          BFT_REALLOC(trava, n_cells_ext, cs_real_3_t);
          cs_mesh_sync_var_vect((cs_real_t *)trava);
        }
        else {
          velk = vela;
        }

      } /* halo != NULL */

    } /* cs_turbomachinery_get_n_couplings() < 1 */

    /* Update the Dirichlet wall boundary conditions for velocity (based on the
     * solid body rotation on the new mesh).
     * Note that the velocity BC update is made only if the user has
     * not specified any specific Dirichlet condition for velocity. */

    cs_real_t *coftur = NULL,  *hfltur = NULL;
    cs_turbomachinery_get_wall_bc_coeffs(&coftur, &hfltur);
    const int *irotce = cs_turbomachinery_get_cell_rotor_num();

    const cs_real_3_t *restrict b_face_u_normal
      = (const cs_real_3_t *restrict )mq->b_face_u_normal;
    const cs_real_3_t *restrict b_face_cog
      = (const cs_real_3_t *restrict)mq->b_face_cog;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      const cs_lnum_t c_id = b_face_cells[face_id];

      if (coftur[face_id] >= cs_math_infinite_r*0.5)
        continue;

      /* Physical Propreties */
      const cs_real_t visclc = viscl[c_id];
      const cs_real_t visctc = visct[c_id];

      /* Geometrical quantities */
      const cs_real_t distbf = mq->b_dist[face_id];

      /* Unit normal */
      const cs_real_t *ufn = b_face_u_normal[face_id];

      cs_real_t hint;
      if (cs_glob_turb_model->itytur == 3)
        hint = visclc / distbf;
      else
        hint = (visclc+visctc) / distbf;

      cs_real_t vr[3];
      cs_rotation_velocity(cs_glob_rotation + irotce[c_id],
                           b_face_cog[face_id],
                           vr);

      /* Gradient boundary conditions (Dirichlet) */
      const cs_real_t vrn = cs_math_3_dot_product(vr, ufn);

      for (cs_lnum_t ii = 0; ii < 3; ii++)
        coefau[face_id][ii] =   (1. - coftur[face_id]) * (vr[ii] - vrn*ufn[ii])
                              + vrn*ufn[ii];

      /* Flux boundary conditions (Dirichlet) */
      for (cs_lnum_t ii = 0; ii < 3; ii++)
        cofafu[face_id][ii] = -hfltur[face_id] * (vr[ii] - vrn*ufn[ii])
                              -hint*vrn*ufn[ii];
    }

    const cs_real_t t2 = cs_timer_wtime();

    rs_ell[1] = t2 - t1;

  } /* (iterns == 1) &&
       (cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT) */

  /* Pressure correction step
     ------------------------ */

  if (eqp_u->iwarni > 0)
    bft_printf(" ** SOLVING CONTINUITY PRESSURE\n");

  cs_real_t *coefa_dp = cs_field_by_name("pressure_increment")->bc_coeffs->a;
  cs_real_t *coefb_dp = cs_field_by_name("pressure_increment")->bc_coeffs->b;

  /* Pointers to BC coefficients */
  coefau = (cs_real_3_t *)CS_F_(vel)->bc_coeffs->a;
  coefbu = (cs_real_33_t *)CS_F_(vel)->bc_coeffs->b;

  /* Pressure correction step */
  if (   cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0
      || cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 3) {
    cs_pressure_correction(iterns,
                           w_condensation->nfbpcd,
                           w_condensation->ncmast,
                           w_condensation->ifbpcd,
                           w_condensation->ltmast,
                           isostd,
                           vel,
                           da_uu,
                           coefau,
                           coefbu,
                           coefa_dp,
                           coefb_dp,
                           w_condensation->spcond,
                           w_condensation->svcond,
                           frcxt,
                           dfrcxt,
                           viscf,
                           viscb);
  }

  /* Bad cells regularisation */
  cs_bad_cells_regularisation_scalar(cvar_pr);

  /* Update local pointers on "cells" fields */
  crom = CS_F_(rho)->val;
  crom_eos = CS_F_(rho)->val;

  /* Update density which may be computed in the pressure step */

  if (   fluid_props->irovar == 1
      && (   vp_model->idilat > 1
          || vof_param->vof_model > 0
          || cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 3)) {

    cpro_rho_mass = cs_field_by_name("density_mass")->val;

    /* Time interpolated density */
    if (eqp_u->thetav < 1.0 && vp_param->itpcol == 0) {

      croma = CS_F_(rho)->val_pre;

      if (cpro_rho_tc != NULL) {
        BFT_FREE(cpro_rho_tc);
        BFT_MALLOC(cpro_rho_tc, n_cells_ext, cs_real_t);
      }
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cpro_rho_tc[c_id] =   eqp_u->thetav * cpro_rho_mass[c_id]
                            + (1.0 - eqp_u->thetav) * croma[c_id];
      }

      cs_mesh_sync_var_scal(cpro_rho_tc);

      crom = cpro_rho_tc;
      cromk1 = cpro_rho_tc; /* rho at time n+1/2,k-1 */
    }

    else
      crom = cpro_rho_mass;

  }

  /* Mesh velocity solving (ALE) */

  if (cs_glob_ale > CS_ALE_NONE) {
    if (itrale > cs_glob_ale_n_ini_f)
      cs_ale_solve_mesh_velocity(iterns, impale, ale_bc_type);
  }

  /* Update of the fluid velocity field
     ---------------------------------- */

  if (   cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0
      || cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 3)
    _update_fluid_vel(m,
                      mq,
                      eqp_p,
                      vof_param,
                      dt,
                      crom,
                      cromk1,
                      imasfl,
                      bmasfl,
                      coefa_dp,
                      vel,
                      dfrcxt,
                      frcxt,
                      dttens,
                      isostd);

  /* Bad cells regularisation */
  cs_bad_cells_regularisation_vector(vel, 1);

  /* Mass flux initialization for VOF algorithm */
  if (vof_param->vof_model > 0) {
    cs_array_real_fill_zero(n_i_faces, imasfl);
    cs_array_real_fill_zero(n_b_faces, bmasfl);
  }

  /* In the ALE framework, we add the mesh velocity */
  if (cs_glob_ale > CS_ALE_NONE)
    _mesh_velocity_mass_flux(m, mq,
                             dt,
                             crom, brom,
                             imasfl, bmasfl);

  /* FIXME for me we should do that before predvv */
  /* Add solid's velocity in convective flux if the mesh is mobile (rigid solid).
   * For turbomachinery, the mesh velocity to add is known exactly */

  if (cs_turbomachinery_get_model() > CS_TURBOMACHINERY_NONE) {
    const cs_real_t t3 = cs_timer_wtime();
    _turbomachinery_mass_flux(m, mq,
                              crom, brom,
                              imasfl, bmasfl);
    rs_ell[1] += cs_timer_wtime() - t3;
  }

  /* VoF: void fraction solving and update the mixture density/viscosity and
   *      mass flux (cs_pressure_correction solved the convective flux of
   *      void fraction, divU)
   * ------------------------------------------------------------------------ */

  if (vof_param->vof_model > 0) {

    /* Void fraction solving */
    cs_solve_void_fraction(dt, iterns);

    /* Halo synchronization */
    cs_real_t *cvar_voidf = cs_field_by_name("void_fraction")->val;
    cs_mesh_sync_var_scal(cvar_voidf);

    /* Update mixture density/viscosity and mass flux */
    cs_vof_update_phys_prop(m);

    /* Logging */
    if (iterns == vp_param->nterup && cs_log_default_is_active())
      cs_vof_log_mass_budget(m, mq);
  }

  /* Update density (which is coherent with the mass) */

  if (   fluid_props->irovar == 1
      && (   vp_model->idilat > 1
          || vof_param->vof_model > 0
          || cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 3)) {
    cs_array_real_copy(n_cells_ext, crom_eos, cpro_rho_mass);
    cs_array_real_copy(n_cells_ext, crom_eos, crom);
    cs_array_real_copy(n_b_faces, brom_eos, bpro_rho_mass);
  }

  /* Compute error estimators for correction step and the global algorithm
     --------------------------------------------------------------------- */

  cs_field_t *iescor = cs_field_by_name_try("est_error_cor_2");
  cs_field_t *iestot = cs_field_by_name_try("est_error_tot_2");

  if (iescor != NULL || iestot != NULL) {

    const cs_real_t *cell_f_vol = mq->cell_f_vol;

    cs_real_t *esflum = NULL, *esflub = NULL;
    BFT_MALLOC(esflum, n_i_faces, cs_real_t);
    BFT_MALLOC(esflub, n_b_faces, cs_real_t);

    cs_mesh_sync_var_vect((cs_real_t *)vel);

    if (iestot != NULL)
      cs_mesh_sync_var_scal(cvar_pr);

    int iflmb0 = 1;
    if (cs_glob_ale > CS_ALE_NONE)
      iflmb0 = 0;

    /* Mass flux based on updated velocity */

    cs_mass_flux(m,
                 mq,
                 CS_F_(vel)->id,
                 1,  /* itypfl */
                 iflmb0,
                 1,  /* init */
                 1,  /* inc */
                 eqp_u->imrgra,
                 eqp_u->nswrgr,
                 eqp_u->imligr,
                 eqp_u->verbosity,
                 eqp_u->epsrgr,
                 eqp_u->climgr,
                 crom, brom,
                 vel,
                 coefau, coefbu,
                 esflum , esflub );

    /* Correction estimator: div(rom * U(n + 1)) - gamma */

    if (iescor != NULL) {
      cs_real_t *c_estim = iescor->val;
      cs_divergence(m, 1, esflum, esflub, c_estim);

      int *itpsmp = NULL;
      cs_lnum_t ncetsm = 0;
      cs_lnum_t *icetsm = NULL;
      cs_real_t *smacel, *gamma = NULL;
      cs_volume_mass_injection_get_arrays(CS_F_(p), &ncetsm, &icetsm, &itpsmp,
                                          &smacel, &gamma);

      if (ncetsm > 0) {
#       pragma omp parallel for if (ncetsm > CS_THR_MIN)
        for (cs_lnum_t c_idx = 0; c_idx < ncetsm; c_idx++) {
          cs_lnum_t c_id = icetsm[c_idx] - 1;
          c_estim[c_id] -= cell_f_vol[c_id] * smacel[c_idx];
        }
      }

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        c_estim[c_id] = fabs(c_estim[c_id]) / mq->cell_f_vol[c_id];
    }

    /* Total estimator */

    if (iestot != NULL) {
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t rovolsdt = crom[c_id] * cell_f_vol[c_id] / dt[c_id];
        for (cs_lnum_t isou = 0; isou < 3; isou++)
          trav[c_id][isou] = rovolsdt * (vela[c_id][isou] - vel[c_id][isou]);
      }

      if (vp_param->staggered == 0)
        cs_velocity_prediction(2,
                               iterns,
                               dt,
                               vel,
                               vel,
                               velk,
                               da_uu,
                               coefau,
                               coefbu,
                               (cs_real_3_t *)(CS_F_(vel)->bc_coeffs->af),
                               (cs_real_33_t *)(CS_F_(vel)->bc_coeffs->bf),
                               frcxt,
                               grdphd,
                               trava,
                               dfrcxt,
                               dttens,
                               trav,
                               viscf,
                               viscb,
                               viscfi,
                               viscbi,
                               secvif,
                               secvib);
    }

    BFT_FREE(esflum);
    BFT_FREE(esflub);

  }

  /* Velocity/pressure inner iterations
     ---------------------------------- */

  if (vp_param->nterup > 1) {

    /* Convergence test on U/P inner iterations, icvrge is 1 if converged */
    *icvrge = 1;

    const cs_real_t *cell_f_vol = mq->cell_f_vol;

    cs_real_t xnrtmp = 0;
#   pragma omp parallel for reduction(+:xnrtmp) if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t xduvw[3] = {vel[c_id][0] - velk[c_id][0],
                            vel[c_id][1] - velk[c_id][1],
                            vel[c_id][2] - velk[c_id][2]};
      xnrtmp += cs_math_3_dot_product(xduvw, xduvw) * cell_f_vol[c_id];
    }
    cs_parall_sum(1, CS_REAL_TYPE, &xnrtmp);
    vp_param->xnrmu = xnrtmp;

    cs_real_t xnr_mu[] = {vp_param->xnrmu};
    for (int numcpl = 1; numcpl < nbrcpl+1; numcpl++) {
        cs_real_t xnrdis[1];
        cs_sat_coupling_array_exchange(numcpl,
                                       1, /* nbrdis */
                                       1, /* nbrloc */
                                       xnr_mu,
                                       xnrdis);
        xnr_mu[0] += xnrdis[0];
    }
    vp_param->xnrmu = sqrt(xnr_mu[0]);

    /* Fixed-point convergence indicator */
    if (vp_param->xnrmu >= vp_param->epsup * vp_param->xnrmu0)
      *icvrge = 0;

  }

  /* Shift pressure field to set its spatial mean value to zero
   * if there is no boundary faces with a Dirichlet condition on the pressure.
   * Number of faces with Dirichlet condition for the pressure is:
   * - ndircl if idiricl = 1
   * - ndircl-1 if idircl = 0 */

  int ndircp = 0;
  if (eqp_p->ndircl == 1)
    ndircp = eqp_p->ndircl;
  else
    ndircp = eqp_p->ndircl - 1;
  if (ndircp <= 0)
    cs_field_set_volume_average(CS_F_(p), fluid_props->pred0);

  /* Compute the total pressure (defined as a post-processed property).
   * For the compressible module, the solved pressure is already the
   * total pressure.
   * Remark: for Eddy Viscosity Models,
   *         TKE might be included in the solved pressure. */

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0)
    cs_solve_navier_stokes_update_total_pressure(m, mq, fluid_props);

  if (eqp_u->verbosity > 0)
    _log_norm(m, mq,
              iterns,
              *icvrge,
              crom, brom,
              imasfl, bmasfl,
              cvar_pr,
              vel);

  if (cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT) {
    if (iterns == vp_param->nterup && cs_log_default_is_active())
      bft_printf("** INFORMATION ON UNSTEADY ROTOR/STATOR TREATMENT\n"
                 "   ----------------------------------------------\n"
                 " Time dedicated to mesh update (s): %10.4lf         \n"
                 " Global time                   (s): %10.4lf\n\n", rs_ell[0],
                 rs_ell[0] + rs_ell[1]);
  }

  BFT_FREE(trav);
  BFT_FREE(da_uu);
  BFT_FREE(dfrcxt);

  if (iterns == vp_param->nterup)
    BFT_FREE(trava);

  BFT_FREE(secvib);
  BFT_FREE(secvif);

  BFT_FREE(grdphd);

  BFT_FREE(bpro_rho_tc);
  BFT_FREE(cpro_rho_tc);

  BFT_FREE(wvisbi);
  BFT_FREE(wvisfi);

  BFT_FREE(uvwk);

  BFT_FREE(viscb);
  BFT_FREE(viscf);

  BFT_FREE(cpro_rho_k1);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
