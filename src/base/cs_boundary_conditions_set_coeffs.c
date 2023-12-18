/*============================================================================
 * Boundary condition management.
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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"

#include "cs_1d_wall_thermal.h"
#include "cs_ale.h"
#include "cs_array.h"
#include "cs_assert.h"
#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_conditions_check.h"
#include "cs_boundary_conditions_set_coeffs_symmetry.h"
#include "cs_boundary_conditions_set_coeffs_turb.h"
#include "cs_boundary_conditions_type.h"
#include "cs_cf_boundary_conditions.h"
#include "cs_coupling.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_gradient_boundary.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_gui_util.h"
#include "cs_ht_convert.h"
#include "cs_internal_coupling.h"
#include "cs_les_inflow.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mobile_structures.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_rad_transfer.h"
#include "cs_rad_transfer_bcs.h"
#include "cs_sat_coupling.h"
#include "cs_syr_coupling.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_turbomachinery.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions_set_coeffs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Macro definitions
 *----------------------------------------------------------------------------*/

#define NOZPPM 2000 /* max number of boundary conditions zone */

/*============================================================================
 * External function prototypes
 *============================================================================*/

/* Bindings to Fortran routines */

int *
cs_f_boundary_conditions_get_bc_type(void);

void
cs_f_ppprcl(int        itypfb[],
            cs_real_t  dt[]);

void
cs_f_tagmri(void);

void
cs_f_cou1di(void);

void
cs_f_mmtycl(const int  *itypfb);

void
cs_f_pptycl(bool        init,
            int        *itypfb,
            const int  *izfppp,
            cs_real_t   dt[]);

void
cs_f_cscloc(void);

void
cs_f_cscfbr(int        *itypfb,
            cs_real_t  *dt);

void
cs_f_cscfbr_init(int *itypfb);

void
cs_f_user_boundary_conditions_wrapper(const cs_lnum_t  itrifb[],
                                      int              itypfb[],
                                      const int        izfppp[],
                                      cs_real_t        dt[]);

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_boundary_conditions_set_coeffs.c
        Boundary condition management.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute boundary condition code for radiative transfer
 *----------------------------------------------------------------------------*/

static void
_boundary_condition_rt_type(const cs_mesh_t             *m,
                            const cs_mesh_quantities_t  *mq,
                            const bool                  init,
                            const int                   bc_type[])
{
  /* Unfinished function */

  CS_UNUSED(m);
  CS_UNUSED(mq);
  CS_UNUSED(init);
  CS_UNUSED(bc_type);

  const cs_lnum_t n_b_faces =  m->n_b_faces;

  cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;

  int nwsgg = rt_params->nwsgg;

  /* Initialization
   * -------------- */

  /* TODO: make a rt_bc_type as ale_bc_type */

  for (int gg_id = 0; gg_id < nwsgg; gg_id++) {

    // cs_real_t *radiance = CS_FI_(radiance, gg_id)->val;
    cs_field_t *f_rad = CS_FI_(radiance, gg_id);

    // int *icodcl_rad = NULL;
    cs_real_t *rcodcl1_rad = NULL;

    if (f_rad->bc_coeffs != NULL) {
      // icodcl_rad = f_rad->bc_coeffs->icodcl;
      rcodcl1_rad = f_rad->bc_coeffs->rcodcl1;
    }

#   pragma omp parallel for  if (n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (rcodcl1_rad[face_id] > cs_math_infinite_r*0.5)
        rcodcl1_rad[face_id] = 0;
    }

  }

}

/*----------------------------------------------------------------------------
 * Compute boundary condition code for ALE
 *----------------------------------------------------------------------------*/

static void
_boundary_condition_ale_type(const cs_mesh_t             *m,
                             const cs_mesh_quantities_t  *mq,
                             const bool                  init,
                             const cs_real_t             dt[],
                             const int                   bc_type[])
{
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const cs_real_t *surfbn = mq->b_face_surf;
  const cs_real_3_t *surfbo = (const cs_real_3_t *)mq->b_face_normal;

  /* Initialization
   * -------------- */

  int *impale = cs_glob_ale_data->impale;
  int *ale_bc_type = cs_glob_ale_data->bc_type;

  cs_real_3_t *disale
    = (cs_real_3_t *)cs_field_by_name("mesh_displacement")->val;

  cs_real_3_t *xyzno0 = (cs_real_3_t *)cs_field_by_name("vtx_coord0")->val;

  /* Set to 0 non specified rcodcl for mesh velocity arrays */

  cs_field_build_bc_codes_all();

  int *icodcl_mesh_u = NULL;
  cs_real_t *_rcodcl1_mesh_u = NULL;
  cs_real_t *rcodcl1_mesh_u = NULL;
  cs_real_t *rcodcl1_vel = CS_F_(vel)->bc_coeffs->rcodcl1;

  if (CS_F_(mesh_u)->bc_coeffs != NULL) {
    icodcl_mesh_u = CS_F_(mesh_u)->bc_coeffs->icodcl;
    rcodcl1_mesh_u = CS_F_(mesh_u)->bc_coeffs->rcodcl1;
  }

  if (cs_glob_ale == CS_ALE_CDO) {
    const int size_uma = (CS_F_(mesh_u)->dim + 1)*n_b_faces;
    BFT_MALLOC(_rcodcl1_mesh_u, size_uma, cs_real_t);
    rcodcl1_mesh_u = _rcodcl1_mesh_u;
  }

# pragma omp parallel for  if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      if (rcodcl1_mesh_u[n_b_faces*ii + face_id] > cs_math_infinite_r*0.5)
        rcodcl1_mesh_u[n_b_faces*ii + face_id] = 0;
    }
  }

  /* Check the consistency of BC types
   * --------------------------------- */

  int ierror[1] = {0};

  /* When using CDO solver, no need for checks. */
  if (cs_glob_ale == CS_ALE_CDO) {

    cs_real_3_t *b_fluid_vel = NULL;
    BFT_MALLOC(b_fluid_vel, n_b_faces, cs_real_3_t);

    cs_array_real_fill_zero(3*n_b_faces, (cs_real_t *)b_fluid_vel);

    cs_ale_update_bcs(ale_bc_type, b_fluid_vel);

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      for (cs_lnum_t ii = 0; ii < 3; ii++)
        rcodcl1_mesh_u[n_b_faces*ii + face_id] = b_fluid_vel[face_id][ii];
    }

    BFT_FREE(b_fluid_vel);

  }

  if (cs_glob_ale != CS_ALE_CDO) {

#   pragma omp parallel for  if (n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   ale_bc_type[face_id] != 0
          && ale_bc_type[face_id] != CS_BOUNDARY_ALE_FIXED
          && ale_bc_type[face_id] != CS_BOUNDARY_ALE_SLIDING
          && ale_bc_type[face_id] != CS_BOUNDARY_ALE_FREE_SURFACE
          && ale_bc_type[face_id] != CS_BOUNDARY_ALE_IMPOSED_VEL) {
        if (ale_bc_type[face_id] > 0)
          ale_bc_type[face_id] = -ale_bc_type[face_id];
        ierror[0]++;
      }
    }

    cs_parall_max(1, CS_INT_TYPE, &ierror[0]);

    if (ierror[0] != 0) {
      bft_printf("ALE METHOD\n"
                 "At least one boundary face has an unknown boundary type.\n"
                 "  The calculation will not be run."
                 "Check boundary condition definitions.");
      cs_boundary_conditions_error(ale_bc_type, NULL);
    }

    /* Conversion into BC and values
     *-------------------------------*/

    const cs_real_3_t *vtx_coord = (const cs_real_3_t *)m->vtx_coord;

    /* If all the nodes of a face have an imposed displacement, rcodcl is
       computed or overwritten, ale bc type is therfore
       CS_BOUNDARY_ALE_IMPOSED_VEL */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      int iecrw = 0, icpt  = 0;
      cs_real_t ddep_xyz[3] = {0, 0, 0};

      const cs_lnum_t s = m->b_face_vtx_idx[face_id];
      const cs_lnum_t e = m->b_face_vtx_idx[face_id+1];

      for (cs_lnum_t ii = s; ii < e; ii++) {
        const cs_lnum_t vtx_id = m->b_face_vtx_lst[ii];
        if (impale[vtx_id] == 0)
          iecrw++;
        icpt++;
        for (cs_lnum_t jj = 0; jj < 3; jj++)
          ddep_xyz[jj]
            += disale[vtx_id][jj] + xyzno0[vtx_id][jj] - vtx_coord[vtx_id][jj];
      }

      if (iecrw == 0 && ale_bc_type[face_id] != CS_BOUNDARY_ALE_SLIDING) {

        const cs_lnum_t c_id = b_face_cells[face_id];
        ale_bc_type[face_id] = CS_BOUNDARY_ALE_IMPOSED_VEL;
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          rcodcl1_mesh_u[n_b_faces*ii + face_id] = ddep_xyz[ii]/dt[c_id]/icpt;
      }
    }

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      int icpt = 0;

      if (ale_bc_type[face_id] == CS_BOUNDARY_ALE_FIXED) {
        icpt = 0;
        if (icodcl_mesh_u[face_id] == 0) {
          icpt ++;
          icodcl_mesh_u[face_id] = 1;
          rcodcl1_mesh_u[face_id] = 0;
        }
        if (icodcl_mesh_u[n_b_faces + face_id] == 0) {
          icpt ++;
          icodcl_mesh_u[n_b_faces + face_id] = 1;
          rcodcl1_mesh_u[n_b_faces + face_id] = 0;
        }
        if (icodcl_mesh_u[n_b_faces*2 + face_id] == 0) {
          icpt ++;
          icodcl_mesh_u[n_b_faces*2 + face_id] = 1;
          rcodcl1_mesh_u[n_b_faces*2 + face_id] = 0;
        }

        if (icpt == 3) {
          const cs_lnum_t s = m->b_face_vtx_idx[face_id];
          const cs_lnum_t e = m->b_face_vtx_idx[face_id+1];
          for (cs_lnum_t ii = s; ii < e; ii++) {
            const cs_lnum_t vtx_id = m->b_face_vtx_lst[ii];
            if (impale[vtx_id] == 0) {
              impale[vtx_id] = 1;
              for (cs_lnum_t jj = 0; jj < 3; jj++)
                disale[vtx_id][jj] = 0;
            }
          }
        }
      }

      /* Sliding face */
      else if (ale_bc_type[face_id] == CS_BOUNDARY_ALE_SLIDING) {
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          if (icodcl_mesh_u[n_b_faces*ii + face_id] == 0)
            icodcl_mesh_u[n_b_faces*ii + face_id] = 4;
      }

      /* Imposed mesh velocity face */
      else if (ale_bc_type[face_id] == CS_BOUNDARY_ALE_IMPOSED_VEL) {
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          if (icodcl_mesh_u[n_b_faces*ii + face_id] == 0)
            icodcl_mesh_u[n_b_faces*ii + face_id] = 1;
      }

      /* Free surface face: the mesh velocity is imposed by the mass flux */
      else if (ale_bc_type[face_id] == CS_BOUNDARY_ALE_FREE_SURFACE) {
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          if (icodcl_mesh_u[n_b_faces*ii + face_id] == 0)
            icodcl_mesh_u[n_b_faces*ii + face_id] = 1;
      }

    }

    /* Check icodcl consistency
     * ------------------------ */

    int irkerr = -1;
    int icoder[2] = {-1, -1};

    /* When called before time loop, some values are not yet available. */
    if (init == true)
      return;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      if (   icodcl_mesh_u[face_id] != 1
          && icodcl_mesh_u[face_id] != 2
          && icodcl_mesh_u[face_id] != 3
          && icodcl_mesh_u[face_id] != 4) {
        if (ale_bc_type[face_id] > 0)
          ale_bc_type[face_id] = -ale_bc_type[face_id];
        ierror[0] ++;
      }

      if (ale_bc_type[face_id] < 0) {
        irkerr = cs_glob_rank_id;
        icoder[0] = -ale_bc_type[face_id];
        icoder[1] = icodcl_mesh_u[face_id];
      }
    }

    cs_parall_max(1, CS_INT_TYPE, ierror);

    if (ierror[0] > 0) {
      cs_parall_max(1, CS_INT_TYPE, &irkerr);
      cs_parall_bcast(irkerr, 2, CS_INT_TYPE, icoder);

      bft_printf(_("ALE method\n\n"
                   "At least one boundary face has the following boundary "
                   "conditions for mesh velocity:\n\n"
                   "  ale_bc_type: %d\n"
                   "  mesh_u->bc_coeffs->icodcl[face_id]: %d\n\n"
                   "The only allowed values for icodcl are\n"
                   "  1: Dirichlet\n"
                   "  2: Convective outlet\n"
                   "  3: Neumann\n"
                   "  4: Slip\n\n"
                   "Check boundary condition definitions.\n"),
                 icoder[0], icoder[1]);

      cs_boundary_conditions_error(ale_bc_type, NULL);
    }

  } /* if (cs_glob_ale != CS_ALE_CDO) */

  /* Fluid velocity BCs for walls and symmetries
   * (due to mesh movement)
   * ------------------------------------------- */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    if (ale_bc_type[face_id] != CS_BOUNDARY_ALE_IMPOSED_VEL)
      continue;

    if (bc_type[face_id] == CS_SYMMETRY)
      for(int ii = 0; ii < 3; ii++)
        rcodcl1_vel[n_b_faces*ii + face_id]
          = rcodcl1_mesh_u[n_b_faces*ii + face_id];

    if (   bc_type[face_id] != CS_SMOOTHWALL
        && bc_type[face_id] != CS_ROUGHWALL)
      continue;

    if (   rcodcl1_vel[n_b_faces*0 + face_id] > cs_math_infinite_r*0.5
        && rcodcl1_vel[n_b_faces*1 + face_id] > cs_math_infinite_r*0.5
        && rcodcl1_vel[n_b_faces*2 + face_id] > cs_math_infinite_r*0.5) {
      for (int ii = 0; ii < 3; ii++)
        rcodcl1_vel[n_b_faces*ii + face_id]
          = rcodcl1_mesh_u[n_b_faces*ii + face_id];
    }
    else {
      for (int ii = 0; ii < 3; ii++)
        if (rcodcl1_vel[n_b_faces*ii + face_id] > cs_math_infinite_r*0.5)
          rcodcl1_vel[n_b_faces*ii + face_id] = 0;

      const cs_real_t srfbnf = surfbn[face_id];
      const cs_real_t rnxyz[3] = {surfbo[face_id][0]/srfbnf,
                                  surfbo[face_id][1]/srfbnf,
                                  surfbo[face_id][2]/srfbnf};

      const cs_real_t rcodcxyz[3] = {rcodcl1_vel[n_b_faces*0 + face_id],
                                     rcodcl1_vel[n_b_faces*1 + face_id],
                                     rcodcl1_vel[n_b_faces*2 + face_id]};
      cs_real_t rcodsn = 0;
      for (int ii = 0; ii < 3; ii++)
        rcodsn
          += (rcodcl1_mesh_u[n_b_faces*ii + face_id] - rcodcxyz[ii])*rnxyz[ii];

      for (int ii = 0; ii < 3; ii++)
        rcodcl1_vel[n_b_faces*ii + face_id] = rcodcxyz[ii] + rcodsn*rnxyz[ii];
    }

  }

  BFT_FREE(_rcodcl1_mesh_u);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Translation of the boundary conditions given by the user
 * in a form that fits to the solver.
 *
 * The values at a boundary face \f$ \fib \f$ stored in the face center
 * \f$ \centf \f$ of the variable \f$ P \f$ and its diffusive flux \f$ Q \f$
 * are written as:
 * \f[
 * P_{\face} = A_P^g + B_P^g P_{\centi}
 * \f]
 * and
 * \f[
 * Q_{\face} = A_P^f + B_P^f P_{\centi}
 * \f]
 * where \f$ P_\centi \f$ is the value of the variable \f$ P \f$ at the
 * neighboring cell.
 *
 * \warning
 * - If we consider an increment of a variable, the boundary conditions
 *   read:
 *   \f[
 *   \delta P_{\face} = B_P^g \delta P_{\centi}
 *   \f]
 *   and
 *   \f[
 *   \delta Q_{\face} = -B_P^f \delta P_{\centi}
 *   \f]
 *
 * - For a vector field such as the velocity \f$ \vect{u} \f$ the boundary
 *   conditions may read:
 *   \f[
 *   \vect{u}_{\face} = \vect{A}_u^g + \tens{B}_u^g \vect{u}_{\centi}
 *   \f]
 *   and
 *   \f[
 *   \vect{Q}_{\face} = \vect{A}_u^f + \tens{B}_u^f \vect{u}_{\centi}
 *   \f]
 *   where \f$ \tens{B}_u^g \f$ and \f$ \tens{B}_u^f \f$ are 3x3 tensor matrix
 *   which coupled velocity components next to a boundary.
 *
 * Please refer to the
 * <a href="../../theory.pdf#boundary"><b>boundary conditions</b></a> section
 * of the theory guide for more informations, as well as the
 * <a href="../../theory.pdf#condli"><b>condli</b></a> section.
 *
 * \param[in]     nvar          total number of variables
 * \param[in]     iterns        iteration number on Navier-Stokes equations
 * \param[in]     isvhb         indicator to save exchange coeffient
 *                               at the walls
 * \param[in]     itrale        ALE iteration number
 * \param[in]     italim        for ALE
 * \param[in]     itrfin        for ALE
 * \param[in]     ineefl        for ALE
 * \param[in]     itrfup        for ALE
 * \param[in,out] isostd        indicator for standard outlet
 *                              and reference face index
 * \param[out]    visvdr        dynamic viscosity after V. Driest damping in
 *                              boundary cells
 * \param[out]    hbord         exchange coefficient at boundary
 * \param[out]    theipb        value of thermal scalar at \f$ \centip \f$
 *                              of boundary cells
 * \param[in]     nftcdt        Global indicator of condensation source terms
 *                              (ie. sum on the processors of nfbpcd) cells
 *                              associated to the face with condensation
 *                              phenomenon
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_coeffs(int        nvar,
                                  int        iterns,
                                  int        isvhb,
                                  int        itrale,
                                  int        italim,
                                  int        itrfin,
                                  int        ineefl,
                                  int        itrfup,
                                  int        isostd[],
                                  cs_real_t  visvdr[],
                                  cs_real_t  hbord[],
                                  cs_real_t  theipb[],
                                  int        nftcdt)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_vertices  = mesh->n_vertices;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_real_3_t *vtx_coord = (const cs_real_3_t *)mesh->vtx_coord;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)mesh->b_face_cells;
  const cs_real_3_t *b_face_normal  = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_real_3_t *b_face_u_normal = (const cs_real_3_t *)fvq->b_face_u_normal;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)fvq->cell_cen;
  const cs_real_3_t *b_face_cog = (const cs_real_3_t *)fvq->b_face_cog;
  const cs_real_3_t *restrict diipb = (const cs_real_3_t *restrict)fvq->diipb;
  const cs_real_t   *b_face_surf    = fvq->b_face_surf;
  const cs_real_t   *b_dist         = fvq->b_dist;
  int               *isympa         = fvq->b_sym_flag;

  cs_real_t *dt = CS_F_(dt)->val;

  const int n_fields = cs_field_n_fields();
  const double cp0 = fluid_props->cp0;

  const int thermal_variable = cs_glob_thermal_model->thermal_variable;
  const int irijrb = cs_glob_turb_rans_model->irijrb;

  const int keysca  = cs_field_key_id("scalar_id");
  const int kivisl  = cs_field_key_id("diffusivity_id");
  const int kturt   = cs_field_key_id("turbulent_flux_model");
  const int kscacp  = cs_field_key_id("is_temperature");
  const int kbfid   = cs_field_key_id("boundary_value_id");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  const int kbmasf  = cs_field_key_id("boundary_mass_flux_id");

  const cs_lnum_t nt_cur = cs_glob_time_step->nt_cur;

  cs_field_t *f_th = cs_thermal_model_field();

  const int *izfppp = cs_glob_bc_pm_info->izfppp;
  int *itrifb = cs_glob_bc_pm_info->itrifb;
  int *bc_type = cs_f_boundary_conditions_get_bc_type();

  int *impale = cs_glob_ale_data->impale;
  int *ale_bc_type = cs_glob_ale_data->bc_type; /* (ialtyb en fortran) */

  /* Global number of boundary faces which are coupled with
   * a 1D wall thermal module */
  cs_gnum_t nfpt1t = cs_get_glob_1d_wall_thermal()->nfpt1t;

  /*--------------------------------------------------------------------------
   * 0) User calls
   *--------------------------------------------------------------------------*/

  cs_boundary_conditions_reset();

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >=  0)
    cs_cf_boundary_conditions_reset();

  if (cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] >=  1)
    cs_f_ppprcl(bc_type, dt);

  /* Base definitions from the GUI
     ----------------------------- */

  cs_gui_boundary_conditions_processing(bc_type);

  cs_boundary_conditions_complete(bc_type);

  /* User-defined functions
     ---------------------- */

  cs_f_user_boundary_conditions_wrapper(itrifb, bc_type, izfppp, dt);

  cs_user_boundary_conditions(cs_glob_domain, bc_type);

  /* Check consistency with GUI definitions */
  cs_gui_boundary_conditions_verify();

  /* BC'based coupling with other code_saturne instances. */

  if (cs_sat_coupling_n_couplings() > 0)
    cs_f_cscfbr(bc_type, dt);

  /* Synthetic Eddy Method for L.E.S. */
  cs_les_inflow_compute();

  /* ALE method (mesh velocity BC and vertices displacement) */
  cs_field_t *f_displ = cs_field_by_name_try("mesh_displacement");

  if (f_displ != NULL) { // (cs_glob_ale >= 1)

    cs_real_3_t *disale = (cs_real_3_t *)(f_displ->val);
    cs_real_3_t *xyzno0 = (cs_real_3_t *)cs_field_by_name("vtx_coord0")->val;

    cs_array_lnum_fill_zero(n_vertices, impale);

    /* GUI and user-defined function-based definitions */

    cs_gui_mobile_mesh_boundary_conditions(ale_bc_type, impale, disale);

    cs_user_boundary_conditions_ale(cs_glob_domain,
                                    bc_type,
                                    ale_bc_type,
                                    impale);

    /* In case the user has modified disale whthout setting impale=1, we restore
       the initial displacement. */

    for (cs_lnum_t ii = 0; ii < n_vertices; ii++) {
      if (impale[ii] == 0) {
        disale[ii][0] = vtx_coord[ii][0] - xyzno0[ii][0];
        disale[ii][1] = vtx_coord[ii][1] - xyzno0[ii][1];
        disale[ii][2] = vtx_coord[ii][2] - xyzno0[ii][2];
      }
    }

    /* In case of structures coupling, compute a predicted displacement */

    cs_mobile_structures_prediction(itrale, italim, ineefl, impale);
  }

  /* Once some BC codes are set by the user, they can be completed by couplings
     at boundaries (such as Syrthes), unless we need to do an additional pass
     so as to centralize what is relative to the Syrthes coupling.y
     We place here the  call relative to the Syrthes volume coupling so as
     to benefit from the lates computed velocity if we loop on U/P.
     The colume coupling must be called before the surface one to follow
     the communication scheme. */

  if (itrfin == 1 && itrfup == 1) {

    cs_syr_coupling_exchange_volume();

    cs_syr_coupling_recv_boundary(nvar, bc_type);

    if (nfpt1t > 0)
      cs_f_cou1di();

    /* Coupling 1D thermal model with condensation modelling
       to take into account the solid temperature evolution over time */
    if (nftcdt > 0)
      cs_f_tagmri();
  }

  /* For internal coupling, set itypfb to wall function by default
     if not set by the user */
  cs_internal_coupling_bcs(bc_type);

  /* Radiative transfer: add contribution to energy BCs. */
  if (cs_glob_rad_transfer_params->type > 0 && itrfin == 1 && itrfup == 1)
    cs_rad_transfer_bcs(bc_type);

  /* Convert temperature to enthalpy for Dirichlet conditions */

  cs_lnum_t nbt2h = 0;
  cs_real_t *vbt2h = NULL;
  cs_lnum_t *lbt2h = NULL;

  if (thermal_variable == CS_THERMAL_MODEL_ENTHALPY) {

    BFT_MALLOC(lbt2h, n_b_faces, cs_lnum_t);
    BFT_MALLOC(vbt2h, n_b_faces, cs_real_t);

    cs_field_t *f_h = CS_F_(h);
    cs_real_t *rcodcl1_h = f_h->bc_coeffs->rcodcl1;
    int       *icodcl_h  = f_h->bc_coeffs->icodcl;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      if (icodcl_h[f_id] < 0) {
        lbt2h[nbt2h]    = f_id;
        nbt2h           = nbt2h + 1;
        icodcl_h[f_id]  = - icodcl_h[f_id];
        vbt2h[f_id]     =   rcodcl1_h[f_id];
      }
      else
        vbt2h[f_id] = 0.0;
    }

    cs_ht_convert_t_to_h_faces_l(nbt2h, lbt2h, vbt2h, rcodcl1_h);
  }

  /*--------------------------------------------------------------------------
   * 1) initializations
   *--------------------------------------------------------------------------*/

  /* Allocate temporary arrays */
  cs_real_3_t *velipb = NULL;
  BFT_MALLOC(velipb, n_b_faces, cs_real_3_t);

  cs_turb_model_type_t iturb  = cs_glob_turb_model->iturb;
  int itytur = cs_glob_turb_model->itytur;

  /* coefa and coefb are required to compute the cell gradients for the wall
     turbulent boundary conditions.
     so, their initial values are kept (note that at the first time step,
     they are initialized to zero flux in inivar.f90) */

  /* velipb stores the velocity in i' of boundary cells */

  /* initialize variables to avoid compiler warnings */

  cs_real_t rinfiv[3] = {cs_math_infinite_r,
                         cs_math_infinite_r,
                         cs_math_infinite_r};

  /* Pointers to y+, t+ and t* if saved */
  cs_real_t *tplusp = NULL, *tstarp = NULL, *yplbr = NULL;

  /* Initialization of the array storing yplus
     which is computed in clptur.f90 and/or clptrg.f90 */

  cs_field_t *yplus = cs_field_by_name_try("yplus");
  if (yplus != NULL) {
    yplbr = yplus->val;
    cs_array_real_fill_zero(n_b_faces, yplbr);
  }

  cs_field_t *itplus = cs_field_by_name_try("tplus");
  if (itplus != NULL) {
    tplusp = itplus->val;
    cs_array_real_fill_zero(n_b_faces, tplusp);
  }

  cs_field_t *itstar = cs_field_by_name_try("tstar");
  if (itstar != NULL) {
    tstarp = itstar->val;
    cs_array_real_fill_zero(n_b_faces, tstarp);
  }

  /* Map field arrays */
  cs_field_t *vel = CS_F_(vel);
  int *icodcl_vel = vel->bc_coeffs->icodcl;
  const cs_real_3_t *var_vela = (const cs_real_3_t *)vel->val_pre;
  cs_equation_param_t *eqp_vel = cs_field_get_equation_param(vel);

  /* Pointers to the mass fluxes */
  const cs_real_t *b_massflux
    = cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf))->val;

  cs_real_t   *bfconv  = NULL, *bhconv = NULL;
  cs_real_3_t *forbr   = NULL;
  const cs_real_6_t *dttens = NULL;

  /* Pointers to specific fields */

  if (cs_glob_rad_transfer_params->type >= 1) {
    bfconv = cs_field_by_name("rad_convective_flux")->val;
    bhconv = cs_field_by_name("rad_exchange_coefficient")->val;
  }

  cs_field_t *f_dttens  = cs_field_by_name_try("dttens");
  if (f_dttens != NULL)
    dttens = (const cs_real_6_t *)f_dttens->val;

  cs_field_t *f_forbr = cs_field_by_name_try("boundary_forces");

  if (f_forbr != NULL && iterns == 1)
    forbr = (cs_real_3_t *)f_forbr->val;

  /*--------------------------------------------------------------------------
   * 2) Treatment of types of bcs given by itypfb
   *--------------------------------------------------------------------------*/

  {
    if (   (   cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] >=  1
            && cs_glob_physical_model_flag[CS_GAS_MIX]             == -1
            && cs_glob_physical_model_flag[CS_JOULE_EFFECT]        == -1
            && cs_glob_physical_model_flag[CS_ELECTRIC_ARCS]       == -1)
        || (   cs_glob_physical_model_flag[CS_COMPRESSIBLE]       >=  0
            && cs_glob_physical_model_flag[CS_GAS_MIX]            >=  0)) {

      cs_f_pptycl(false,
                  bc_type,
                  izfppp,
                  dt);
    }

    if (cs_glob_ale > CS_ALE_NONE)
      _boundary_condition_ale_type(mesh,
                                   fvq,
                                   false,
                                   dt,
                                   bc_type);

    if (cs_glob_rad_transfer_params->type != CS_RAD_TRANSFER_NONE)
      _boundary_condition_rt_type(mesh,
                                  fvq,
                                  false,
                                  bc_type);

    if (cs_turbomachinery_get_model() != CS_TURBOMACHINERY_NONE)
      cs_f_mmtycl(bc_type);

    cs_boundary_conditions_type(false,
                                bc_type,
                                itrifb,
                                isostd);
  }

  /*--------------------------------------------------------------------------
   * 3) Check the consistency of the bcs
   *--------------------------------------------------------------------------*/

  cs_boundary_conditions_check(bc_type,
                               ale_bc_type);

  /*--------------------------------------------------------------------------
   * 4) Variables
   *--------------------------------------------------------------------------*/

  /* Physical quantities */

  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;

  /*--------------------------------------------------------------------------
   * 5) Compute the temperature or the enthalpy in i' for boundary cells
   *    (thanks to the formula: fi + grad(fi).ii')
   *
   * For the coupling with syrthes
   *  theipb is used by cs_syr_coupling_send_boundary after condli.
   * For the coupling with the 1d wall thermal module
   *  theipb is used by cou1do after this function.
   * For the radiation module
   *  theipb is used to compute the required flux in raypar.
   *
   * This could be done outside the loop.
   *--------------------------------------------------------------------------*/

  {
    /* For the Syrthes coupling or 1d thermal module
       ---------------------------------------------
       Here we do an extraneous loop (we od something only foricpsyr = 1).
       This is to prepare the eventual handling of multiple temperatures
       (i.e. multiple simultaneous Syrthes couplings, for example for sensitivity
       analysis). Note that even in this case, one temperature only is received
       from each coupling; in multiphase cases, the enthalpies need to be
       reconstructed...
       Here, there can be only one scalar with icpsyr = 1 and then only
       of there is actually a coupling with Syrthes.

       For the 1D module, we use the thermal variable.

       For the radiative module
       ------------------------
       We compute the value at I' if there is a thermal variable.

       We search for the only scalar which matches;
       it can be t, h, or e (in the compressible case)

       Compute the boundary value of required scalars

       Check for boundary values
    */

    for (int ii = 0; ii < n_fields; ii++) {

      cs_field_t *f_scal = cs_field_by_id(ii);

      if (!(f_scal->type & CS_FIELD_VARIABLE))
        continue;
      if (cs_field_get_key_int(f_scal, keysca) <= 0)
        continue;

      cs_field_t  *f_scal_b = NULL;
      cs_real_t   *bvar_s = NULL;
      cs_real_3_t *bvar_v = NULL;
      cs_real_6_t *bvar_t = NULL;

      int b_f_id = cs_field_get_key_int(f_scal, kbfid);

      if (b_f_id > -1)
        f_scal_b = cs_field_by_id(b_f_id);
      else {
        /* if thermal variable has no boundary but temperature does, use it */
        if (f_scal == f_th && f_scal == CS_F_(h))
          f_scal_b = cs_field_by_name_try("boundary_temperature");
      }

      if (f_scal_b == NULL && f_scal != f_th)
        continue; /* nothing to do for this scalar */

      if (f_scal_b != NULL) {
        if (f_scal_b->dim == 1)
          bvar_s = f_scal_b->val;
        else if (f_scal_b->dim == 3)
          bvar_v = (cs_real_3_t *)f_scal_b->val;
        else if (f_scal_b->dim == 6)
          bvar_t = (cs_real_6_t *)f_scal_b->val;
      }

      cs_equation_param_t *eqp_scal = cs_field_get_equation_param(f_scal);

      if (f_scal->dim == 1) {

        if (cs_glob_space_disc->itbrrb == 1 && eqp_scal->ircflu == 1) {
          cs_real_t *var_iprime = theipb;
          if (f_scal_b != NULL)
            var_iprime = bvar_s;

          cs_field_gradient_boundary_iprime_scalar(f_scal,
                                                   true, /* use_previous_t */
                                                   n_b_faces,
                                                   NULL,
                                                   var_iprime);
        }
        else { /* itbrrb, ircflu */
          const cs_real_t *cvara_s = f_scal->val_pre;

          if (f_scal_b != NULL) {
            for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
              cs_lnum_t c_id = b_face_cells[f_id];
              bvar_s[f_id] = cvara_s[c_id];
            }
          }
          else {
            for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
              cs_lnum_t c_id = b_face_cells[f_id];
              theipb[f_id] = cvara_s[c_id];
            }
          }
        }

        /* Copy bvar_s to theipb if both theipb and bvar_s present */

        if (f_scal_b != NULL && f_th == f_scal) {
          for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
            theipb[f_id] = bvar_s[f_id];
        }

      }

      else if (f_scal_b != NULL) {

        if (f_scal->dim == 3) {

          if (cs_glob_space_disc->itbrrb == 1 && eqp_scal->ircflu == 1) {
            const cs_real_3_t *cvar_v
              = (const cs_real_3_t *)f_scal->val;

            cs_real_33_t *gradv = NULL;
            BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

            const int inc = 1;
            const bool iprev = true;

            cs_field_gradient_vector(f_scal, iprev, inc, gradv);

            for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
              cs_lnum_t c_id = b_face_cells[f_id];
              for (int isou = 0; isou < 3; isou++) {
                bvar_v[f_id][isou] =   cvar_v[c_id][isou]
                                     + cs_math_3_dot_product(gradv[c_id][isou],
                                                             diipb[f_id]);
              }
            }

            BFT_FREE(gradv);
          }
          else {
            const cs_real_3_t *cvara_v
              = (const cs_real_3_t *)f_scal->val_pre;

            for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
              cs_lnum_t c_id = b_face_cells[f_id];
              for (int isou = 0; isou < 3; isou++)
                bvar_v[f_id][isou] = cvara_v[c_id][isou];
            }
          }

        }
        else if (f_scal->dim == 6) {

          if (cs_glob_space_disc->itbrrb == 1 && eqp_scal->ircflu == 1) {
            const cs_real_6_t *cvar_t
              = (const cs_real_6_t *)f_scal->val;

            cs_real_63_t *gradt = NULL;
            BFT_MALLOC(gradt, n_cells_ext, cs_real_63_t);

            const int inc = 1;
            const bool iprev = true;

            cs_field_gradient_tensor(f_scal, iprev, inc, gradt);

            for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
              cs_lnum_t c_id = b_face_cells[f_id];
              for (int isou = 0; isou < 6; isou++) {
                bvar_t[f_id][isou] =   cvar_t[c_id][isou]
                                     + cs_math_3_dot_product(gradt[c_id][isou],
                                                             diipb[f_id]);
              }
            }

            BFT_FREE(gradt);
          }
          else {
            const cs_real_6_t *cvara_t
              = (const cs_real_6_t *)f_scal->val_pre;

            for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
              cs_lnum_t c_id = b_face_cells[f_id];
              for (int isou = 0; isou < 6; isou++)
                bvar_t[f_id][isou] = cvara_t[c_id][isou];
            }
          }

        }

        else {
          cs_assert(0);
        }
      }

    } /* end of loop on scalar fields */
  }

  /*--------------------------------------------------------------------------
   * 6) Compute the velocity and Reynolds stesses tensor in i' for boundary
   *    cells (thanks to the formula: fi + grad(fi).ii') if there are
   *    symmetry or wall faces with wall functions boundary conditions
   *--------------------------------------------------------------------------*/

  /* Indicator for symmetries or wall with wall functions */

  int iclsym = 0, ipatur = 0, ipatrg = 0;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    if (icodcl_vel[f_id] == 4)
      iclsym = 1;
    else if (icodcl_vel[f_id] == 5)
      ipatur = 1;
    else if (icodcl_vel[f_id] == 6)
      ipatrg = 1;

    if (iclsym != 0 && ipatur != 0 && ipatrg != 0)
      break;
  }

  int have_bc_flag[3] = {iclsym, ipatur, ipatrg};
  cs_parall_max(3, CS_INT_TYPE, have_bc_flag);
  iclsym = have_bc_flag[0];
  ipatur = have_bc_flag[1];
  ipatrg = have_bc_flag[2];

  /* Compute the velocity in i' for boundary cells */

  if (iclsym != 0 || ipatur != 0 || ipatrg != 0 || f_forbr != NULL) {

    if (nt_cur > 1 && eqp_vel->ircflu == 1) {
      cs_field_gradient_boundary_iprime_vector(vel,
                                               true, /* use_previous_t */
                                               n_b_faces,
                                               NULL,
                                               velipb);
    }

    /* nb: at the first time step, coefa and coefb are unknown, so the walue
           in i is stored instead of the value in i' */
    else {
      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        cs_lnum_t c_id = b_face_cells[f_id];
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          velipb[f_id][isou] = var_vela[c_id][isou];
        }
      }
    }
  }

  /* Compute rij in i' for boundary cells */

  cs_real_6_t *rijipb = NULL;
  if ((iclsym != 0 || ipatur != 0 || ipatrg != 0) && itytur == 3) {

    /* Allocate a work array to store rij values at boundary faces */
    BFT_MALLOC(rijipb, n_b_faces, cs_real_6_t);

    cs_equation_param_t *eqp_rij = cs_field_get_equation_param(CS_F_(rij));

    if (nt_cur > 1 && irijrb == 1 && eqp_rij->ircflu == 1) {
      cs_field_gradient_boundary_iprime_tensor(vel,
                                               true, /* use_previous_t */
                                               n_b_faces,
                                               NULL,
                                               rijipb);
    }

    /* nb: at the first time step, coefa and coefb are unknown, so the value
           in i is stored instead of the value in i' */
    else {
      const cs_real_6_t *cvara_ts = (const cs_real_6_t *)CS_F_(rij)->val_pre;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        cs_lnum_t c_id = b_face_cells[f_id];
        for (int isou = 0; isou < 6; isou++) {
          rijipb[f_id][isou] = cvara_ts[c_id][isou];
        }
      }
    }
  }

  /*--------------------------------------------------------------------------
   * 6) turbulence at walls:
   *    (velocity, k, epsilon, rij, temperature)
   *
   * We need velipb and rijipb (and theipb for radiation).
   *
   * Initialize visvdr to -999.d0.
   * In clptur, we damp the turbulent viscosity at wall cells if Van Driest
   * is activated. The final value is stored in visvdr.
   * Later on, in distyp, the viscosity at wall cells will be damped again;
   * we use visvdr to restore the correct value.
   *--------------------------------------------------------------------------*/

  if (itytur == 4 && cs_glob_turb_les_model->idries == 1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      visvdr[c_id] = -999.0;
  }

  if (ipatur != 0 || ipatrg != 0)
    /* Smooth and rough wall laws */
    cs_boundary_conditions_set_coeffs_turb(isvhb,
                                           velipb,
                                           rijipb,
                                           visvdr,
                                           hbord,
                                           theipb);

  /*--------------------------------------------------------------------------
   * 6) Symmetry for vectors and tensors
   *    (velocity, rij)
   *
   * We need velipb and rijipb.
   *--------------------------------------------------------------------------*/

  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    isympa[i] = 1;

  if (iclsym != 0)
    cs_boundary_conditions_set_coeffs_symmetry(velipb, rijipb);

  BFT_FREE(rijipb);

  /*--------------------------------------------------------------------------
   * 9) velocity: outlet, Dirichlet and Neumann and convective outlet
   *--------------------------------------------------------------------------*/

  { /* outlet: in case of incoming mass flux, the mass flux is set to zero. */

    cs_lnum_t isoent = 0;
    cs_lnum_t isorti = 0;

    cs_real_3_t  *coefa_vel = (cs_real_3_t  *)vel->bc_coeffs->a;
    cs_real_33_t *coefb_vel = (cs_real_33_t *)vel->bc_coeffs->b;
    cs_real_3_t  *cofaf_vel = (cs_real_3_t  *)vel->bc_coeffs->af;
    cs_real_33_t *cofbf_vel = (cs_real_33_t *)vel->bc_coeffs->bf;

    for (int f_id = 0; f_id < n_b_faces; f_id++) {

      if (icodcl_vel[f_id] == 9) {

        const cs_real_t flumbf = b_massflux[f_id];

        /* physical properties */
        const cs_lnum_t c_id = b_face_cells[f_id];
        const cs_real_t visclc  = viscl[c_id];
        const cs_real_t visctc  = visct[c_id];
        cs_real_t hint = 0.0;

        /* geometric quantities */
        const cs_real_t distbf = b_dist[f_id];

        if (itytur == 3)
          hint = visclc / distbf;
        else
          hint = (visclc + visctc) / distbf;

        isorti = isorti + 1;

        if (flumbf < - cs_math_epzero) {

          /* Dirichlet boundary condition
             ---------------------------- */

          cs_real_3_t pimpv = {0., 0., 0.};

          /* coupled solving of the velocity components */

          cs_boundary_conditions_set_dirichlet_vector(coefa_vel[f_id],
                                                      cofaf_vel[f_id],
                                                      coefb_vel[f_id],
                                                      cofbf_vel[f_id],
                                                      pimpv,
                                                      hint,
                                                      rinfiv);

            isoent = isoent + 1;
        }
        else {

          /* Neumann boundary conditions
             --------------------------- */

          cs_real_3_t qimpv = {0., 0., 0.};

          /* coupled solving of the velocity components */

          cs_boundary_conditions_set_neumann_vector(coefa_vel[f_id],
                                                    cofaf_vel[f_id],
                                                    coefb_vel[f_id],
                                                    cofbf_vel[f_id],
                                                    qimpv,
                                                    hint);

        }
      }
    }

    if (nt_cur%cs_glob_log_frequency == 0 || eqp_vel->verbosity >= 0) {
      cs_gnum_t isocpt[2] = {isoent, isorti};
      cs_parall_sum(2, CS_GNUM_TYPE, isocpt);
      if (isocpt[1] > 0 && (eqp_vel->verbosity >= 2 || isocpt[0] > 0))
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("Incoming flow detained for %llu out of %llu outlet faces\n"),
           (unsigned long long)isocpt[0],
           (unsigned long long)isocpt[1]);
    }

    /* Dirichlet and Neumann */

    cs_real_t *rcodcl1_vel = vel->bc_coeffs->rcodcl1;
    cs_real_t *rcodcl2_vel = vel->bc_coeffs->rcodcl2;
    cs_real_t *rcodcl3_vel = vel->bc_coeffs->rcodcl3;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      const cs_lnum_t c_id = b_face_cells[f_id];

      /* physical properties */
      const cs_real_t visclc = viscl[c_id];
      const cs_real_t visctc = visct[c_id];
      cs_real_t hint = 0.0;

      /* geometric quantities */
      const cs_real_t distbf = b_dist[f_id];

      if (itytur == 3)
        hint = visclc / distbf;
      else
        hint = (visclc + visctc) / distbf;

      cs_real_t qimpv[3]   = {0., 0., 0.};
      cs_real_t hextv[3]   = {0., 0., 0.};
      cs_real_t pimpv[3]   = {0., 0., 0.};
      cs_real_t cflv[3]    = {0., 0., 0.};

      /* Dirichlet boundary conditions
         ----------------------------- */

      if (icodcl_vel[f_id] == 1) {

        for (cs_lnum_t k = 0; k < 3; k++)
          pimpv[k] = rcodcl1_vel[n_b_faces*k + f_id];

        for (cs_lnum_t k = 0; k < 3; k++)
          hextv[k] = rcodcl2_vel[n_b_faces*k + f_id];

        cs_boundary_conditions_set_dirichlet_vector(coefa_vel[f_id],
                                                    cofaf_vel[f_id],
                                                    coefb_vel[f_id],
                                                    cofbf_vel[f_id],
                                                    pimpv,
                                                    hint,
                                                    hextv);

      }

      /* Neumann boundary conditions
         --------------------------- */

      else if (icodcl_vel[f_id] == 3) {

        /* coupled solving of the velocity components */

        for (cs_lnum_t k = 0; k < 3; k++)
          qimpv[k] = rcodcl3_vel[n_b_faces*k + f_id];

        cs_boundary_conditions_set_neumann_vector(coefa_vel[f_id],
                                                  cofaf_vel[f_id],
                                                  coefb_vel[f_id],
                                                  cofbf_vel[f_id],
                                                  qimpv,
                                                  hint);
      }

      /* Convective boundary conditions
         ------------------------------ */

      else if (icodcl_vel[f_id] == 2 && iterns <= 1) {

        /* Coupled solving of the velocity components */
        for (cs_lnum_t k = 0; k < 3; k++)
          pimpv[k] = rcodcl1_vel[n_b_faces*k + f_id];

        for (cs_lnum_t k = 0; k < 3; k++)
          cflv[k] = rcodcl2_vel[n_b_faces*k + f_id];

        cs_boundary_conditions_set_convective_outlet_vector(coefa_vel[f_id],
                                                            cofaf_vel[f_id],
                                                            coefb_vel[f_id],
                                                            cofbf_vel[f_id],
                                                            pimpv,
                                                            cflv,
                                                            hint);
      }

      /* Imposed value for the convection operator, imposed flux for diffusion
         --------------------------------------------------------------------- */

      else if (icodcl_vel[f_id] == 13) {

        for (cs_lnum_t k = 0; k < 3; k++)
          pimpv[k] = rcodcl1_vel[n_b_faces*k + f_id];

        for (cs_lnum_t k = 0; k < 3; k++)
          qimpv[k] = rcodcl3_vel[n_b_faces*k + f_id];

        cs_boundary_conditions_set_dirichlet_conv_neumann_diff_vector
          (coefa_vel[f_id], cofaf_vel[f_id],
           coefb_vel[f_id], cofbf_vel[f_id],
           pimpv, qimpv);

      }

      /* Convective boundary for Marangoni effects
         (generalized symmetry condition)
         ----------------------------------------- */

      else if (icodcl_vel[f_id] == 14) {

        for (cs_lnum_t k = 0; k < 3; k++)
          pimpv[k] = rcodcl1_vel[n_b_faces*k + f_id];

        for (cs_lnum_t k = 0; k < 3; k++)
          qimpv[k] = rcodcl3_vel[n_b_faces*k + f_id];

        /* Coupled solving of the velocity components */

        cs_boundary_conditions_set_generalized_sym_vector(coefa_vel[f_id],
                                                          cofaf_vel[f_id],
                                                          coefb_vel[f_id],
                                                          cofbf_vel[f_id],
                                                          pimpv,
                                                          qimpv,
                                                          hint,
                                                          b_face_u_normal[f_id]);
      }

      /* Neumann on the normal component, Dirichlet on tangential components
         ------------------------------------------------------------------- */

      else if (icodcl_vel[f_id] == 11) {

        /* Dirichlet to impose on the tangential components */
        for (cs_lnum_t k = 0; k < 3; k++)
          pimpv[k] = rcodcl1_vel[n_b_faces*k + f_id];

        /* Flux to impose on the normal component */
        for (cs_lnum_t k = 0; k < 3; k++)
          qimpv[k] = rcodcl3_vel[n_b_faces*k + f_id];

        /* coupled solving of the velocity components */

        cs_boundary_conditions_set_generalized_dirichlet_vector
          (coefa_vel[f_id], cofaf_vel[f_id],
           coefb_vel[f_id], cofbf_vel[f_id],
           pimpv, qimpv, hint, b_face_u_normal[f_id]);
      }
    }
  }

  /*--------------------------------------------------------------------------
   * 10) Pressure: Dirichlet and Neumann and convective outlet
   *--------------------------------------------------------------------------*/

  {
    cs_field_t *p = CS_F_(p);

    cs_real_t *coefa_p = p->bc_coeffs->a;
    cs_real_t *coefb_p = p->bc_coeffs->b;
    cs_real_t *cofaf_p = p->bc_coeffs->af;
    cs_real_t *cofbf_p = p->bc_coeffs->bf;

    const int *icodcl_p = (const int *)p->bc_coeffs->icodcl;
    const cs_real_t *rcodcl1_p = (const cs_real_t *)p->bc_coeffs->rcodcl1;
    const cs_real_t *rcodcl2_p = (const cs_real_t *)p->bc_coeffs->rcodcl2;
    const cs_real_t *rcodcl3_p = (const cs_real_t *)p->bc_coeffs->rcodcl3;

    cs_equation_param_t *eqp_p = cs_field_get_equation_param(p);
    cs_real_t *crom = NULL;

    if (cs_glob_vof_parameters->vof_model > 0)
      crom = CS_F_(rho)->val; // FIXME consistency with correction step

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      const cs_lnum_t c_id = b_face_cells[f_id];
      cs_real_t hint  = 0.0;

      /* geometric quantities */
      const cs_real_t distbf = b_dist[f_id];
      cs_real_t visci[3][3], dist[3];
      const cs_real_t surf = b_face_surf[f_id];
      const cs_real_t *n = b_face_normal[f_id];

      /* if a flux dt.grad p (w/m2) is set in cs_user_boundary_conditions */
      if (eqp_p->idften & CS_ISOTROPIC_DIFFUSION) {
        hint = dt[c_id]/distbf;

        if (cs_glob_vof_parameters->vof_model > 0)
          hint = hint / crom[c_id];
      }
      else if (eqp_p->idften & CS_ORTHOTROPIC_DIFFUSION) {

        hint = (dttens[c_id][0] * cs_math_pow2(n[0])
              + dttens[c_id][1] * cs_math_pow2(n[1])
              + dttens[c_id][2] * cs_math_pow2(n[2]))
              / (cs_math_pow2(surf) * distbf);

        if (cs_glob_vof_parameters->vof_model > 0)
          hint = hint / crom[c_id];
      }

      /* symmetric tensor diffusivity */
      else if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION) {

        visci[0][0] = dttens[c_id][0];
        visci[1][1] = dttens[c_id][1];
        visci[2][2] = dttens[c_id][2];
        visci[0][1] = dttens[c_id][3];
        visci[1][0] = dttens[c_id][3];
        visci[1][2] = dttens[c_id][4];
        visci[2][1] = dttens[c_id][4];
        visci[0][2] = dttens[c_id][5];
        visci[2][0] = dttens[c_id][5];

        dist[0] = b_face_cog[f_id][0] - cell_cen[c_id][0];
        dist[1] = b_face_cog[f_id][1] - cell_cen[c_id][1];
        dist[2] = b_face_cog[f_id][2] - cell_cen[c_id][2];

        // ||ki.s||^2
        const cs_real_t viscis = cs_math_pow2(  visci[0][0]*n[0]
                                              + visci[1][0]*n[1]
                                              + visci[2][0]*n[2])
                               + cs_math_pow2(  visci[0][1]*n[0]
                                              + visci[1][1]*n[1]
                                              + visci[2][1]*n[2])
                               + cs_math_pow2(  visci[0][2]*n[0]
                                              + visci[1][2]*n[1]
                                              + visci[2][2]*n[2]);

        // if.ki.s
        cs_real_t fikis
          = (  cs_math_3_dot_product(dist, visci[0]) * n[0]
             + cs_math_3_dot_product(dist, visci[1]) * n[1]
             + cs_math_3_dot_product(dist, visci[2]) * n[2]);

        /* take i" so that i"f= eps*||fi||*ki.n when j" is in cell rji
           nb: eps =1.d-1 must be consistent with vitens.f90 */
        fikis = cs_math_fmax(fikis, 1.e-1*sqrt(viscis)*distbf);

        hint = viscis / surf / fikis;
        if (cs_glob_vof_parameters->vof_model > 0)
          hint = hint / crom[c_id];
      }

      /* We must modify the Dirichlet pressure value again so as to obtain P*
         Because in cs_boundary_conditions_type we have used the total pressure
         provided by the user: ptotal= p*+ rho.g.r.
         In the compressible case, we leave rcodcl as such. */

      /* Dirichlet boundary condition
         ----------------------------- */

      if (icodcl_p[f_id] == 1) {

        const cs_real_t hext = rcodcl2_p[f_id];
        const cs_real_t pimp = rcodcl1_p[f_id];

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_p[f_id],
                                                    &cofaf_p[f_id],
                                                    &coefb_p[f_id],
                                                    &cofbf_p[f_id],
                                                    pimp,
                                                    hint,
                                                    hext);
      }

      /* Neumann boundary conditions
         ---------------------------- */

      if (icodcl_p[f_id] == 3) {

        const cs_real_t dimp = rcodcl3_p[f_id];

        cs_boundary_conditions_set_neumann_scalar(&coefa_p[f_id],
                                                  &cofaf_p[f_id],
                                                  &coefb_p[f_id],
                                                  &cofbf_p[f_id],
                                                  dimp,
                                                  hint);
      }

      /* Convective boundary conditions
         ------------------------------ */

      else if (icodcl_p[f_id] == 2) {

        const cs_real_t pimp = rcodcl1_p[f_id];
        const cs_real_t cfl  = rcodcl2_p[f_id];

        cs_boundary_conditions_set_convective_outlet_scalar(&coefa_p[f_id],
                                                            &cofaf_p[f_id],
                                                            &coefb_p[f_id],
                                                            &cofbf_p[f_id],
                                                            pimp,
                                                            cfl,
                                                            hint);

      }

      /* Boundary value proportional to boundary cell value
         -------------------------------------------------- */

      else if (icodcl_p[f_id] == 10) {

        const cs_real_t pinf  = rcodcl1_p[f_id];
        const cs_real_t ratio = rcodcl2_p[f_id];

        cs_boundary_conditions_set_affine_function_scalar(&coefa_p[f_id],
                                                          &cofaf_p[f_id],
                                                          &coefb_p[f_id],
                                                          &cofbf_p[f_id],
                                                          pinf,
                                                          ratio,
                                                          hint);
      }

      /* Imposed value for the convection operator is proportional to boundary
         cell value, imposed flux for diffusion
         --------------------------------------------------------------------- */

      else if (icodcl_p[f_id] == 12) {

        const cs_real_t pinf  = rcodcl1_p[f_id];
        const cs_real_t ratio = rcodcl2_p[f_id];
        const cs_real_t dimp  = rcodcl3_p[f_id];

        cs_boundary_conditions_set_affine_function_conv_neumann_diff_scalar
          (&coefa_p[f_id], &cofaf_p[f_id],
           &coefb_p[f_id], &cofbf_p[f_id],
           pinf, ratio, dimp);
      }

      /* Imposed value for the convection operator, imposed flux for diffusion
         --------------------------------------------------------------------- */

      else if (icodcl_p[f_id] == 13) {

        const cs_real_t pimp = rcodcl1_p[f_id];
        const cs_real_t dimp = rcodcl3_p[f_id];

        cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
          (&coefa_p[f_id], &cofaf_p[f_id],
           &coefb_p[f_id], &cofbf_p[f_id],
           pimp, dimp);
      }

      /* Neumann for the convection operator, zero flux for diffusion
         ------------------------------------------------------------ */

      else if (icodcl_p[f_id] == 15) {

        const cs_real_t dimp = rcodcl3_p[f_id];

        cs_boundary_conditions_set_neumann_conv_h_neumann_diff_scalar
          (&coefa_p[f_id], &cofaf_p[f_id],
           &coefb_p[f_id], &cofbf_p[f_id],
           dimp, hint);
      }
    }
  } /* pressure */

  /*--------------------------------------------------------------------------
   * 11) void fraction (VOF): Dirichlet and Neumann and convective outlet
   *--------------------------------------------------------------------------*/

  if (cs_glob_vof_parameters->vof_model > 0) {

    cs_field_t *volf2 = CS_F_(void_f);

    cs_real_t *coefa_vol = volf2->bc_coeffs->a;
    cs_real_t *coefb_vol = volf2->bc_coeffs->b;
    cs_real_t *cofaf_vol = volf2->bc_coeffs->af;
    cs_real_t *cofbf_vol = volf2->bc_coeffs->bf;

    int *icodcl_vol = volf2->bc_coeffs->icodcl;
    cs_real_t *rcodcl1_vol = volf2->bc_coeffs->rcodcl1;
    cs_real_t *rcodcl2_vol = volf2->bc_coeffs->rcodcl2;
    cs_real_t *rcodcl3_vol = volf2->bc_coeffs->rcodcl3;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      /* hint is unused since there is no diffusion for the void fraction */
      const cs_real_t hint = 1.0;

      /* Dirichlet boundary condition
         ---------------------------- */

      if (icodcl_vol[f_id] == 1) {

        const cs_real_t pimp = rcodcl1_vol[f_id];
        const cs_real_t hext = rcodcl2_vol[f_id];

        cs_boundary_conditions_set_dirichlet_scalar(&coefa_vol[f_id],
                                                    &cofaf_vol[f_id],
                                                    &coefb_vol[f_id],
                                                    &cofbf_vol[f_id],
                                                    pimp,
                                                    hint,
                                                    hext);
      }

      /* Neumann boundary conditions
         --------------------------- */

      if (icodcl_vol[f_id] == 3) {

        const cs_real_t dimp = rcodcl3_vol[f_id];

        cs_boundary_conditions_set_neumann_scalar(&coefa_vol[f_id],
                                                  &cofaf_vol[f_id],
                                                  &coefb_vol[f_id],
                                                  &cofbf_vol[f_id],
                                                  dimp,
                                                  hint);
      }

      /* Convective boundary conditions
         ------------------------------ */

      else if (icodcl_vol[f_id] == 2) {

        const cs_real_t pimp = rcodcl1_vol[f_id];
        const cs_real_t cfl  = rcodcl2_vol[f_id];

        cs_boundary_conditions_set_convective_outlet_scalar
          (&coefa_vol[f_id], &cofaf_vol[f_id],
           &coefb_vol[f_id], &cofbf_vol[f_id],
           pimp, cfl, hint);
      }
    }
  } /* VOF */

  /*----------------------------------------------------------------------------
    12. turbulent quantities: Dirichlet and Neumann and convective outlet
    --------------------------------------------------------------------------*/

  { /* k-epsilon and k-omega */

    if (itytur == 2 || iturb == CS_TURB_K_OMEGA) {

      cs_field_t *turb = NULL;
      cs_real_t sigma = 0.0;

      for (int ii = 0; ii < 2; ii++) {

        /* For k-omega, use sigma_k2 and sigma_w2 values as this term
           is in practice only for inlets (no issue at walls or with 0 flux). */

        if (ii == 0) {
          turb  = CS_F_(k);
          if (itytur == 2)
            sigma = cs_field_get_key_double(turb, ksigmas);
          else if (iturb == CS_TURB_K_OMEGA) {
            sigma = cs_turb_ckwsk2; /* FIXME: not consistent with the model */
          }
        }
        else {
          if (itytur == 2) {
            turb  = CS_F_(eps);
            sigma = cs_field_get_key_double(turb, ksigmas);
          }
          else {
            turb  = CS_F_(omg);
            sigma = cs_turb_ckwsw2; /* FIXME: not consistent with the model */
          }
        }

        cs_real_t *coefa_turb = turb->bc_coeffs->a;
        cs_real_t *coefb_turb = turb->bc_coeffs->b;
        cs_real_t *cofaf_turb = turb->bc_coeffs->af;
        cs_real_t *cofbf_turb = turb->bc_coeffs->bf;

        int *icodcl_turb = turb->bc_coeffs->icodcl;
        cs_real_t *rcodcl1_turb = turb->bc_coeffs->rcodcl1;
        cs_real_t *rcodcl2_turb = turb->bc_coeffs->rcodcl2;
        cs_real_t *rcodcl3_turb = turb->bc_coeffs->rcodcl3;

        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

          const cs_lnum_t c_id = b_face_cells[f_id];

          /* physical properties */
          const cs_real_t visclc = viscl[c_id];
          const cs_real_t visctc = visct[c_id];

          /* geometric quantities */
          const cs_real_t distbf = b_dist[f_id];
          const cs_real_t hint   = (visclc + visctc / sigma) / distbf;

          /* Dirichlet boundary condition
             ----------------------------- */

          if (icodcl_turb[f_id] == 1) {

            const cs_real_t pimp = rcodcl1_turb[f_id];
            const cs_real_t hext = rcodcl2_turb[f_id];

            cs_boundary_conditions_set_dirichlet_scalar(&coefa_turb[f_id],
                                                        &cofaf_turb[f_id],
                                                        &coefb_turb[f_id],
                                                        &cofbf_turb[f_id],
                                                        pimp,
                                                        hint,
                                                        hext);

          }

          /* Neumann boundary conditions
             ---------------------------- */

          if (icodcl_turb[f_id] == 3) {

            const cs_real_t dimp = rcodcl3_turb[f_id];

            cs_boundary_conditions_set_neumann_scalar(&coefa_turb[f_id],
                                                      &cofaf_turb[f_id],
                                                      &coefb_turb[f_id],
                                                      &cofbf_turb[f_id],
                                                      dimp,
                                                      hint);
          }

          /* convective boundary conditions
             ------------------------------- */

          else if (icodcl_turb[f_id] == 2) {

            const cs_real_t pimp = rcodcl1_turb[f_id];
            const cs_real_t cfl  = rcodcl2_turb[f_id];

            cs_boundary_conditions_set_convective_outlet_scalar
              (&coefa_turb[f_id], &cofaf_turb[f_id],
               &coefb_turb[f_id], &cofbf_turb[f_id],
               pimp, cfl, hint);
          }

          /* Imposed value for the convection operator,
             imposed flux for diffusion
             ---------------------------------------- */

          else if (icodcl_turb[f_id] == 13) {

            const cs_real_t pimp = rcodcl1_turb[f_id];
            const cs_real_t dimp = rcodcl3_turb[f_id];

            cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
              (&coefa_turb[f_id], &cofaf_turb[f_id],
               &coefb_turb[f_id], &cofbf_turb[f_id],
               pimp, dimp);
          }
        }
      }
    }

    /* Rij-epsilon */

    else if (itytur == 3) {

      cs_field_t *rij = CS_F_(rij);

      cs_real_6_t  *coefa_ts = (cs_real_6_t  *)rij->bc_coeffs->a;
      cs_real_66_t *coefb_ts = (cs_real_66_t *)rij->bc_coeffs->b;
      cs_real_6_t  *cofaf_ts = (cs_real_6_t  *)rij->bc_coeffs->af;
      cs_real_66_t *cofbf_ts = (cs_real_66_t *)rij->bc_coeffs->bf;
      cs_real_6_t  *cofad_ts = (cs_real_6_t  *)rij->bc_coeffs->ad;
      cs_real_66_t *cofbd_ts = (cs_real_66_t *)rij->bc_coeffs->bd;

      int *icodcl_ts = rij->bc_coeffs->icodcl;
      cs_real_t *rcodcl1_ts = rij->bc_coeffs->rcodcl1;
      cs_real_t *rcodcl2_ts = rij->bc_coeffs->rcodcl2;
      cs_real_t *rcodcl3_ts = rij->bc_coeffs->rcodcl3;

      cs_equation_param_t *eqp_ts = cs_field_get_equation_param(rij);
      cs_field_t *f_a_t_visc = NULL;
      cs_real_6_t *visten = NULL;

      if (eqp_ts->idften & CS_ANISOTROPIC_DIFFUSION) {
        f_a_t_visc = cs_field_by_name("anisotropic_turbulent_viscosity");
        visten = (cs_real_6_t *)f_a_t_visc->val;
      }

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

        cs_real_t pimpts[6] = {0., 0., 0., 0., 0., 0.};
        cs_real_t hextts[6] = {0., 0., 0., 0., 0., 0.};
        cs_real_t qimpts[6] = {0., 0., 0., 0., 0., 0.};
        cs_real_t cflts[6]  = {0., 0., 0., 0., 0., 0.};

        const cs_lnum_t c_id = b_face_cells[f_id];

        /* physical properties */
        const cs_real_t visclc = viscl[c_id];
        cs_real_t hint = 0.;

        /* geometric quantities */
        const cs_real_t distfi = b_dist[f_id];
        const cs_real_t surf = b_face_surf[f_id];
        const cs_real_t *n = b_face_normal[f_id];
        cs_real_t visci[3][3], dist[3];

        dist[0] = b_face_cog[f_id][0] - cell_cen[c_id][0];
        dist[1] = b_face_cog[f_id][1] - cell_cen[c_id][1];
        dist[2] = b_face_cog[f_id][2] - cell_cen[c_id][2];

        /* symmetric tensor diffusivity (Daly Harlow - GGDH) TODO */
        if (eqp_ts->idften & CS_ANISOTROPIC_RIGHT_DIFFUSION) {

          visci[0][0] = visclc + visten[c_id][0];
          visci[1][1] = visclc + visten[c_id][1];
          visci[2][2] = visclc + visten[c_id][2];
          visci[0][1] =          visten[c_id][3];
          visci[1][0] =          visten[c_id][3];
          visci[1][2] =          visten[c_id][4];
          visci[2][1] =          visten[c_id][4];
          visci[0][2] =          visten[c_id][5];
          visci[2][0] =          visten[c_id][5];

          /* ||ki.s||^2 */
          const cs_real_t viscis = cs_math_pow2(  visci[0][0]*n[0]
                                                + visci[1][0]*n[1]
                                                + visci[2][0]*n[2])
                                 + cs_math_pow2(  visci[0][1]*n[0]
                                                + visci[1][1]*n[1]
                                                + visci[2][1]*n[2])
                                 + cs_math_pow2(  visci[0][2]*n[0]
                                                + visci[1][2]*n[1]
                                                + visci[2][2]*n[2]);

          /* if.ki.s */
          cs_real_t fikis
            = (  cs_math_3_dot_product(dist, visci[0]) * n[0]
               + cs_math_3_dot_product(dist, visci[1]) * n[1]
               + cs_math_3_dot_product(dist, visci[2]) * n[2]);

          /* take i" so that i"f= eps*||fi||*ki.n when j" is in cell rji
             nb: eps =1.d-1 must be consistent with vitens.f90 */
          fikis = cs_math_fmax(fikis, 1.e-1*sqrt(viscis)*distfi);

          hint = viscis / surf / fikis;
        }

        /* scalar diffusivity */
        else {
          const cs_real_t visctc = visct[c_id];
          hint = (visclc + visctc * cs_turb_csrij / cs_turb_cmu) / distfi;
        }

        /* Dirichlet Boundary Condition
           ---------------------------- */

        if (icodcl_ts[f_id] == 1) {

          for (cs_lnum_t k = 0; k < 6; k++)
            pimpts[k] = rcodcl1_ts[n_b_faces*k + f_id];

          for (cs_lnum_t k = 0; k < 6; k++)
            hextts[k] = rcodcl2_ts[n_b_faces*k + f_id];

          cs_boundary_conditions_set_dirichlet_tensor(coefa_ts[f_id],
                                                      cofaf_ts[f_id],
                                                      coefb_ts[f_id],
                                                      cofbf_ts[f_id],
                                                      pimpts,
                                                      hint,
                                                      hextts);

          /* Boundary conditions for the momentum equation */
          for (cs_lnum_t isou = 0; isou < 6; isou++) {
            cofad_ts[f_id][isou]       = coefa_ts[f_id][isou];
            cofbd_ts[f_id][isou][isou] = coefb_ts[f_id][isou][isou];
          }
        }

        /* Neumann Boundary Condition
           -------------------------- */

        else if (icodcl_ts[f_id] == 3) {

          for (cs_lnum_t k = 0; k < 6; k++)
            qimpts[k] = rcodcl3_ts[n_b_faces*k + f_id];

          cs_boundary_conditions_set_neumann_tensor(coefa_ts[f_id],
                                                    cofaf_ts[f_id],
                                                    coefb_ts[f_id],
                                                    cofbf_ts[f_id],
                                                    qimpts,
                                                    hint);

          /* Boundary conditions for the momentum equation */
          for (cs_lnum_t isou = 0; isou < 6; isou++) {
            cofad_ts[f_id][isou]       = coefa_ts[f_id][isou];
            cofbd_ts[f_id][isou][isou] = coefb_ts[f_id][isou][isou];
          }
        }

        /* Convective Boundary Condition
           ----------------------------- */

        else if (icodcl_ts[f_id] == 2) {

          for (cs_lnum_t k = 0; k < 6; k++)
            pimpts[k] = rcodcl1_ts[n_b_faces*k + f_id];

          for (cs_lnum_t k = 0; k < 6; k++)
            cflts[k] = rcodcl2_ts[n_b_faces*k + f_id];

          cs_boundary_conditions_set_convective_outlet_tensor
            (coefa_ts[f_id], cofaf_ts[f_id],
             coefb_ts[f_id], cofbf_ts[f_id],
             pimpts, cflts, hint);

          /* Boundary conditions for the momentum equation */
          for (cs_lnum_t isou = 0; isou < 6; isou++) {
            cofad_ts[f_id][isou]       = coefa_ts[f_id][isou];
            cofbd_ts[f_id][isou][isou] = coefb_ts[f_id][isou][isou];
          }
        }

        /* Imposed value for the convection operator,
           imposed flux for diffusion
           ------------------------------------------ */

        else if (icodcl_ts[f_id] == 13) {

          for (cs_lnum_t k = 0; k < 6; k++)
            pimpts[k] = rcodcl1_ts[n_b_faces*k + f_id];

          for (cs_lnum_t k = 0; k < 6; k++)
            qimpts[k] = rcodcl3_ts[n_b_faces*k + f_id];

          cs_boundary_conditions_set_dirichlet_conv_neumann_diff_tensor
            (coefa_ts[f_id], cofaf_ts[f_id],
             coefb_ts[f_id], cofbf_ts[f_id],
             pimpts, qimpts);

          /* Boundary conditions for the momentum equation */
          for (cs_lnum_t isou = 0; isou < 6; isou++) {
            cofad_ts[f_id][isou]       = coefa_ts[f_id][isou];
            cofbd_ts[f_id][isou][isou] = coefb_ts[f_id][isou][isou];
          }
        }

      }

      /* epsilon */

      cs_field_t *eps = CS_F_(eps);

      cs_real_t *coefa_eps = eps->bc_coeffs->a;
      cs_real_t *coefb_eps = eps->bc_coeffs->b;
      cs_real_t *cofaf_eps = eps->bc_coeffs->af;
      cs_real_t *cofbf_eps = eps->bc_coeffs->bf;

      int *icodcl_eps = eps->bc_coeffs->icodcl;
      cs_real_t *rcodcl1_eps = eps->bc_coeffs->rcodcl1;
      cs_real_t *rcodcl2_eps = eps->bc_coeffs->rcodcl2;
      cs_real_t *rcodcl3_eps = eps->bc_coeffs->rcodcl3;

      cs_real_t sigmae = cs_field_get_key_double(eps, ksigmas);

      cs_equation_param_t *eqp_eps = cs_field_get_equation_param(eps);
      f_a_t_visc = NULL, visten = NULL;

      if (eqp_eps->idften & CS_ANISOTROPIC_DIFFUSION) {
        f_a_t_visc = cs_field_by_name("anisotropic_turbulent_viscosity");
        visten = (cs_real_6_t *)f_a_t_visc->val;
      }

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

        const cs_lnum_t c_id = b_face_cells[f_id];
        const cs_real_t surf = b_face_surf[f_id];

        /* --- Physical Properties */
        const cs_real_t visclc = viscl[c_id];
        const cs_real_t visctc = visct[c_id];
        cs_real_t hint = 0.;

        /* Geometric quantities */
        const cs_real_t distfi = b_dist[f_id];
        cs_real_t visci[3][3], dist[3];
        const cs_real_t *n = b_face_normal[f_id];

        dist[0] = b_face_cog[f_id][0] - cell_cen[c_id][0];
        dist[1] = b_face_cog[f_id][1] - cell_cen[c_id][1];
        dist[2] = b_face_cog[f_id][2] - cell_cen[c_id][2];

        /* Symmetric tensor diffusivity (Daly Harlow - GGDH) */
        if (eqp_eps->idften & CS_ANISOTROPIC_DIFFUSION) {

          visci[0][0] = visclc + visten[c_id][0]/sigmae;
          visci[1][1] = visclc + visten[c_id][1]/sigmae;
          visci[2][2] = visclc + visten[c_id][2]/sigmae;
          visci[0][1] =          visten[c_id][3]/sigmae;
          visci[1][0] =          visten[c_id][3]/sigmae;
          visci[1][2] =          visten[c_id][4]/sigmae;
          visci[2][1] =          visten[c_id][4]/sigmae;
          visci[0][2] =          visten[c_id][5]/sigmae;
          visci[2][0] =          visten[c_id][5]/sigmae;

          /* ||Ki.S||^2 */
          const cs_real_t viscis = cs_math_pow2(  visci[0][0]*n[0]
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
          fikis = cs_math_fmax(fikis, 1.e-1*sqrt(viscis)*distfi);

          hint = viscis / surf / fikis;
        }

        /* Scalar diffusivity */

        else
          hint = (visclc + visctc / sigmae) / distfi;

        /* Dirichlet Boundary Condition
           ---------------------------- */

        if (icodcl_eps[f_id] == 1) {

          const cs_real_t pimp = rcodcl1_eps[f_id];
          const cs_real_t hext = rcodcl2_eps[f_id];

          cs_boundary_conditions_set_dirichlet_scalar(&coefa_eps[f_id],
                                                      &cofaf_eps[f_id],
                                                      &coefb_eps[f_id],
                                                      &cofbf_eps[f_id],
                                                      pimp,
                                                      hint,
                                                      hext);
        }

        /* Neumann Boundary Condition
           -------------------------- */

        else if (icodcl_eps[f_id] == 3) {

          const cs_real_t dimp = rcodcl3_eps[f_id];

          cs_boundary_conditions_set_neumann_scalar(&coefa_eps[f_id],
                                                    &cofaf_eps[f_id],
                                                    &coefb_eps[f_id],
                                                    &cofbf_eps[f_id],
                                                    dimp,
                                                    hint);
        }

        /* Convective Boundary Condition
           ----------------------------- */

        else if (icodcl_eps[f_id] == 2) {

          const cs_real_t pimp = rcodcl1_eps[f_id];
          const cs_real_t cfl =  rcodcl2_eps[f_id];

          cs_boundary_conditions_set_convective_outlet_scalar
            (&coefa_eps[f_id], &cofaf_eps[f_id],
             &coefb_eps[f_id], &cofbf_eps[f_id],
             pimp, cfl, hint);
        }

        /* Imposed value for the convection operator,
           imposed flux for diffusion
           ----------------------------------------- */

        else if (icodcl_eps[f_id] == 13) {

          const cs_real_t pimp = rcodcl1_eps[f_id];
          const cs_real_t dimp = rcodcl3_eps[f_id];

          cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
            (&coefa_eps[f_id], &cofaf_eps[f_id],
             &coefb_eps[f_id], &cofbf_eps[f_id],
             pimp, dimp);
        }
      }

      /* Alpha for the EBRSM */

      if (iturb == CS_TURB_RIJ_EPSILON_EBRSM) {

        cs_field_t *alpha = CS_F_(alp_bl);

        cs_real_t *coefa_alp = alpha->bc_coeffs->a;
        cs_real_t *coefb_alp = alpha->bc_coeffs->b;
        cs_real_t *cofaf_alp = alpha->bc_coeffs->af;
        cs_real_t *cofbf_alp = alpha->bc_coeffs->bf;

        int *icodcl_alp = alpha->bc_coeffs->icodcl;
        cs_real_t *rcodcl1_alp = alpha->bc_coeffs->rcodcl1;
        cs_real_t *rcodcl2_alp = alpha->bc_coeffs->rcodcl2;
        cs_real_t *rcodcl3_alp = alpha->bc_coeffs->rcodcl3;

        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

          const cs_real_t distbf = b_dist[f_id];
          const cs_real_t hint = 1. / distbf;

          /* Dirichlet Boundary Condition
             ---------------------------- */

          if (icodcl_alp[f_id] == 1) {

            const cs_real_t pimp = rcodcl1_alp[f_id];
            const cs_real_t hext = rcodcl2_alp[f_id];

            cs_boundary_conditions_set_dirichlet_scalar(&coefa_alp[f_id],
                                                        &cofaf_alp[f_id],
                                                        &coefb_alp[f_id],
                                                        &cofbf_alp[f_id],
                                                        pimp,
                                                        hint,
                                                        hext);

          }

          /* Neumann Boundary Condition
             -------------------------- */

          else if (icodcl_alp[f_id] == 3) {

            const cs_real_t dimp = rcodcl3_alp[f_id];

            cs_boundary_conditions_set_neumann_scalar(&coefa_alp[f_id],
                                                      &cofaf_alp[f_id],
                                                      &coefb_alp[f_id],
                                                      &cofbf_alp[f_id],
                                                      dimp,
                                                      hint);
          }

          /* Convective Boundary Condition
             ----------------------------- */

          else if (icodcl_alp[f_id] == 2) {

            const cs_real_t pimp = rcodcl1_alp[f_id];
            const cs_real_t cfl  = rcodcl2_alp[f_id];

            cs_boundary_conditions_set_convective_outlet_scalar
              (&coefa_alp[f_id], &cofaf_alp[f_id],
               &coefb_alp[f_id], &cofbf_alp[f_id],
               pimp, cfl, hint);
          }

          /* Imposed value for the convection operator,
             imposed flux for diffusion
             ------------------------------------------ */

          else if (icodcl_alp[f_id] == 13) {

            const cs_real_t pimp = rcodcl1_alp[f_id];
            const cs_real_t dimp = rcodcl3_alp[f_id];

            cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
              (&coefa_alp[f_id], &cofaf_alp[f_id],
               &coefb_alp[f_id], &cofbf_alp[f_id],
               pimp, dimp);
          }
        }
      }
    }

    /* v2f type models (phi_bar and Bl-v2/k) */

    else if (itytur == 5) {

      /* k, epsilon  and phi */

      cs_field_t *v2f = NULL;
      cs_real_t sigma = 0.;

      for (int ii = 0; ii < 3; ii++) {

        if (ii == 1) {
          v2f = CS_F_(k);
          sigma = cs_field_get_key_double(v2f, ksigmas);
        }
        else if (ii == 2) {
          v2f = CS_F_(eps);
          sigma = cs_field_get_key_double(v2f, ksigmas);
        }
        else {
          v2f = CS_F_(phi);
          sigma = cs_field_get_key_double(v2f, ksigmas);
        }

        cs_real_t *coefa_v2f = v2f->bc_coeffs->a;
        cs_real_t *coefb_v2f = v2f->bc_coeffs->b;
        cs_real_t *cofaf_v2f = v2f->bc_coeffs->af;
        cs_real_t *cofbf_v2f = v2f->bc_coeffs->bf;

        int *icodcl_v2f = v2f->bc_coeffs->icodcl;
        cs_real_t *rcodcl1_v2f = v2f->bc_coeffs->rcodcl1;
        cs_real_t *rcodcl2_v2f = v2f->bc_coeffs->rcodcl2;
        cs_real_t *rcodcl3_v2f = v2f->bc_coeffs->rcodcl3;

        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

          const cs_lnum_t c_id = b_face_cells[f_id];

          /* physical properties */
          const cs_real_t visclc = viscl[c_id];
          const cs_real_t visctc = visct[c_id];

          /* geometric quantities */
          const cs_real_t distbf = b_dist[f_id];

          const cs_real_t hint = (visclc + visctc / sigma) / distbf;

          /* Dirichlet Boundary Condition
             ---------------------------- */

          if (icodcl_v2f[f_id] == 1) {

            const cs_real_t pimp = rcodcl1_v2f[f_id];
            const cs_real_t hext = rcodcl2_v2f[f_id];

            cs_boundary_conditions_set_dirichlet_scalar(&coefa_v2f[f_id],
                                                        &cofaf_v2f[f_id],
                                                        &coefb_v2f[f_id],
                                                        &cofbf_v2f[f_id],
                                                        pimp,
                                                        hint,
                                                        hext);
          }

          /* Neumann Boundary Condition
             -------------------------- */

          if (icodcl_v2f[f_id] == 3) {

            const cs_real_t dimp = rcodcl3_v2f[f_id];

            cs_boundary_conditions_set_neumann_scalar(&coefa_v2f[f_id],
                                                      &cofaf_v2f[f_id],
                                                      &coefb_v2f[f_id],
                                                      &cofbf_v2f[f_id],
                                                      dimp,
                                                      hint);
          }

          /* Convective Boundary Condition
             ----------------------------- */

          else if (icodcl_v2f[f_id] == 2) {

            const cs_real_t pimp = rcodcl1_v2f[f_id];
            const cs_real_t cfl  = rcodcl2_v2f[f_id];

            cs_boundary_conditions_set_convective_outlet_scalar
              (&coefa_v2f[f_id], &cofaf_v2f[f_id],
               &coefb_v2f[f_id], &cofbf_v2f[f_id],
               pimp, cfl, hint);
          }

          /* Imposed value for the convection operator,
             imposed flux for diffusion
             ----------------------------------------- */

          else if (icodcl_v2f[f_id] == 13) {

            const cs_real_t pimp = rcodcl1_v2f[f_id];
            const cs_real_t dimp = rcodcl3_v2f[f_id];

            cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
              (&coefa_v2f[f_id], &cofaf_v2f[f_id],
               &coefb_v2f[f_id], &cofbf_v2f[f_id],
               pimp, dimp);
          }
        }
      }

      if (iturb == CS_TURB_V2F_PHI) {

        /* FB */

        cs_field_t *f_bar = CS_F_(f_bar);

        cs_real_t *coefa_fb = f_bar->bc_coeffs->a;
        cs_real_t *coefb_fb = f_bar->bc_coeffs->b;
        cs_real_t *cofaf_fb = f_bar->bc_coeffs->af;
        cs_real_t *cofbf_fb = f_bar->bc_coeffs->bf;

        int *icodcl_fb = f_bar->bc_coeffs->icodcl;
        cs_real_t *rcodcl1_fb = f_bar->bc_coeffs->rcodcl1;
        cs_real_t *rcodcl2_fb = f_bar->bc_coeffs->rcodcl2;
        cs_real_t *rcodcl3_fb = f_bar->bc_coeffs->rcodcl3;

        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

          /* Physical Properties */
          const cs_real_t visclc = 1.0;

          /* Geometric quantities */
          const cs_real_t distbf = b_dist[f_id];

          const cs_real_t hint = visclc / distbf;

          /* Dirichlet Boundary Condition
             ---------------------------- */

          if (icodcl_fb[f_id] == 1) {

            const cs_real_t pimp = rcodcl1_fb[f_id];
            const cs_real_t hext = rcodcl2_fb[f_id];

            cs_boundary_conditions_set_dirichlet_scalar(&coefa_fb[f_id],
                                                        &cofaf_fb[f_id],
                                                        &coefb_fb[f_id],
                                                        &cofbf_fb[f_id],
                                                        pimp,
                                                        hint,
                                                        hext);
          }

          /* Neumann Boundary Condition
             -------------------------- */

          if (icodcl_fb[f_id] == 3) {

            const cs_real_t dimp = rcodcl3_fb[f_id];

            cs_boundary_conditions_set_neumann_scalar(&coefa_fb[f_id],
                                                      &cofaf_fb[f_id],
                                                      &coefb_fb[f_id],
                                                      &cofbf_fb[f_id],
                                                      dimp,
                                                      hint);
          }

          /* Convective Boundary Condition
             ------------------------------ */

          else if (icodcl_fb[f_id] == 2) {

            const cs_real_t pimp = rcodcl1_fb[f_id];
            const cs_real_t cfl  = rcodcl2_fb[f_id];

            cs_boundary_conditions_set_convective_outlet_scalar
              (&coefa_fb[f_id], &cofaf_fb[f_id],
               &coefb_fb[f_id],        &cofbf_fb[f_id],
               pimp, cfl, hint);
          }

          /* Imposed value for the convection operator,
             imposed flux for diffusion
             ------------------------------------------ */

          else if(icodcl_fb[f_id] == 13) {

            const cs_real_t pimp = rcodcl1_fb[f_id];
            const cs_real_t dimp = rcodcl3_fb[f_id];

            cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
              (&coefa_fb[f_id], &cofaf_fb[f_id],
               &coefb_fb[f_id], &cofbf_fb[f_id],
               pimp, dimp);
          }
        }
      }

      else if (iturb == CS_TURB_V2F_BL_V2K) {

        /* alpha */

        cs_field_t *alpha = CS_F_(alp_bl);

        cs_real_t *coefa_alp = alpha->bc_coeffs->a;
        cs_real_t *coefb_alp = alpha->bc_coeffs->b;
        cs_real_t *cofaf_alp = alpha->bc_coeffs->af;
        cs_real_t *cofbf_alp = alpha->bc_coeffs->bf;

        int *icodcl_alp = alpha->bc_coeffs->icodcl;
        cs_real_t *rcodcl1_alp = alpha->bc_coeffs->rcodcl1;
        cs_real_t *rcodcl2_alp = alpha->bc_coeffs->rcodcl2;
        cs_real_t *rcodcl3_alp = alpha->bc_coeffs->rcodcl3;

        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

          /* physical properties */
          const cs_real_t visclc = 1.0;

          /* geometric quantities */
          const cs_real_t distbf = b_dist[f_id];
          const cs_real_t hint = visclc / distbf;

          /* Dirichlet Boundary Condition
             ---------------------------- */

          if (icodcl_alp[f_id] == 1) {

            const cs_real_t pimp = rcodcl1_alp[f_id];
            const cs_real_t hext = rcodcl2_alp[f_id];

            cs_boundary_conditions_set_dirichlet_scalar(&coefa_alp[f_id],
                                                        &cofaf_alp[f_id],
                                                        &coefb_alp[f_id],
                                                        &cofbf_alp[f_id],
                                                        pimp,
                                                        hint,
                                                        hext);
          }

          /* Neumann Boundary Condition
             --------------------------- */

          else if (icodcl_alp[f_id] == 3) {

            const cs_real_t dimp = rcodcl3_alp[f_id];

            cs_boundary_conditions_set_neumann_scalar(&coefa_alp[f_id],
                                                      &cofaf_alp[f_id],
                                                      &coefb_alp[f_id],
                                                      &cofbf_alp[f_id],
                                                      dimp,
                                                      hint);
          }

          /* Convective Boundary Condition
             ----------------------------- */

          else if (icodcl_alp[f_id] == 2) {

            const cs_real_t pimp = rcodcl1_alp[f_id];
            const cs_real_t cfl  = rcodcl2_alp[f_id];

            cs_boundary_conditions_set_convective_outlet_scalar
              (&coefa_alp[f_id], &cofaf_alp[f_id],
               &coefb_alp[f_id], &cofbf_alp[f_id],
               pimp, cfl, hint);
          }

          /* Imposed value for the convection operator,
             imposed flux for diffusion
             ------------------------------------------ */

          else if (icodcl_alp[f_id] == 13) {

            const cs_real_t pimp = rcodcl1_alp[f_id];
            const cs_real_t dimp = rcodcl3_alp[f_id];

            cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
              (&coefa_alp[f_id], &cofaf_alp[f_id],
               &coefb_alp[f_id], &cofbf_alp[f_id],
               pimp, dimp);
          }
        }
      }
    }

    /* Spalart Allmaras */

    else if (iturb == CS_TURB_SPALART_ALLMARAS) {

      cs_field_t *nusa = CS_F_(nusa);

      cs_real_t *coefa_nusa = nusa->bc_coeffs->a;
      cs_real_t *coefb_nusa = nusa->bc_coeffs->b;
      cs_real_t *cofaf_nusa = nusa->bc_coeffs->af;
      cs_real_t *cofbf_nusa = nusa->bc_coeffs->bf;

      int *icodcl_nusa = nusa->bc_coeffs->icodcl;
      cs_real_t *rcodcl1_nusa = nusa->bc_coeffs->rcodcl1;
      cs_real_t *rcodcl2_nusa = nusa->bc_coeffs->rcodcl2;
      cs_real_t *rcodcl3_nusa = nusa->bc_coeffs->rcodcl3;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

        const cs_lnum_t c_id = b_face_cells[f_id];

        /* physical properties */
        const cs_real_t visclc = viscl[c_id];

        /* geometric quantities */
        const cs_real_t distbf = b_dist[f_id];
        const cs_real_t hint = visclc / distbf;

        /* Dirichlet Boundary Condition
           ---------------------------- */

        if (icodcl_nusa[f_id] == 1) {

          const cs_real_t pimp = rcodcl1_nusa[f_id];
          const cs_real_t hext = rcodcl2_nusa[f_id];

          cs_boundary_conditions_set_dirichlet_scalar(&coefa_nusa[f_id],
                                                      &cofaf_nusa[f_id],
                                                      &coefb_nusa[f_id],
                                                      &cofbf_nusa[f_id],
                                                      pimp,
                                                      hint,
                                                      hext);
        }

        /* Neumann Boundary Condition
           -------------------------- */

        if (icodcl_nusa[f_id] == 3) {

          const cs_real_t dimp = rcodcl3_nusa[f_id];

          cs_boundary_conditions_set_neumann_scalar(&coefa_nusa[f_id],
                                                    &cofaf_nusa[f_id],
                                                    &coefb_nusa[f_id],
                                                    &cofbf_nusa[f_id],
                                                    dimp,
                                                    hint);
        }

        /* Convective Boundary Condition
           ----------------------------- */

        else if (icodcl_nusa[f_id] == 2) {

          const cs_real_t pimp = rcodcl1_nusa[f_id];
          const cs_real_t cfl  = rcodcl2_nusa[f_id];

          cs_boundary_conditions_set_convective_outlet_scalar
            (&coefa_nusa[f_id], &cofaf_nusa[f_id],
             &coefb_nusa[f_id], &cofbf_nusa[f_id],
             pimp, cfl, hint);
        }

        /* Imposed value for the convection operator,
           imposed flux for diffusion
           ----------------------------------------- */

        else if (icodcl_nusa[f_id] == 13) {

          const cs_real_t pimp = rcodcl1_nusa[f_id];
          const cs_real_t dimp = rcodcl3_nusa[f_id];

          cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
            (&coefa_nusa[f_id],        &cofaf_nusa[f_id],
             &coefb_nusa[f_id], &cofbf_nusa[f_id],
             pimp, dimp);
        }
      }
    }
  }

  /*--------------------------------------------------------------------------
   * 13) Other scalars (except variances):
   *     Dirichlet and Neumann and convective outlet
   *--------------------------------------------------------------------------*/

  {
    const cs_real_t *cpro_cv = NULL, *cpro_cp = NULL;

    if (fluid_props->icp >= 0)
      cpro_cp = CS_F_(cp)->val;

    cs_field_t *f_id_cv = cs_field_by_name_try("isobaric_heat_capacity");
    if (f_id_cv != NULL)
      cpro_cv = f_id_cv->val;

    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0
        && fluid_props->icv >= 0)
      cpro_cv = cs_field_by_id(fluid_props->icv)->val;

    for (int ii = 0; ii < n_fields; ii++) {

      cs_field_t *f_scal = cs_field_by_id(ii);

      if (!(f_scal->type & CS_FIELD_VARIABLE))
        continue;
      if (cs_field_get_key_int(f_scal, keysca) <= 0)
        continue;

      cs_lnum_t isvhbl = 0;
      if (ii == isvhb)
        isvhbl = isvhb;

      const int ifcvsl = cs_field_get_key_int(f_scal, kivisl);
      const cs_real_t *viscls = NULL;
      if (ifcvsl >= 0)
        viscls = cs_field_by_id(ifcvsl)->val;

      /* Get the turbulent flux model for the scalar */
      int turb_flux_model = cs_field_get_key_int(f_scal, kturt);
      int turb_flux_model_type = turb_flux_model / 10;

      /* --- Indicateur de prise en compte de Cp ou non
         (selon si le scalaire (scalaire associe pour une fluctuation)
         doit etre ou non traite comme une temperature)
         Si le scalaire est une variance et que le
         scalaire associe n'est pas resolu, on suppose alors qu'il
         doit etre traite comme un scalaire passif (defaut IHCP = 0)*/
      cs_lnum_t ihcp = 0;
      const int kscavr = cs_field_key_id("first_moment_id");
      const int iscavr = cs_field_get_key_int(f_scal, kscavr);

      /* Reference diffusivity */
      const int kvisl0 = cs_field_key_id("diffusivity_ref");
      cs_real_t visls_0 = 0.;
      cs_lnum_t iscacp = 0;
      cs_field_t *f = NULL;

      if (iscavr > 0) {
        f = cs_field_by_id(iscavr);
        visls_0 = cs_field_get_key_double(f, kvisl0);
        iscacp  = cs_field_get_key_int(cs_field_by_id(iscavr), kscacp);
      }
      else {
        f = f_scal;
        visls_0 = cs_field_get_key_double(f_scal, kvisl0);
        iscacp  = cs_field_get_key_int(f_scal, kscacp);
      }

      if (iscacp == 1) {
        if (fluid_props->icp >= 0)
          ihcp = 2;
        else
          ihcp = 1;
      }
      else if (iscacp == 2) {
        if (fluid_props->icp >= 0)
          ihcp = 4;
        else
          ihcp = 3;
      }

      cs_equation_param_t *eqp_scal = cs_field_get_equation_param(f_scal);
      cs_field_t *f_vis = NULL;
      cs_real_t ctheta = 0.;
      cs_real_6_t *visten = NULL;
      cs_field_t  *f_a_t_visc = NULL;

      if ((eqp_scal->idften & CS_ANISOTROPIC_DIFFUSION)
          || turb_flux_model_type == 3) {

        if (iturb != CS_TURB_RIJ_EPSILON_EBRSM || turb_flux_model_type == 3) {
          f_a_t_visc = cs_field_by_name("anisotropic_turbulent_viscosity");
          visten = (cs_real_6_t *)f_a_t_visc->val;
        }
        else { /* EBRSM and (GGDH or AFM) */
          f_vis = cs_field_by_name("anisotropic_turbulent_viscosity_scalar");

          visten = (cs_real_6_t *)f_vis->val;
        }
        const int kctheta = cs_field_key_id("turbulent_flux_ctheta");
        ctheta = cs_field_get_key_double(f_scal, kctheta);
      }

      cs_real_t turb_schmidt = cs_field_get_key_double(f_scal, ksigmas);

      /* Get boundary value (for post-processing) */
      int b_f_id = cs_field_get_key_int(f_scal, kbfid);
      cs_lnum_t f_dim = f->dim;

      /* Scalar transported quantity */
      if (f_dim == 1) {

        cs_real_t *coefa_sc = f_scal->bc_coeffs->a;
        cs_real_t *coefb_sc = f_scal->bc_coeffs->b;
        cs_real_t *cofaf_sc = f_scal->bc_coeffs->af;
        cs_real_t *cofbf_sc = f_scal->bc_coeffs->bf;

        int *icodcl_sc = f_scal->bc_coeffs->icodcl;
        cs_real_t *rcodcl1_sc = f_scal->bc_coeffs->rcodcl1;
        cs_real_t *rcodcl2_sc = f_scal->bc_coeffs->rcodcl2;
        cs_real_t *rcodcl3_sc = f_scal->bc_coeffs->rcodcl3;

        cs_equation_param_t *eqp_sc = cs_field_get_equation_param(f_scal);
        cs_field_t  *f_scal_b = NULL;
        cs_real_t *bvar_s = NULL;

        if (b_f_id > -1)
          f_scal_b = cs_field_by_id(b_f_id);
        else {
          /* if thermal variable has no boundary but temperature does, use it */
          if (f_scal == f_th && f_scal == CS_F_(h))
            f_scal_b = cs_field_by_name_try("boundary_temperature");
        }

        if (f_scal_b != NULL)
          bvar_s = f_scal_b->val;

        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

          const cs_lnum_t c_id = b_face_cells[f_id];
          const cs_real_t surf = b_face_surf[f_id];

          /* Physical Properties */
          const cs_real_t visctc = visct[c_id];
          const cs_real_t visclc = viscl[c_id];
          cs_real_t hint = 0.0;

          /* Geometric quantities */
          const cs_real_t distbf = b_dist[f_id];
          cs_real_t visci[3][3], dist[3];
          const cs_real_t *n = b_face_normal[f_id];

          dist[0] = b_face_cog[f_id][0] - cell_cen[c_id][0];
          dist[1] = b_face_cog[f_id][1] - cell_cen[c_id][1];
          dist[2] = b_face_cog[f_id][2] - cell_cen[c_id][2];

          /* --- Prise en compte de Cp ou CV
             (dans le Cas compressible ihcp=0) */

          cs_real_t cpp = 1.0;
          if (ihcp == 1)
            cpp = cp0;
          else if (ihcp == 2)
            cpp = cpro_cp[c_id];
          else if (ihcp >= 3)
            cpp = cpro_cv[c_id];

          cs_real_t rkl;
          /* --- Viscosite variable ou non */
          if (ifcvsl < 0)
            rkl = visls_0;
          else
            rkl = viscls[c_id];

          /* Scalar diffusivity */
          if (eqp_sc->idften & CS_ISOTROPIC_DIFFUSION)
            hint = (rkl+eqp_sc->idifft*cpp*visctc/turb_schmidt)/distbf;

          /* Symmetric tensor diffusivity */
          else if (eqp_sc->idften & CS_ANISOTROPIC_DIFFUSION) {
            const cs_real_t temp = eqp_sc->idifft*cpp*ctheta/cs_turb_csrij;
            visci[0][0] = rkl + temp*visten[c_id][0];
            visci[1][1] = rkl + temp*visten[c_id][1];
            visci[2][2] = rkl + temp*visten[c_id][2];
            visci[0][1] =       temp*visten[c_id][3];
            visci[1][0] =       temp*visten[c_id][3];
            visci[1][2] =       temp*visten[c_id][4];
            visci[2][1] =       temp*visten[c_id][4];
            visci[0][2] =       temp*visten[c_id][5];
            visci[2][0] =       temp*visten[c_id][5];

            /* ||Ki.S||^2 */
            const cs_real_t viscis = cs_math_pow2(  visci[0][0]*n[0]
                                                 + visci[1][0]*n[1]
                                                 + visci[2][0]*n[2])
                                  + cs_math_pow2(  visci[0][1]*n[0]
                                                 + visci[1][1]*n[1]
                                                 + visci[2][1]*n[2])
                                  + cs_math_pow2(  visci[0][2]*n[0]
                                                 + visci[1][2]*n[1]
                                                 + visci[2][2]*n[2]);

            /* if.ki.s */
            cs_real_t fikis
              = (  cs_math_3_dot_product(dist, visci[0]) * n[0]
                 + cs_math_3_dot_product(dist, visci[1]) * n[1]
                 + cs_math_3_dot_product(dist, visci[2]) * n[2]);

            const cs_real_t distfi = b_dist[f_id];

            /* Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
               NB: eps =1.d-1 must be consistent with vitens.f90 */
            fikis = cs_math_fmax(fikis, 1.e-1*sqrt(viscis)*distfi);

            hint = viscis / surf / fikis;

          }

          /* Dirichlet Boundary Condition
             ---------------------------- */

          if (icodcl_sc[f_id] == 1) {

            const cs_real_t pimp = rcodcl1_sc[f_id];
            const cs_real_t hext = rcodcl2_sc[f_id];

            cs_boundary_conditions_set_dirichlet_scalar(&coefa_sc[f_id],
                                                        &cofaf_sc[f_id],
                                                        &coefb_sc[f_id],
                                                        &cofbf_sc[f_id],
                                                        pimp,
                                                        hint,
                                                        hext);

            /* Store boundary value */
            if (f_scal_b != NULL)
              bvar_s[f_id] = coefa_sc[f_id] + coefb_sc[f_id] * bvar_s[f_id];
          }

          /* Neumann Boundary Conditions
             ---------------------------- */

          if (icodcl_sc[f_id] == 3) {

            const cs_real_t dimp = rcodcl3_sc[f_id];

            cs_boundary_conditions_set_neumann_scalar(&coefa_sc[f_id],
                                                      &cofaf_sc[f_id],
                                                      &coefb_sc[f_id],
                                                      &cofbf_sc[f_id],
                                                      dimp,
                                                      hint);

            /* Store boundary value only for faces
               for which it was not previously computed
               in clptur.f90 */
            if (icodcl_vel[f_id] != 5) {
              if (f_scal_b != NULL)
                bvar_s[f_id] = coefa_sc[f_id] + coefb_sc[f_id] * bvar_s[f_id];
            }
          }

          /* Convective Boundary Condition
             ------------------------------ */

          else if (icodcl_sc[f_id] == 2 && iterns <= 1) {

            const cs_real_t pimp = rcodcl1_sc[f_id];
            const cs_real_t cfl  = rcodcl2_sc[f_id];

            cs_boundary_conditions_set_convective_outlet_scalar
              (&coefa_sc[f_id], &cofaf_sc[f_id],
               &coefb_sc[f_id], &cofbf_sc[f_id],
               pimp, cfl, hint);

            /* Store boundary value */
            if (f_scal_b != NULL)
              bvar_s[f_id] = coefa_sc[f_id] + coefb_sc[f_id] * bvar_s[f_id];
          }

          /* Set total flux as a Robin condition
             ----------------------------------- */

          else if (icodcl_sc[f_id] == 12) {

            const cs_real_t hext = rcodcl2_sc[f_id];
            const cs_real_t dimp = rcodcl3_sc[f_id];

            cs_boundary_conditions_set_total_flux(&coefa_sc[f_id],
                                                  &cofaf_sc[f_id],
                                                  &coefb_sc[f_id],
                                                  &cofbf_sc[f_id],
                                                  hext,
                                                  dimp);

            /* Store boundary value */
            if (f_scal_b != NULL)
              bvar_s[f_id] = coefa_sc[f_id] + coefb_sc[f_id] * bvar_s[f_id];
          }

          /* Imposed value for the convection operator,
             imposed flux for diffusion
             ------------------------------------------ */

          else if (icodcl_sc[f_id] == 13) {

            const cs_real_t pimp = rcodcl1_sc[f_id];
            const cs_real_t dimp = rcodcl3_sc[f_id];

            cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
              (&coefa_sc[f_id], &cofaf_sc[f_id],
               &coefb_sc[f_id], &cofbf_sc[f_id],
               pimp, dimp);

            /* Store boundary value */
            if (f_scal_b != NULL)
              bvar_s[f_id] = coefa_sc[f_id] + coefb_sc[f_id] * bvar_s[f_id];
          }

          /* Store the thermal exchange coefficient
             (conversion in case of energy or enthalpy)
             the exchange coefficient is in W/(m2 K)
             Useful for thermal coupling or radiative transfer */

          if (icodcl_sc[f_id] == 1 || icodcl_sc[f_id] == 3) {
            cs_real_t exchange_coef = 0.;
            if ((cs_glob_rad_transfer_params->type >= 1 &&
                 f_th == f_scal) || isvhbl > 0) {

              /* Enthalpy */
              if (thermal_variable == CS_THERMAL_MODEL_ENTHALPY) {
                /* If Cp is variable */
                if (fluid_props->icp >= 0)
                  exchange_coef = hint*cpro_cp[c_id];
                else
                  exchange_coef = hint*cp0;
              }

              /* Total energy (compressible module) */
              else if (thermal_variable == CS_THERMAL_MODEL_TOTAL_ENERGY) {
                /* If Cv is variable */
                if (fluid_props->icv >= 0)
                  exchange_coef = hint*cpro_cv[c_id];
                else
                  exchange_coef = hint*fluid_props->cv0;
              }

              /* Temperature */
              else if (iscacp > 0)
                exchange_coef = hint;
            }

            /* ---> Thermal coupling, store hint = lambda/d */
            if (isvhbl > 0)
              hbord[f_id] = exchange_coef;

            /* ---> Radiative transfer */
            if (cs_glob_rad_transfer_params->type >= 1 && f_th == f_scal) {
              bhconv[f_id] = exchange_coef;

              /* The outgoing flux is stored (Q = h(Ti'-Tp): negative if
                 gain for the fluid) in W/m2 */
              bfconv[f_id] = cofaf_sc[f_id] + cofbf_sc[f_id] * theipb[f_id];
            }
          }

          /* Thermal heat flux boundary conditions */
          if (turb_flux_model_type == 3) {

            cs_field_t *f_tf
              = cs_field_by_composite_name_try(f_scal->name, "turbulent_flux");

            cs_real_3_t  *coefa_tf = (cs_real_3_t  *)f_tf->bc_coeffs->a;
            cs_real_33_t *coefb_tf = (cs_real_33_t *)f_tf->bc_coeffs->b;
            cs_real_3_t  *cofaf_tf = (cs_real_3_t  *)f_tf->bc_coeffs->af;
            cs_real_33_t *cofbf_tf = (cs_real_33_t *)f_tf->bc_coeffs->bf;
            cs_real_3_t  *cofar_tf = (cs_real_3_t  *)f_tf->bc_coeffs->ad;
            cs_real_33_t *cofbr_tf = (cs_real_33_t *)f_tf->bc_coeffs->bd;

            int *icodcl_tf = f_tf->bc_coeffs->icodcl;
            cs_real_t *rcodcl1_tf = f_tf->bc_coeffs->rcodcl1;
            cs_real_t *rcodcl2_tf = f_tf->bc_coeffs->rcodcl2;
            cs_real_t *rcodcl3_tf = f_tf->bc_coeffs->rcodcl3;

            cs_real_3_t pimpv = {0., 0., 0.};
            cs_real_3_t hextv = {0., 0., 0.};
            cs_real_3_t qimpv = {0., 0., 0.};
            cs_real_3_t cflv  = {0., 0., 0.};
            cs_real_6_t hintt = {0., 0. ,0., 0., 0., 0.};

            if (ifcvsl < 0)
              rkl = visls_0/cpp;
            else
              rkl = viscls[c_id]/cpp;

            hintt[0] = 0.5*(visclc+rkl)/distbf
                     + visten[c_id][0]*ctheta/distbf/cs_turb_csrij;

            hintt[1] = 0.5*(visclc+rkl)/distbf
                     + visten[c_id][1]*ctheta/distbf/cs_turb_csrij;

            hintt[2] = 0.5*(visclc+rkl)/distbf
                     + visten[c_id][2]*ctheta/distbf/cs_turb_csrij;

            hintt[3] = visten[c_id][3]*ctheta/distbf/cs_turb_csrij;
            hintt[4] = visten[c_id][4]*ctheta/distbf/cs_turb_csrij;
            hintt[5] = visten[c_id][5]*ctheta/distbf/cs_turb_csrij;

            /* Dirichlet Boundary Condition
               ---------------------------- */

            if (icodcl_tf[f_id] == 1) {

              for (cs_lnum_t k = 0; k < 3; k++)
                pimpv[k] = rcodcl1_tf[n_b_faces*k + f_id];

              for (cs_lnum_t k = 0; k < 3; k++)
                hextv[k] = rcodcl2_tf[n_b_faces*k + f_id];

              cs_boundary_conditions_set_dirichlet_vector_aniso
                (coefa_tf[f_id], cofaf_tf[f_id],
                 coefb_tf[f_id], cofbf_tf[f_id],
                 pimpv, hintt, hextv);

              /* Boundary conditions for thermal transport equation */
              for (int isou = 0; isou < 3; isou++) {
                cofar_tf[f_id][isou] = coefa_tf[f_id][isou];
                for (int jsou = 0; jsou < 3; jsou++)
                  cofbr_tf[f_id][isou][jsou] = coefb_tf[f_id][isou][jsou];
              }
            }

            /* Neumann Boundary Conditions
               ---------------------------- */

            else if (icodcl_tf[f_id] == 3) {

              for (cs_lnum_t k = 0; k < 3; k++)
                qimpv[k] = rcodcl3_tf[n_b_faces*k + f_id];

              cs_boundary_conditions_set_neumann_vector_aniso
                (coefa_tf[f_id], cofaf_tf[f_id],
                 coefb_tf[f_id], cofbf_tf[f_id],
                 qimpv, hintt);

              /* Boundary conditions for thermal transport equation */
              for (int isou = 0; isou < 3; isou++) {
                cofar_tf[f_id][isou] = coefa_tf[f_id][isou];
                for (int jsou = 0; jsou < 3; jsou++)
                  cofbr_tf[f_id][isou][jsou] = coefb_tf[f_id][isou][jsou];
              }
            }

            /* Convective Boundary Conditions
               ------------------------------- */

            else if (icodcl_tf[f_id] == 2) {

              for (cs_lnum_t k = 0; k < 3; k++)
                pimpv[k] = rcodcl1_tf[n_b_faces*k + f_id];

              for (cs_lnum_t k = 0; k < 3; k++)
                cflv[k] = rcodcl2_tf[n_b_faces*k + f_id];

              cs_boundary_conditions_set_convective_outlet_vector_aniso
                (coefa_tf[f_id], cofaf_tf[f_id],
                 coefb_tf[f_id], cofbf_tf[f_id],
                 pimpv,        cflv, hintt);

              /* Boundary conditions for thermal transport equation */
              for (int isou = 0; isou < 3; isou++) {
                cofar_tf[f_id][isou] = coefa_tf[f_id][isou];
                for (int jsou = 0; jsou < 3; jsou++)
                  cofbr_tf[f_id][isou][jsou] = coefb_tf[f_id][isou][jsou];
              }
            }
          }
        } /* end f_id < n_b_faces */
      }/* end if f->dim = 1 */

      /* Vector transported quantity (dimension may be greater than 3) */
      else {

        cs_real_3_t *bvar_v = NULL;
        if (b_f_id >= 0)
          bvar_v = (cs_real_3_t *)cs_field_by_id(b_f_id)->val;

        cs_field_t *vtq  = f_scal;

        cs_real_3_t  *coefa_vtq = (cs_real_3_t  *)vtq->bc_coeffs->a;
        cs_real_33_t *coefb_vtq = (cs_real_33_t *)vtq->bc_coeffs->b;
        cs_real_3_t  *cofaf_vtq = (cs_real_3_t  *)vtq->bc_coeffs->af;
        cs_real_33_t *cofbf_vtq = (cs_real_33_t *)vtq->bc_coeffs->bf;

        int *icodcl_vtq = vtq->bc_coeffs->icodcl;
        cs_real_t *rcodcl1_vtq = vtq->bc_coeffs->rcodcl1;
        cs_real_t *rcodcl2_vtq = vtq->bc_coeffs->rcodcl2;
        cs_real_t *rcodcl3_vtq = vtq->bc_coeffs->rcodcl3;

        cs_equation_param_t *eqp_vtq = cs_field_get_equation_param(vtq);

        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

          const cs_lnum_t c_id  = b_face_cells[f_id];

          /* physical properties */
          const cs_real_t visctc = visct[c_id];

          /* geometric quantities */
          const cs_real_t distbf = b_dist[f_id];

          cs_real_t hintt[6]   = {0., 0., 0., 0., 0., 0.};
          cs_real_t pimpv[3]   = {0., 0., 0};
          cs_real_t hextv[3]   = {0., 0., 0};
          cs_real_t qimpv[3]   = {0., 0., 0};
          cs_real_t cflv[3]    = {0., 0., 0};
          cs_real_t b_pvari[3] = {0., 0., 0};

          /* Account for Cp or CV
             (in the compressible case ihcp=0) */

          cs_real_t cpp = 1.;
          if (ihcp == 1)
            cpp = cp0;
          else if (ihcp == 2)
            cpp = cpro_cp[c_id];
          else if (ihcp == 3)
            cpp = cp0 - fluid_props->r_pg_cnst; /* TODO: humid air */
          else if (ihcp == 4)
            cpp = cpro_cp[c_id] - fluid_props->r_pg_cnst;

          /* variable or constant viscosity */
          cs_real_t rkl, hint;
          if (ifcvsl < 0)
            rkl = visls_0;
          else
            rkl = viscls[c_id];

          /* Scalar diffusivity */
          if (eqp_vtq->idften & CS_ISOTROPIC_DIFFUSION) {
            /* FIXME */
            hint = (rkl+eqp_vtq->idifft*cpp*visctc/turb_schmidt)/distbf;

            hintt[0] = hint;
            hintt[1] = hint;
            hintt[2] = hint;
            hintt[3] = 0.0;
            hintt[4] = 0.0;
            hintt[5] = 0.0;
          }

          /* Symmetric tensor diffusivity */
          else if (eqp_vtq->idften & CS_ANISOTROPIC_DIFFUSION) {
            const cs_real_t temp = eqp_vtq->idifft*cpp*ctheta/cs_turb_csrij;
            hintt[0] = (rkl + temp*visten[c_id][0])/distbf;
            hintt[1] = (rkl + temp*visten[c_id][1])/distbf;
            hintt[2] = (rkl + temp*visten[c_id][2])/distbf;
            hintt[3] =        temp*visten[c_id][3] /distbf;
            hintt[4] =        temp*visten[c_id][4] /distbf;
            hintt[5] =        temp*visten[c_id][5] /distbf;
          }

          /* Dirichlet Boundary Condition
             ---------------------------- */

          if (icodcl_vtq[f_id] == 1) {

            for (cs_lnum_t k = 0; k < 3; k++)
              pimpv[k] = rcodcl1_vtq[n_b_faces*k + f_id];

            for (cs_lnum_t k = 0; k < 3; k++)
              hextv[k] = rcodcl2_vtq[n_b_faces*k + f_id];

            cs_boundary_conditions_set_dirichlet_vector_aniso
              (coefa_vtq[f_id], cofaf_vtq[f_id],
               coefb_vtq[f_id], cofbf_vtq[f_id],
               pimpv, hintt, hextv);

            /* Store boundary value */
            if (b_f_id >= 0) {
              /* B_ij. Pj(I) */
              cs_math_33_3_product(coefb_vtq[f_id], bvar_v[f_id], b_pvari);
              for (int isou = 0; isou < 3; isou++)
                bvar_v[f_id][isou] = coefa_vtq[f_id][isou] + b_pvari[isou];
            }
          }

          /* Neumann Boundary Condition
             -------------------------- */

          if (icodcl_vtq[f_id] == 3) {

            for (cs_lnum_t k = 0; k < 3; k++)
              qimpv[k] = rcodcl3_vtq[n_b_faces*k + f_id];

            cs_boundary_conditions_set_neumann_vector_aniso
              (coefa_vtq[f_id], cofaf_vtq[f_id],
               coefb_vtq[f_id], cofbf_vtq[f_id],
               qimpv, hintt);

            /* Store boundary value */
            if (b_f_id >= 0) {
              // B_ij. Pj(I)
              cs_math_33_3_product(coefb_vtq[f_id], bvar_v[f_id], b_pvari);
              for (int isou = 0; isou < 3; isou++)
                bvar_v[f_id][isou] = coefa_vtq[f_id][isou] + b_pvari[isou];
            }
          }

          /* Convective Boundary Condition
             ----------------------------- */

          else if (icodcl_vtq[f_id] == 2) {

            for (cs_lnum_t k = 0; k < 3; k++)
              pimpv[k] = rcodcl1_vtq[n_b_faces*k + f_id];

            for (cs_lnum_t k = 0; k < 3; k++)
              cflv[k] = rcodcl2_vtq[n_b_faces*k + f_id];

            cs_boundary_conditions_set_convective_outlet_vector_aniso
              (coefa_vtq[f_id], cofaf_vtq[f_id],
               coefb_vtq[f_id], cofbf_vtq[f_id],
               pimpv, cflv, hintt);

            /* Store boundary value */
            if (b_f_id >= 0) {
              /* B_ij. Pj(I) */
              cs_math_33_3_product(coefb_vtq[f_id], bvar_v[f_id], b_pvari);
              for (int isou = 0; isou < 3; isou++)
                bvar_v[f_id][isou] = coefa_vtq[f_id][isou] + b_pvari[isou];
            }
          }

          /* Imposed value for the convection operator,
             imposed flux for diffusion
             ---------------------------------------- */

          else if (icodcl_vtq[f_id] == 13) {

            for (cs_lnum_t k = 0; k < 3; k++)
              pimpv[k] = rcodcl1_vtq[n_b_faces*k + f_id];

            for (cs_lnum_t k = 0; k < 3; k++)
              qimpv[k] = rcodcl3_vtq[n_b_faces*k + f_id];

            cs_boundary_conditions_set_dirichlet_conv_neumann_diff_vector
              (coefa_vtq[f_id], cofaf_vtq[f_id],
               coefb_vtq[f_id], cofbf_vtq[f_id],
               pimpv, qimpv);

            /* Store boundary value */
            if (b_f_id >= 0) {
              /* B_ij. Pj(I) */
              cs_math_33_3_product(coefb_vtq[f_id], bvar_v[f_id], b_pvari);
              for (int isou = 0; isou < 3; isou++)
                bvar_v[f_id][isou] = coefa_vtq[f_id][isou] + b_pvari[isou];
            }
          }

          /* Convective boundary for Marangoni effects
             (generalized symmetry condition)
             ------------------------------------------ */

          else if (icodcl_vtq[f_id] == 14) {

            for (cs_lnum_t k = 0; k < 3; k++)
              pimpv[k] = rcodcl1_vtq[n_b_faces*k + f_id];

            for (cs_lnum_t k = 0; k < 3; k++)
              qimpv[k] = rcodcl3_vtq[n_b_faces*k + f_id];

            /* Coupled solving of the velocity components */

            cs_boundary_conditions_set_generalized_sym_vector_aniso
              (coefa_vtq[f_id], cofaf_vtq[f_id],
               coefb_vtq[f_id], cofbf_vtq[f_id],
               pimpv, qimpv, hintt, b_face_u_normal[f_id]);

            /* Store boundary value */
            if (b_f_id >= 0) {
              /* B_ij. Pj(I) */
              cs_math_33_3_product(coefb_vtq[f_id], bvar_v[f_id], b_pvari);
              for (int isou = 0; isou < 3; isou++)
                bvar_v[f_id][isou] = coefa_vtq[f_id][isou] + b_pvari[isou];
            }
          }

          /* Neumann on the normal component,
             Dirichlet on tangential components
             ---------------------------------- */

          else if (icodcl_vtq[f_id] == 11) {

            /* Dirichlet to impose on the tangential components */
            for (cs_lnum_t k = 0; k < 3; k++)
              pimpv[k] = rcodcl1_vtq[n_b_faces*k + f_id];

            /* Flux to impose on the normal component */
            for (cs_lnum_t k = 0; k < 3; k++)
              qimpv[k] = rcodcl3_vtq[n_b_faces*k + f_id];

            /* coupled solving of the velocity components */

            cs_boundary_conditions_set_generalized_dirichlet_vector_aniso
              (coefa_vtq[f_id], cofaf_vtq[f_id],
               coefb_vtq[f_id], cofbf_vtq[f_id],
               pimpv, qimpv, hintt, b_face_u_normal[f_id]);

            /* Store boundary value */
            if (b_f_id >= 0) {
              // B_ij. Pj(I)
              cs_math_33_3_product(coefb_vtq[f_id], bvar_v[f_id], b_pvari);
              for (int isou = 0; isou < 3; isou++)
                bvar_v[f_id][isou] = coefa_vtq[f_id][isou] + b_pvari[isou];
            }
          }
        } /* End of loop on faces */
      } /* End of vector transported quantities */

      /* EB-GGDH/AFM/DFM alpha boundary conditions */
      if (   turb_flux_model == 11
          || turb_flux_model == 21
          || turb_flux_model == 31) {

        cs_field_t *f_al = cs_field_by_composite_name_try(f_scal->name, "alpha");
        cs_real_t *coefa_al = f_al->bc_coeffs->a;
        cs_real_t *coefb_al = f_al->bc_coeffs->b;
        cs_real_t *cofaf_al = f_al->bc_coeffs->af;
        cs_real_t *cofbf_al = f_al->bc_coeffs->bf;

        int *icodcl_al = f_al->bc_coeffs->icodcl;
        cs_real_t *rcodcl1_al = f_al->bc_coeffs->rcodcl1;
        cs_real_t *rcodcl2_al = f_al->bc_coeffs->rcodcl2;
        cs_real_t *rcodcl3_al = f_al->bc_coeffs->rcodcl3;

        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

          const cs_real_t distbf = b_dist[f_id];
          const cs_real_t hint = 1. / distbf;

          /* Dirichlet Boundary Condition
             ---------------------------- */

          if (icodcl_al[f_id] == 1) {

            const cs_real_t pimp = rcodcl1_al[f_id];
            const cs_real_t hext = rcodcl2_al[f_id];

            cs_boundary_conditions_set_dirichlet_scalar(&coefa_al[f_id],
                                                        &cofaf_al[f_id],
                                                        &coefb_al[f_id],
                                                        &cofbf_al[f_id],
                                                        pimp,
                                                        hint,
                                                        hext);

          }

          /* Neumann Boundary Condition
             -------------------------- */

          if (icodcl_al[f_id] == 3) {

            const cs_real_t dimp = rcodcl3_al[f_id];

            cs_boundary_conditions_set_neumann_scalar(&coefa_al[f_id],
                                                      &cofaf_al[f_id],
                                                      &coefb_al[f_id],
                                                      &cofbf_al[f_id],
                                                      dimp,
                                                      hint);
          }

          /* Radiative Boundary Condition
             ---------------------------- */

          else if (icodcl_al[f_id] == 2) {

            const cs_real_t pimp = rcodcl1_al[f_id];
            const cs_real_t cfl  = rcodcl2_al[f_id];

            cs_boundary_conditions_set_convective_outlet_scalar
              (&coefa_al[f_id], &cofaf_al[f_id],
               &coefb_al[f_id], &cofbf_al[f_id],
               pimp, cfl, hint);
          }

          /* Imposed value for the convection operator,
             imposed flux for diffusion
             ------------------------------------------- */

          else if (icodcl_al[f_id] == 13) {

            const cs_real_t pimp = rcodcl1_al[f_id];
            const cs_real_t dimp = rcodcl3_al[f_id];

            cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar
              (&coefa_al[f_id], &cofaf_al[f_id],
               &coefb_al[f_id], &cofbf_al[f_id],
               pimp, dimp);

          }
        } /* End of loop on face */
      }
    } /* End of loop on scalars */
  } /* End other scalars */

  /*--------------------------------------------------------------------------
   * 14) Mesh velocity (ALE module):
   *     Dirichlet and Neumann and convective outlet
   *--------------------------------------------------------------------------*/

  if (cs_glob_ale == CS_ALE_LEGACY) {

    cs_field_t *m_vel = cs_field_by_name("mesh_velocity");

    cs_real_3_t  *claale = (cs_real_3_t  *)m_vel->bc_coeffs->a;
    cs_real_33_t *clbale = (cs_real_33_t *)m_vel->bc_coeffs->b;
    cs_real_3_t  *cfaale = (cs_real_3_t  *)m_vel->bc_coeffs->af;
    cs_real_33_t *cfbale = (cs_real_33_t *)m_vel->bc_coeffs->bf;

    int *icodcl_displ = m_vel->bc_coeffs->icodcl;
    cs_real_t *rcodcl1_displ = m_vel->bc_coeffs->rcodcl1;
    cs_real_t *rcodcl2_displ = m_vel->bc_coeffs->rcodcl2;
    cs_real_t *rcodcl3_displ = m_vel->bc_coeffs->rcodcl3;

    cs_equation_param_t *eqp_displ = cs_field_get_equation_param(m_vel);

    const cs_real_t   *cpro_visma_s = NULL;
    const cs_real_6_t *cpro_visma_v = NULL;

    if (eqp_displ->idften & CS_ISOTROPIC_DIFFUSION)
      cpro_visma_s = CS_F_(vism)->val;
    else if (eqp_displ->idften & CS_ANISOTROPIC_DIFFUSION)
      cpro_visma_v = (const cs_real_6_t *)CS_F_(vism)->val;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      const cs_lnum_t c_id = b_face_cells[f_id];
      const cs_real_t distbf  = b_dist[f_id];

      cs_real_3_t pimpv  = {0., 0., 0.};
      cs_real_3_t hextv  = {0., 0., 0.};
      cs_real_3_t qimpv  = {0., 0., 0.};
      cs_real_3_t cflv   = {0., 0., 0.};
      cs_real_6_t hintt  = {0., 0., 0., 0., 0., 0.};

      if (eqp_displ->idften & CS_ISOTROPIC_DIFFUSION) {

        hintt[0] = cpro_visma_s[c_id]/distbf;
        hintt[1] = cpro_visma_s[c_id]/distbf;
        hintt[2] = cpro_visma_s[c_id]/distbf;
        hintt[3] = 0.0;
        hintt[4] = 0.0;
        hintt[5] = 0.0;

      }
      else if (eqp_displ->idften & CS_ANISOTROPIC_DIFFUSION) {

        hintt[0] = cpro_visma_v[c_id][0]/distbf;
        hintt[1] = cpro_visma_v[c_id][1]/distbf;
        hintt[2] = cpro_visma_v[c_id][2]/distbf;
        hintt[3] = cpro_visma_v[c_id][3]/distbf;
        hintt[4] = cpro_visma_v[c_id][4]/distbf;
        hintt[5] = cpro_visma_v[c_id][5]/distbf;

      }

      /* Dirichlet Boundary Condition
         ---------------------------- */

      if (icodcl_displ[f_id] == 1) {

        for (cs_lnum_t k = 0; k < 3; k++)
          pimpv[k] = rcodcl1_displ[n_b_faces*k + f_id];

        for (cs_lnum_t k = 0; k < 3; k++)
          hextv[k] = rcodcl2_displ[n_b_faces*k + f_id];

        cs_boundary_conditions_set_dirichlet_vector_aniso(claale[f_id],
                                                          cfaale[f_id],
                                                          clbale[f_id],
                                                          cfbale[f_id],
                                                          pimpv,
                                                          hintt,
                                                          hextv);
      }

      /* Neumann Boundary Condition
         -------------------------- */

      else if (icodcl_displ[f_id] == 3) {

        /* Coupled solving of the velocity components */

        for (cs_lnum_t k = 0; k < 3; k++)
          qimpv[k] = rcodcl3_displ[n_b_faces*k + f_id];

        cs_boundary_conditions_set_neumann_vector_aniso(claale[f_id],
                                                        cfaale[f_id],
                                                        clbale[f_id],
                                                        cfbale[f_id],
                                                        qimpv,
                                                        hintt);
      }

      /* Convective Boundary Condition
         ----------------------------- */

      else if (icodcl_displ[f_id] == 2) {

        /* Coupled solving of the velocity components */

        for (cs_lnum_t k = 0; k < 3; k++)
          pimpv[k] = rcodcl1_displ[n_b_faces*k + f_id];

        for (cs_lnum_t k = 0; k < 3; k++)
          cflv[k] = rcodcl2_displ[n_b_faces*k + f_id];

        cs_boundary_conditions_set_convective_outlet_vector_aniso
          (claale[f_id], cfaale[f_id],
           clbale[f_id], cfbale[f_id],
           pimpv, cflv, hintt);
      }
    }
  }

  /*--------------------------------------------------------------------------
   * 15) Compute stresses at boundary (step 1 of 5)
   *--------------------------------------------------------------------------*/

  if (f_forbr != NULL && iterns == 1) {

    cs_real_3_t  *cofaf_vel = (cs_real_3_t  *)vel->bc_coeffs->af;
    cs_real_33_t *cofbf_vel = (cs_real_33_t *)vel->bc_coeffs->bf;

    /* Coupled solving of the velocity components */
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      const cs_real_t srfbnf = b_face_surf[f_id];

      /* The implicit term is added after having updated the velocity */

      forbr[f_id][0] =   srfbnf * (cofaf_vel[f_id][0]
                       + cofbf_vel[f_id][0][0]*velipb[f_id][0]
                       + cofbf_vel[f_id][1][0]*velipb[f_id][1]
                       + cofbf_vel[f_id][2][0]*velipb[f_id][2]);

      forbr[f_id][1] =   srfbnf * (cofaf_vel[f_id][1] +
                       + cofbf_vel[f_id][0][1]*velipb[f_id][0]
                       + cofbf_vel[f_id][1][1]*velipb[f_id][1]
                       + cofbf_vel[f_id][2][1]*velipb[f_id][2]);

      forbr[f_id][2] =   srfbnf * (cofaf_vel[f_id][2] +
                       + cofbf_vel[f_id][0][2]*velipb[f_id][0]
                       + cofbf_vel[f_id][1][2]*velipb[f_id][1]
                       + cofbf_vel[f_id][2][2]*velipb[f_id][2]);

    }
  }

  /* Free memory */
  BFT_FREE(velipb);

  /*--------------------------------------------------------------------------
   * 16) Update of boundary temperature when saved and not a variable.
   *--------------------------------------------------------------------------*/

  if (thermal_variable == CS_THERMAL_MODEL_ENTHALPY) {

    cs_field_t *f_b_temp = cs_field_by_name_try("boundary_temperature");

    if (f_b_temp != NULL) {

      cs_real_t *btemp_s = (cs_real_t *)f_b_temp->val;

      /* If we also have a boundary value field for the thermal
         scalar, copy its values first.

         If we do not have a boundary value field for the thermal scalar,
         boundary values for the thermal scalar were directly
         saved to the boundary temperature field, so no copy is needed. */

      int b_f_id = cs_field_get_key_int(f_th, kbfid);
      cs_real_t *bvar_s = NULL;

      if (b_f_id > -1)
        bvar_s = cs_field_by_id(b_f_id)->val;
      else {
        BFT_MALLOC(bvar_s, n_b_faces, cs_real_t);
        for (int f_id = 0; f_id < n_b_faces; f_id++) {
          bvar_s[f_id] = btemp_s[f_id];
        }
      }

      cs_ht_convert_h_to_t_faces(bvar_s, btemp_s);

      if (b_f_id < 0)
        BFT_FREE(bvar_s);

      /* In case of assigned temperature values, overwrite computed
         wall temperature with prescribed one to avoid issues due to
         enthalpy -> temperature conversion precision
         (T -> H -> T at the boundary does not preserve T) */
      for (cs_lnum_t ii = 0; ii < nbt2h; ii++) {
        const int f_id = lbt2h[ii];
        btemp_s[f_id] = vbt2h[f_id];
      }
    }
    BFT_FREE(lbt2h);
    BFT_FREE(vbt2h);
  }
}

/*---------------------------------------------------------------------------- */
/*!
 * \brief  Initialization of boundary condition arrays.
 */
/*---------------------------------------------------------------------------- */

void
cs_boundary_conditions_set_coeffs_init(void)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;
  const cs_lnum_t n_vertices  = mesh->n_vertices;
  const cs_real_3_t *vtx_coord = (const cs_real_3_t *)mesh->vtx_coord;

  int *bc_type = cs_f_boundary_conditions_get_bc_type();
  int *itrifb = cs_glob_bc_pm_info->itrifb;
  const int *izfppp = cs_glob_bc_pm_info->izfppp;
  cs_real_t *dt = CS_F_(dt)->val;

  const cs_lnum_t nt_cur  = cs_glob_time_step->nt_cur;
  const cs_lnum_t nt_prev = cs_glob_time_step->nt_prev;

  int *impale = cs_glob_ale_data->impale;
  int *ale_bc_type = cs_glob_ale_data->bc_type;

  cs_field_build_bc_codes_all();

  cs_boundary_conditions_reset();

  /* User calls
     ---------- */

  if (cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] >=  1)
    cs_f_ppprcl(bc_type, dt);

  /* NB. BC zones: we temporarily use specific physical model zones, even without
     the associated models.
     -> will be modified when restructuring BC zones. */

  cs_gui_boundary_conditions_processing(bc_type);

  /* User-defined function settings */

  cs_f_user_boundary_conditions_wrapper(itrifb, bc_type, izfppp, dt);

  cs_user_boundary_conditions(cs_glob_domain, bc_type);

  /* ALE BCs (mesh velocity and nodal displacement) */

  cs_field_t *f_displ = cs_field_by_name_try("mesh_displacement");
  if (f_displ != NULL) {

    cs_real_3_t *disale = (cs_real_3_t *)(f_displ->val);

    const cs_real_3_t *xyzno0
      = (const cs_real_3_t *)cs_field_by_name("vtx_coord0")->val;

    cs_array_lnum_fill_zero(n_vertices, impale);

    /* GUI and user-defined function-based definitions */

    cs_gui_mobile_mesh_boundary_conditions(ale_bc_type, impale, disale);

    cs_user_boundary_conditions_ale(cs_glob_domain,
                                    bc_type,
                                    ale_bc_type,
                                    impale);

    /* In case the user has modified disale whthout setting impale=1, we restore
       the initial displacement. */

    for (cs_lnum_t ii = 0; ii < n_vertices; ii++) {
      if (impale[ii] == 0) {
        disale[ii][0] = vtx_coord[ii][0] - xyzno0[ii][0];
        disale[ii][1] = vtx_coord[ii][1] - xyzno0[ii][1];
        disale[ii][2] = vtx_coord[ii][2] - xyzno0[ii][2];
      }
    }
  }

  /* For internal coupling, set itypfb to wall function by default
     if not set by the user. */

  cs_internal_coupling_bcs(bc_type);

  /* Treatment of types of bcs given by bc_type
     ----------------------------------------- */

  if (cs_glob_ale > CS_ALE_NONE)
    _boundary_condition_ale_type(mesh,
                                 cs_glob_mesh_quantities,
                                 true,
                                 dt,
                                 bc_type);

  if (cs_glob_rad_transfer_params->type != CS_RAD_TRANSFER_NONE)
    _boundary_condition_rt_type(mesh,
                                cs_glob_mesh_quantities,
                                true,
                                bc_type);

  if (cs_turbomachinery_get_model() != CS_TURBOMACHINERY_NONE)
    cs_f_mmtycl(bc_type);

  // Locate internal BC-based coupling
  if (cs_sat_coupling_n_couplings() > 0) {
    cs_f_cscloc();
    cs_f_cscfbr_init(bc_type);
  }

  if (   (   cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] >=  1
          && cs_glob_physical_model_flag[CS_GAS_MIX]             == -1
          && cs_glob_physical_model_flag[CS_JOULE_EFFECT]        == -1
          && cs_glob_physical_model_flag[CS_ELECTRIC_ARCS]       == -1)
      || (   cs_glob_physical_model_flag[CS_COMPRESSIBLE]       >=  0
          && cs_glob_physical_model_flag[CS_GAS_MIX]            >=  0)) {

    cs_f_pptycl(true,
                bc_type,
                izfppp,
                dt);
  }

  int *isostd;
  BFT_MALLOC(isostd, n_b_faces+1, int);

  cs_boundary_conditions_type(true,
                              bc_type,
                              itrifb,
                              isostd);

  BFT_FREE(isostd);

  /* Check the consistency of the BCs
     -------------------------------- */

  /* When called before time loop, some values are not yet available. */

  if (nt_cur > nt_prev) {
    cs_boundary_conditions_check(bc_type,
                                 ale_bc_type);
  }

  cs_field_free_bc_codes_all();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set convective oulet boundary condition for a scalar.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pimp   flux value to impose
 * \param[in]   cfl    local Courant number used to convect
 * \param[in]   hint   internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_convective_outlet_scalar(cs_real_t *a ,
                                                    cs_real_t *af,
                                                    cs_real_t *b,
                                                    cs_real_t *bf,
                                                    cs_real_t  pimp,
                                                    cs_real_t  cfl,
                                                    cs_real_t  hint)
{
  /* Gradient BCs */
  *b = cfl / (1.0 + cfl);
  *a = (1.0 - *b) * pimp;

  /* Flux BCs */
  *af = - hint * *a;
  *bf =   hint * (1.0 - *b);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized BC for an anisotropic symmetric vector for a given
 *         face.
 *
 * \param[out]  a       explicit BC coefficient for gradients
 * \param[out]  af      explicit BC coefficient for diffusive flux
 * \param[out]  b       implicit BC coefficient for gradients
 * \param[out]  bf      implicit BC coefficient for diffusive flux
 * \param[in]   pimpv   Dirichlet value to impose on the normal component
 * \param[in]   qimpv   flux value to impose on the tangential components
 * \param[in]   hint    internal exchange coefficient
 * \param[in]   normal  normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_generalized_sym_vector_aniso
  (cs_real_t        a[3],
   cs_real_t        af[3],
   cs_real_t        b[3][3],
   cs_real_t        bf[3][3],
   const cs_real_t  pimpv[3],
   const cs_real_t  qimpv[3],
   const cs_real_t  hint[6],
   const cs_real_t  normal[3])
{
  cs_real_t m[6] = {0., 0., 0., 0., 0., 0.};

  m[0] = hint[1]*hint[2] - hint[4]*hint[4];
  m[1] = hint[0]*hint[2] - hint[5]*hint[5];
  m[2] = hint[0]*hint[1] - hint[3]*hint[3];
  m[3] = hint[4]*hint[5] - hint[3]*hint[2];
  m[4] = hint[3]*hint[5] - hint[0]*hint[4];
  m[5] = hint[3]*hint[4] - hint[1]*hint[5];

  const cs_real_t invdet = 1.0/(hint[0]*m[0] + hint[3]*m[3] + hint[5]*m[5]);

  cs_real_t invh[6] = {0., 0., 0., 0., 0., 0.};
  invh[0] = m[0] * invdet;
  invh[1] = m[1] * invdet;
  invh[2] = m[2] * invdet;
  invh[3] = m[3] * invdet;
  invh[4] = m[4] * invdet;
  invh[5] = m[5] * invdet;

  cs_real_t qshint[3] = {0., 0., 0.};
  cs_real_t hintpv[3] = {0., 0., 0.};
  cs_real_t hintnm[3] = {0., 0., 0.};

  cs_math_sym_33_3_product(invh, qimpv,  qshint);
  cs_math_sym_33_3_product(hint, pimpv,  hintpv);
  cs_math_sym_33_3_product(hint, normal, hintnm);

  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    a[isou] = - qshint[isou];
    /* "[1 -n(x)n] Qimp / hint" is divided into two */
    for (int jsou = 0; jsou < 3; jsou++) {

      a[isou] = a[isou] + normal[isou]*normal[jsou]
        * (pimpv[jsou] + qshint[jsou]);

      if (jsou == isou)
        b[isou][jsou] = 1.0 - normal[isou]*normal[jsou];
      else
        b[isou][jsou] = - normal[isou]*normal[jsou];
    }

    /* Flux BCs */
    af[isou] = qimpv[isou];
    /* "[1 -n(x)n] Qimp" is divided into two */
    for (int jsou = 0; jsou < 3; jsou++){
      af[isou] = af[isou] - normal[isou]*normal[jsou]
                  * (hintpv[jsou] + qimpv[jsou]);

      bf[isou][jsou] = hintnm[isou] * normal[jsou];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for an anisotropic vector for a given
 *         face.
 *
 * \param[out]  a        explicit BC coefficient for gradients
 * \param[out]  af       explicit BC coefficient for diffusive flux
 * \param[out]  b        implicit BC coefficient for gradients
 * \param[out]  bf       implicit BC coefficient for diffusive flux
 * \param[in]   pimpv    Dirichlet value to impose on the tangential components
 * \param[in]   qimpv    flux value to impose on the normal component
 * \param[in]   hint     internal exchange coefficient
 * \param[in]   normal   normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_generalized_dirichlet_vector_aniso
  (cs_real_t        a[3],
   cs_real_t        af[3],
   cs_real_t        b[3][3],
   cs_real_t        bf[3][3],
   const cs_real_t  pimpv[3],
   const cs_real_t  qimpv[3],
   const cs_real_t  hint[6],
   const cs_real_t  normal[3])
{
  cs_real_t m[6] = {0., 0., 0., 0., 0., 0.};
  m[0] = hint[1]*hint[2] - hint[4]*hint[4];
  m[1] = hint[0]*hint[2] - hint[5]*hint[5];
  m[2] = hint[0]*hint[1] - hint[3]*hint[3];
  m[3] = hint[4]*hint[5] - hint[3]*hint[2];
  m[4] = hint[3]*hint[5] - hint[0]*hint[4];
  m[5] = hint[3]*hint[4] - hint[1]*hint[5];

  const cs_real_t invdet = 1.0/(hint[0]*m[0] + hint[3]*m[3] + hint[5]*m[5]);

  cs_real_t invh[6] = {0., 0., 0., 0., 0., 0.};
  invh[0] = m[0] * invdet;
  invh[1] = m[1] * invdet;
  invh[2] = m[2] * invdet;
  invh[3] = m[3] * invdet;
  invh[4] = m[4] * invdet;
  invh[5] = m[5] * invdet;

  cs_real_t qshint[3] = {0., 0., 0.};
  cs_real_t hintpv[3] = {0., 0., 0.};
  cs_real_t hintnm[3] = {0., 0., 0.};

  cs_math_sym_33_3_product(invh, qimpv,  qshint);
  cs_math_sym_33_3_product(hint, pimpv,  hintpv);
  cs_math_sym_33_3_product(hint, normal, hintnm);

  for (int isou = 0; isou < 3; isou ++) {

    /* Gradient BCs */
    /* "[1 -n(x)n] Pimp" is divided into two */
    a[isou] = pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++) {

      a[isou] = a[isou] - normal[isou] * normal[jsou]
                  * (pimpv[jsou] + qshint[jsou]);

      b[isou][jsou] = normal[isou] * normal[jsou];
    }

    /* Flux BCs */
    /* "[1 -n(x)n] Pimp" is divided into two */
    af[isou] = -hintpv[isou];
    for (int jsou = 0; jsou < 3; jsou++) {

      af[isou] = af[isou] + normal[isou]*normal[jsou]
        *(qimpv[jsou]+hintpv[jsou]);

      if (jsou == isou)
        bf[isou][jsou] = hint[isou]-hintnm[isou]*normal[jsou];
      else
        bf[isou][jsou] = -hintnm[isou]*normal[jsou];
    }
  }
}

/*---------------------------------------------------------------------------- */

END_C_DECLS
