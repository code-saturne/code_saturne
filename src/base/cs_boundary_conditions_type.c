/*============================================================================
 * Handle boundary condition type codes (bc_type).
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

#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_array.h"
#include "cs_atmo.h"
#include "cs_boundary_conditions.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_intprf.h"
#include "cs_lagr.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_sat_coupling.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_velocity_pressure.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions_type.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *============================================================================*/

static bool _initialized = false;

/* Closest free standard outlet face to xyzp0 (icodcl not modified)
   (or closest free inlet) */
static cs_lnum_t _ifrslb = -1;

/* Max of _ifrslb on all ranks, standard outlet face presence indicator */
static cs_lnum_t _itbslb = -1;

/* Rank associated with _itbslb */
static int _irangd = 0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return count of faces with a given boundary condition type.
 *
 * parameters:
 *   bc_type      <-- boundary type
 *   bc_type_idx  <-- boundary type index
 *
 * return:
 *   updated bit mask
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_n_bc_faces(int        bc_type,
            const int  bc_type_idx[])
{
  cs_lnum_t n = bc_type_idx[bc_type] - bc_type_idx[bc_type-1];

  return n;
}

/*============================================================================
 * Public function definitions.
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Handle boundary condition type code (\ref bc_type).
 *
 * \param[in]      init      partial treatment (before time loop) if true
 * \param[in,out]  bc_type   type per boundary face
 * \param[out]     itrifb    indirection for faces ordering
 * \param[out]     isostd    standard output indicator + reference face number
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_type(bool  init,
                            int   bc_type[],
                            int   itrifb[],
                            int   isostd[])
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)mesh->b_face_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_real_3_t *cell_cen  = (const cs_real_3_t *)fvq->cell_cen;
  const cs_real_3_t *b_face_cog = (const cs_real_3_t *)fvq->b_face_cog;
  cs_velocity_pressure_param_t *vp_param
    = cs_get_glob_velocity_pressure_param();

  const cs_lnum_t nt_cur = cs_glob_time_step->nt_cur;
  const cs_lnum_t nt_max = cs_glob_time_step->nt_max;
  const cs_lnum_t nt_prev = cs_glob_time_step->nt_prev;

  cs_fluid_properties_t *fluid_props = cs_get_glob_fluid_properties();
  const cs_real_t ro0 = fluid_props->ro0;
  const cs_real_t p0 = fluid_props->p0;
  const cs_real_t *gxyz = cs_get_glob_physical_constants()->gravity;
  cs_real_t *xyzp0 = fluid_props->xyzp0;

  const int kturt  = cs_field_key_id("turbulent_flux_model");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const int keysca = cs_field_key_id("scalar_id");
  const int meteo_profile = cs_glob_atmo_option->meteo_profile;
  const int n_fields = cs_field_n_fields();

  int *icodcl_p = NULL, *icodcl_vel = NULL;
  cs_real_t *rcodcl1_p = NULL, *rcodcl2_p = NULL, *rcodcl3_p = NULL;
  cs_real_t *rcodcl1_vel = NULL, *rcodcl2_vel = NULL, *rcodcl3_vel = NULL;

  /* Type ids and names (general case) */

  const int type_id[11] = {CS_INLET, CS_SMOOTHWALL, CS_ROUGHWALL,
                           CS_SYMMETRY, CS_OUTLET, CS_FREE_INLET,
                           CS_CONVECTIVE_INLET,
                           CS_COUPLED, CS_COUPLED_FD,
                           CS_FREE_SURFACE, CS_INDEF};

  const char *info[11] = {"Inlet", "Smooth wall", "Rough wall",
                          "Symmetry", "Free outlet", "Free inlet",
                          "Convective inlet",
                          "cs/cs coupling", "cs/cs coupling",
                          "Free surface", "Undefined"};

  /* Type ids and names (incompressible case) */

  const int type_id_c[9] = {CS_EQHCF, CS_EPHCF, CS_ESICF,
                            CS_SOPCF, CS_SSPCF,
                            CS_SMOOTHWALL, CS_ROUGHWALL,
                            CS_SYMMETRY, CS_INDEF};

  const char *info_c[9] = {"Sub. enth. inlet", "Ptot, Htot",
                           "Imp inlet/outlet", "Subsonic outlet",
                           "Supersonic outlet",
                           "Smooth wall", "Rough wall",
                           "Symmetry", "Undefined"};

  /* Initialization
     ============== */

  if (CS_F_(p) != NULL) {
    icodcl_p = CS_F_(p)->bc_coeffs->icodcl;
    rcodcl1_p = CS_F_(p)->bc_coeffs->rcodcl1;
    rcodcl2_p = CS_F_(p)->bc_coeffs->rcodcl2;
    rcodcl3_p = CS_F_(p)->bc_coeffs->rcodcl3;
  }

  if (CS_F_(vel) != NULL) {
    icodcl_vel = CS_F_(vel)->bc_coeffs->icodcl;
    rcodcl1_vel = CS_F_(vel)->bc_coeffs->rcodcl1;
    rcodcl2_vel = CS_F_(vel)->bc_coeffs->rcodcl2;
    rcodcl3_vel = CS_F_(vel)->bc_coeffs->rcodcl3;
  }

  /* Allocate temporary arrays */
  cs_real_t *pripb = NULL;
  BFT_MALLOC(pripb, n_b_faces, cs_real_t);

  cs_lnum_t *bc_type_idx = NULL;
  BFT_MALLOC(bc_type_idx, CS_MAX_BC_TYPE+1, cs_lnum_t);

  /* Initialize variables to avoid compiler warnings */
  cs_array_real_fill_zero(n_b_faces, pripb);
  cs_real_t pref = 0.;

  const cs_real_t *cvara_pr = (const cs_real_t *)CS_F_(p)->val_pre;
  const cs_equation_param_t *eqp_vel
    = cs_field_get_equation_param_const(CS_F_(vel));

  cs_real_t *b_head_loss = cs_boundary_conditions_get_b_head_loss(false);

  /* Check consistency of types given in cs_user_boundary_conditions
     =============================================================== */

  {
    int iok = 0;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      int type = bc_type[f_id];
      if (type <= 0 || type > CS_MAX_BC_TYPE) {
        bc_type[f_id] = 0;
        iok = 1;
      }
    }
    cs_parall_max(1, CS_INT_TYPE, &iok);

    if (iok != 0)
      cs_boundary_conditions_error(bc_type, NULL);
  }

  /* Sort boundary faces
     =================== */

  {
    /* Count faces of each type, then transform count to index */

    for (int i = 0; i < CS_MAX_BC_TYPE+1; i++)
      bc_type_idx[i] = 0;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      const int bc_type_id = bc_type[f_id] - 1;  /* 1-based */
      bc_type_idx[bc_type_id+1] += 1;
    }

    for (int i = 0; i < CS_MAX_BC_TYPE; i++)
      bc_type_idx[i+1] += bc_type_idx[i];

    /* Build list of face ids sorted by bc type: itrifb */

    cs_lnum_t *bc_type_count = NULL;
    BFT_MALLOC(bc_type_count, CS_MAX_BC_TYPE, cs_lnum_t);
    for (int j = 0; j < CS_MAX_BC_TYPE; j++)
      bc_type_count[j] = 0;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      const int bc_type_id = bc_type[f_id] - 1;
      const int k = bc_type_idx[bc_type_id] + bc_type_count[bc_type_id];
      itrifb[k] = f_id;
      bc_type_count[bc_type_id] += 1;
    }

    BFT_FREE(bc_type_count);

    /* Number of boundary faces classified by type */

#if defined(DEBUG) && !defined(NDEBUG)

    cs_gnum_t isum = 0;
    for (int ii = 0; ii < CS_MAX_BC_TYPE; ii++)
      isum += bc_type_idx[ii+1] - bc_type_idx[ii];

    cs_parall_sum(1, CS_GNUM_TYPE, &isum);

    if (isum != cs_glob_mesh->n_g_b_faces)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error ordering of boundary faces.\n\n"
                  "The number of faces classified by type = %llu does not \n"
                  "correspond to the global number of boundary faces = %llu."),
                __func__, (unsigned long long)isum,
                (unsigned long long)cs_glob_mesh->n_g_b_faces);

#endif
  }

  /* Write boundary type with corresponding code and faces number */

  if (_initialized == false || eqp_vel->iwarni > 1) {

    _initialized = true;

    cs_log_printf
      (CS_LOG_DEFAULT,
       _("   ** INFORMATION ON BOUNDARY FACES TYPE\n"
         "      ----------------------------------\n\n"));

    cs_log_separator(CS_LOG_DEFAULT);
    cs_log_printf(CS_LOG_DEFAULT,
                  _("Boundary type          Code    Nb faces\n"));
    cs_log_separator(CS_LOG_DEFAULT);

    /* Global counts par BC type;
       Note: we do not need so many BC user types, so we could make
       CS_MAX_BC_TYPE much smaller, reducing the size of
       the global counter; using zone-based definitions
       and flags from cs_boundary.h, we could even avoid this count. */

    cs_gnum_t inb[CS_MAX_BC_TYPE];
    char is_user_type[CS_MAX_BC_TYPE];

    for (int ii = 0; ii < CS_MAX_BC_TYPE; ii++) {
      inb[ii] = bc_type_idx[ii+1] - bc_type_idx[ii];
      is_user_type[ii] = 1;
    }

    cs_parall_counter(inb, CS_MAX_BC_TYPE);

    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {

      /* Common type ids.
         We will log counts based on these type ids first, in that
         order, so as to make comparisons with prior code versions
         easier. */

      int log_always[11] = {1, 1, 1,
                            1, 1, 1,
                            1,
                            0, 0,
                            1, 1};

      if (cs_sat_coupling_n_couplings() >= 1) {
        if (cs_glob_sat_coupling_face_interpolation_type == 0)
          log_always[7] = 1;
        else
          log_always[8] = 1;
      }

      for (int jj = 0; jj < 11; jj++) {
        const int ii = type_id[jj] - 1;
        if (log_always[jj] || inb[ii] > 0)
          cs_log_printf(CS_LOG_DEFAULT, _("%-17s  %8d%12llu\n"),
                        info[jj], ii+1, (unsigned long long)inb[ii]);
        is_user_type[ii] = 0;
      }

      /* Presence of free inlet face(s) */

      if (inb[CS_FREE_INLET - 1] > 0) {
        vp_param->iifren = 1;
        b_head_loss = cs_boundary_conditions_get_b_head_loss(true);
      }
      else
        vp_param->iifren = 0;

    }
    else {

      for (int jj = 0; jj < 9; jj++) {
        const int ii = type_id_c[jj] - 1;
        cs_log_printf(CS_LOG_DEFAULT, _("%-17s  %8d%12llu\n"),
                      info_c[jj], ii+1, (unsigned long long)inb[ii]);
        is_user_type[ii] = 0;
      }

    } /* End of test on compressible/incompressible case */

    /* User types */

    for (int ii = 0; ii < CS_MAX_BC_TYPE; ii++) {
      if (is_user_type[ii] != 1)
        continue;
      cs_gnum_t inb_user = bc_type_idx[ii+1] - bc_type_idx[ii];
      if (inb_user > 0)
        cs_log_printf
          (CS_LOG_DEFAULT,_("%-17s %8d %12llu\n"),
           "User type", ii+1,  (unsigned long long)inb_user);
    }

    cs_log_separator(CS_LOG_DEFAULT);
  }

  /* rcodcl1 arrrays have been initialized as cs_math_infinite_r so as to
     check what the user has modified. Those not modified are reset to zero
     here. CS_OUTLET, CS_FREE_INLET and CS_INLET are handled later.
     ====================================================================== */

  {
    for (int field_id = 0; field_id < n_fields; field_id++) {
      cs_field_t *f = cs_field_by_id(field_id);

      if (!(f->type & CS_FIELD_VARIABLE))
        continue;
      if (f->type & CS_FIELD_CDO)
        continue;

      cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        if (   (bc_type[f_id] != CS_OUTLET)
            && (bc_type[f_id] != CS_FREE_INLET)
            && (bc_type[f_id] != CS_FREE_SURFACE)
            && (bc_type[f_id] != CS_CONVECTIVE_INLET)
            && (bc_type[f_id] != CS_INLET)
            && (rcodcl1[f_id] > 0.5*cs_math_infinite_r)) {
          for (cs_lnum_t k = 0; k < f->dim; k++)
            rcodcl1[k*n_b_faces+f_id] = 0.;
        }
      }
    }
  }

  /* Compute pressure at boundary (in pripb(*))
     (if we need it, that is if there are outlet boudary faces)
     ========================================================== */

  /* ifrslb = closest free standard outlet face to xyzp0 (icodcl not modified)
     (or closest free inlet)
     itbslb = max of ifrslb on all ranks,
     standard outlet face presence indicator */

  /* Even when the user has not chosen xyzp0 (and it is thus at the
     origin), we choose the face whose center is closest to it, so
     as to be mesh numbering (and partitioning) independent. */

  if (nt_cur == nt_prev+1) {

    cs_real_t d0min = cs_math_infinite_r;

    const int o_type_id[2] = {CS_OUTLET-1,
                              CS_FREE_INLET-1};

    for (int jj = 0; jj < 2; jj++) {
      cs_lnum_t s_id = bc_type_idx[o_type_id[jj]];
      cs_lnum_t e_id = bc_type_idx[o_type_id[jj] + 1];
      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t f_id = itrifb[ii];
        if (icodcl_p[f_id] == 0) {
          const cs_real_t d0 = cs_math_3_square_distance(xyzp0,
                                                         b_face_cog[f_id]);
          if (d0 < d0min) {
            _ifrslb = f_id;
            d0min = d0;
          }
        }
      }
    }

    /* If we have free outlet faces, irangd and itbslb will
       contain respectively the rank having the boundary face whose
       center is closest to xyzp0, and the local number of that face
       on that rank (also equal to ifrslb on that rank).
       If we do not have free outlet faces, than itbslb = 0
       (as it was initialized that way on all ranks). */

    _itbslb = _ifrslb;
    _irangd = cs_glob_rank_id;

    cs_parall_min_id_rank_r(&_itbslb, &_irangd, d0min);
  }

  if (vp_param->iphydr == 1 || vp_param->iifren == 1) {

    /* If the hydrostatic pressure is taken into account,
       we fill the isostd array.
       0 -> not a standard outlet face (i.e. not outlet or outlet with
            modified pressure BC)
       1 -> free outlet face with automatic pressure BC.
       the reference face number is stored in isostd[]
       which is first initialized to -1 (i.e. no standard output face) */

    isostd[n_b_faces] = -1;
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      isostd[f_id] = 0;
      if (   (bc_type[f_id] == CS_OUTLET || bc_type[f_id] == CS_FREE_INLET)
          && (icodcl_p[f_id] == 0))
        isostd[f_id] = 1;
    }
  }

  /* Reference pressure (unique, even if there are multiple outlets)
     In case we account for the hydrostatic pressure, we search for the
     reference face.

     Determine a unique P I' pressure in parallel
     if there are free outlet faces, we have determined that the rank
     with the outlet face closest to xyzp0 is irangd.

     We also retrieve the coordinates of the reference point, so as to
     calculate pref later on. */

  cs_real_t xyzref[3] = {0., 0., 0.};

  if (_itbslb > -1) {

    /* If irangd is the local rank, we assign PI' to xyzref(4)
       (this is always the case in serial mode) */

    if (cs_glob_rank_id == _irangd) {
      xyzref[0] = b_face_cog[_ifrslb][0];
      xyzref[1] = b_face_cog[_ifrslb][1];
      xyzref[2] = b_face_cog[_ifrslb][2];
      if (vp_param->iphydr == 1 || vp_param->iifren == 1)
        isostd[n_b_faces] = _ifrslb;
    }

    /* Broadcast PI' and pressure reference
       from irangd to all other ranks. */

    cs_parall_bcast(_irangd, 3, CS_REAL_TYPE, xyzref);

    /* If the user has not specified anything, we set ixyzp0 to 2 so as
       to update the reference point. */

    if (fluid_props->ixyzp0 == -1)
      fluid_props->ixyzp0 = 2;

  }
  else if (fluid_props->ixyzp0 < 0 && nt_cur == nt_prev + 1) {

    /* If there are no outlet faces, we search for possible Dirichlets
       specified by the user so as to locate the reference point.
       As before, we chose the face closest to xyzp0 so as to
       be mesh numbering (and partitioning) independent. */

    cs_lnum_t d0min = cs_math_infinite_r;

    cs_lnum_t ifadir = -1;
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      if (abs(icodcl_p[f_id]) == 1) {
        const cs_real_t d0 = cs_math_3_square_distance(xyzp0,
                                                       b_face_cog[f_id]);
        if (d0 < d0min) {
          ifadir = f_id;
          d0min = d0;
        }
      }
    }

    _irangd = cs_glob_rank_id;
    cs_parall_min_id_rank_r(&ifadir, &_irangd, d0min);

    if (ifadir > -1)
      /* We set ixyzp0 to 2 to update the reference point */
      fluid_props->ixyzp0 = 2;

    if (cs_glob_rank_id == _irangd) {
      xyzref[0] = b_face_cog[ifadir][0];
      xyzref[1] = b_face_cog[ifadir][1];
      xyzref[2] = b_face_cog[ifadir][2];
    }

    /* Broadcast xyzref from irangd to all other ranks. */
    cs_parall_bcast(_irangd, 3, CS_REAL_TYPE, xyzref);
  }

  /* If the reference point has not been specified by the user,
     we change it and shift coefu if there are outputs.
     Total pressure is also shifted (a priori useless exept if
     the user uses it in cs_user_source_terms for example) */

  if (fluid_props->ixyzp0 == 2) {

    fluid_props->ixyzp0 = 1;

    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {

      cs_real_t *cpro_prtot = cs_field_by_name_try("total_pressure")->val;
      if (cpro_prtot != NULL) {
        cs_real_t prtot_shift = -ro0 * cs_math_3_distance_dot_product(xyzp0,
                                                                      xyzref,
                                                                      gxyz);
        for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
          cpro_prtot[c_id] += prtot_shift;
      }

    }

    xyzp0[0] = xyzref[0];
    xyzp0[1] = xyzref[1];
    xyzp0[2] = xyzref[2];

    if (_itbslb > -1)
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n"
           "Boundary faces with free inlet/outlet detected.\n"
           "Update of reference point for total pressure.\n"
           "  xyzp0 = (%14.5e,%14.5e,%14.5e)"),
         xyzp0[0], xyzp0[1], xyzp0[2]);

    else
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n"
           "Boundary faces with pressure Dirichlet condition detected.\n"
           "Update of reference point for total pressure\n"
           "  xyzp0 = (%14.5e,%14.5e,%14.5e)"),
         xyzp0[0], xyzp0[1], xyzp0[2]);
  }

  else if (fluid_props->ixyzp0 == -1) {

    /* There are no outputs and no Dirichlet, and the user has not
       specified anything. We set IXYZP0 to 0 so as not to touch it again,
       in contrast to the =1 case, which will require subsequent writing. */
    fluid_props->ixyzp0 = 0;

  }

  /* No need to compute pressure gradient for frozen field computations */

  if (   _itbslb > -1
      && cs_glob_lagr_time_scheme->iilagr != CS_LAGR_FROZEN_CONTINUOUS_PHASE) {

    cs_real_3_t *frcxt = NULL;
    {
      cs_field_t *f_vf = cs_field_by_name_try("volume_forces");
      if (f_vf != NULL)
        frcxt = (cs_real_3_t *)f_vf->val;
    }

    /* Allocate a work array for the gradient calculation */
    cs_real_3_t *grad;
    BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

    cs_gradient_porosity_balance(1);

    cs_field_gradient_potential(CS_F_(p),
                                true,
                                1,
                                vp_param->iphydr,
                                frcxt,
                                grad);

    /* Put in pripb the value at F of the
       total pressure, computed from P* shifted from the ro0*g.(X-Xref) */

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      const cs_lnum_t c_id = b_face_cells[f_id];

      /* IF: Direction of projection of the pressure gradient */
      cs_real_t proj_dir[3];
      proj_dir[0] = b_face_cog[f_id][0] - cell_cen[c_id][0];
      proj_dir[1] = b_face_cog[f_id][1] - cell_cen[c_id][1];
      proj_dir[2] = b_face_cog[f_id][2] - cell_cen[c_id][2];

      if (meteo_profile > 0)  {

        cs_real_t vel_dir[3];
        if (meteo_profile == 1) {
          const int met_1d_nlevels_t = cs_glob_atmo_option->met_1d_nlevels_t;
          const int met_1d_ntimes = cs_glob_atmo_option->met_1d_ntimes;

          const cs_real_t xuent = cs_intprf(met_1d_nlevels_t,
                                            met_1d_ntimes,
                                            cs_glob_atmo_option->z_temp_met,
                                            cs_glob_atmo_option->time_met,
                                            cs_glob_atmo_option->u_met,
                                            cell_cen[c_id][2],
                                            cs_glob_time_step->t_cur);

          const cs_real_t xvent = cs_intprf(met_1d_nlevels_t,
                                            met_1d_ntimes,
                                            cs_glob_atmo_option->z_temp_met,
                                            cs_glob_atmo_option->time_met,
                                            cs_glob_atmo_option->v_met,
                                            cell_cen[c_id][2],
                                            cs_glob_time_step->t_cur);

          const cs_real_t xwent = 0.; /* 2D Profile */

          /* Meteo Velocity direction */
          vel_dir[0] = xuent;
          vel_dir[1] = xvent;
          vel_dir[2] = xwent;
        }
        else if (meteo_profile == 2) {
          /* Meteo Velocity direction */
          cs_real_3_t *cpro_met_vel
            = (cs_real_3_t *)cs_field_by_name("meteo_velocity")->val;
          vel_dir[0] = cpro_met_vel[c_id][0];
          vel_dir[1] = cpro_met_vel[c_id][1];
          vel_dir[2] = cpro_met_vel[c_id][2];
        }

        /* Velocity direction normalized */
        cs_math_3_normalize(vel_dir, vel_dir);

        /* (IF.n) n */
        const cs_real_t vs = cs_math_3_distance_dot_product(cell_cen[c_id],
                                                            b_face_cog[f_id],
                                                            vel_dir);

        proj_dir[0] = vs*vel_dir[0];
        proj_dir[1] = vs*vel_dir[1];
        proj_dir[2] = vs*vel_dir[2];
      }

      pripb[f_id] = cvara_pr[c_id] + cs_math_3_dot_product(proj_dir,
                                                           grad[c_id]);
    }

    BFT_FREE(grad);

    if (cs_glob_rank_id == _irangd)
      pref = pripb[_ifrslb];

    cs_parall_bcast(_irangd, 1, CS_REAL_TYPE, &pref);
  }

  /* Convert to rcodcl and icodcl
     (if this has not already been set by the user)

     First, process variables for which a specific treatement is done
     (pressure, velocity, ...)
     ================================================================ */

  /* Translate total pressure P_tot given by the user into solved pressure P
     P = P_tot - p0 - rho0 ( g . (z-z0))

     If icodcl = -1, then the BC Dirichlet value is given in solved pressure P,
     no need of transformation from P_tot to P */

  if (CS_F_(p) != NULL) {

    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {
      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        if (icodcl_p[f_id] == -1)
          icodcl_p[f_id] = 1;
        else if (icodcl_p[f_id] != 0)
          rcodcl1_p[f_id] += - ro0 * cs_math_3_distance_dot_product
                                       (xyzp0,
                                        b_face_cog[f_id],
                                        gxyz)
                             - p0;
      }
    }

  }

  /* Inlet + Convective Inlet
     ------------------------ */

  /* (The pressure has a Neumann processing, the other Dirichlet
     will be processed later) */

  if (CS_F_(p) != NULL) {

    const int i_type_id[2] = {CS_INLET-1,
                              CS_CONVECTIVE_INLET-1};

    for (int jj = 0; jj < 2; jj++) {

      cs_lnum_t s_id = bc_type_idx[i_type_id[jj]];
      cs_lnum_t e_id = bc_type_idx[i_type_id[jj] + 1];

      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t f_id = itrifb[ii];
        if (icodcl_p[f_id] == 0) {
          icodcl_p[f_id] = 3;
          rcodcl1_p[f_id] = 0.;
          rcodcl2_p[f_id] = cs_math_infinite_r;
          rcodcl3_p[f_id] = 0.;
        }
      }

    }

  }

  /* OUTLET (free inlet-outlet)
     -------------------------- */

  /* Pressure has a Dirichlet condition, velocity code 9
     (the rest Neumann, or Dirichlet if user data,
     will be handled later) */

  /* Free inlet/outlet */

  {
    const cs_lnum_t s_id = bc_type_idx[CS_OUTLET-1];
    const cs_lnum_t e_id = bc_type_idx[CS_OUTLET];

    if (CS_F_(p) != NULL) {

      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t f_id = itrifb[ii];
        if (icodcl_p[f_id] == 0) {
          icodcl_p[f_id] = 1;
          rcodcl1_p[f_id] = pripb[f_id] - pref;
          rcodcl2_p[f_id] = cs_math_infinite_r;
          rcodcl3_p[f_id] = 0.;
        }
      }

    }

    if (CS_F_(vel) != NULL) {

      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t f_id = itrifb[ii];
        if (icodcl_vel[f_id] == 0) {
          icodcl_vel[f_id] = 9;

          for (cs_lnum_t k = 0; k < 3; k++) {
            rcodcl1_vel[n_b_faces*k + f_id] = 0.;
            rcodcl2_vel[n_b_faces*k + f_id] = cs_math_infinite_r;
            rcodcl3_vel[n_b_faces*k + f_id] = 0.;
          }
        }
      }

    }
  }

  /* Free inlet (Bernoulli relation), std free outlet
     (a specific treatment is performed on the pressure increment) */

  {
    const cs_lnum_t s_id = bc_type_idx[CS_FREE_INLET-1];
    const cs_lnum_t e_id = bc_type_idx[CS_FREE_INLET];

    /* Pressure */
    if (CS_F_(p) != NULL) {

      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t f_id = itrifb[ii];

        if (icodcl_p[f_id] == 0) {

          /* If the user has given a value of boundary head loss */
          if (rcodcl2_p[f_id] <= 0.5*cs_math_infinite_r)
            b_head_loss[f_id] = rcodcl2_p[f_id];
          else
            b_head_loss[f_id] = 0.;

          /* Std outlet */
          icodcl_p[f_id] = 1;
          rcodcl1_p[f_id] = pripb[f_id] - pref;
          rcodcl2_p[f_id] = cs_math_infinite_r;
          rcodcl3_p[f_id] = 0.;
        }
      }
    }

    /* Velocity */
    if (CS_F_(vel) != NULL) {

      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t f_id = itrifb[ii];
        /* Homogeneous Neumann */
        if (icodcl_vel[f_id] == 0) {
          icodcl_vel[f_id] = 3;

          for (cs_lnum_t k = 0; k < 3; k++) {
            rcodcl1_vel[n_b_faces*k + f_id] = 0.;
            rcodcl2_vel[n_b_faces*k + f_id] = cs_math_infinite_r;
            rcodcl3_vel[n_b_faces*k + f_id] = 0.;
          }
        }
      }

    }

  }

  /* Free surface: Dirichlet on the pressure */

  {
    const cs_lnum_t s_id = bc_type_idx[CS_FREE_SURFACE-1];
    const cs_lnum_t e_id = bc_type_idx[CS_FREE_SURFACE];

    if (CS_F_(p) != NULL) {

      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t f_id = itrifb[ii];
        if (icodcl_p[f_id] == 0) {
          icodcl_p[f_id] = 1;
          rcodcl1_p[f_id]
            = - ro0 * cs_math_3_distance_dot_product(xyzp0,
                                                     b_face_cog[f_id],
                                                     gxyz);
          rcodcl2_p[f_id] = cs_math_infinite_r;
          rcodcl3_p[f_id] = 0.;
        }
      }

    }

    if (CS_F_(vel) != NULL) {

      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t f_id = itrifb[ii];
        /* Homogeneous Neumann */
        if (icodcl_vel[f_id] == 0) {
          icodcl_vel[f_id] = 3;

          for (cs_lnum_t k = 0; k < 3; k++) {
            rcodcl1_vel[n_b_faces*k + f_id] = 0.;
            rcodcl2_vel[n_b_faces*k + f_id] = cs_math_infinite_r;
            rcodcl3_vel[n_b_faces*k + f_id] = 0.;
          }
        }
      }
    }

  }

  /* Free memory */
  BFT_FREE(pripb);

  /* Symmetry
     -------- */

  /* Vectors and tensors have a special treatment.
     Scalars have an homogeneous Neumann. */

  if (_n_bc_faces(CS_SYMMETRY, bc_type_idx) > 0) {

    const cs_lnum_t s_id = bc_type_idx[CS_SYMMETRY-1];
    const cs_lnum_t e_id = bc_type_idx[CS_SYMMETRY];

    for (int field_id = 0; field_id < n_fields; field_id++) {

      cs_field_t *f = cs_field_by_id(field_id);

      if (!(f->type & CS_FIELD_VARIABLE))
        continue;
      if (f->type & CS_FIELD_CDO)
        continue;

      bool is_uncoupled_rij = false;
      if (   cs_glob_turb_rans_model->irijco == 0
          && f == CS_F_(rij))
        is_uncoupled_rij = true;

      int *icodcl = f->bc_coeffs->icodcl;
      cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1;
      cs_real_t *rcodcl2 = f->bc_coeffs->rcodcl2;
      cs_real_t *rcodcl3 = f->bc_coeffs->rcodcl3;

      /* Loop over faces */
      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t f_id = itrifb[ii];

        /* Special treatment for uncoupled version of Rij models */
        if (icodcl[f_id] == 0 && is_uncoupled_rij) {
          icodcl[f_id] = 4;
          for (cs_lnum_t k = 0; k < 6; k++) {
            rcodcl1[k*n_b_faces+f_id] = 0.;
            rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
            rcodcl3[k*n_b_faces+f_id] = 0.;
          }
        }

        /* Homogeneous Neumann on scalars */
        if (f->dim == 1 && icodcl[f_id] == 0) {
          icodcl[f_id] = 3;
          rcodcl1[f_id] = 0.;
          rcodcl2[f_id] = cs_math_infinite_r;
          rcodcl3[f_id] = 0.;
        }

        /* Symmetry BC if nothing is set by the user on vector and tensors */
        else if (icodcl[f_id] == 0) {
          icodcl[f_id] = 4;
          for (cs_lnum_t k = 0; k < f->dim; k++) {
            rcodcl1[k*n_b_faces+f_id] = 0.;
            rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
            rcodcl3[k*n_b_faces+f_id] = 0.;
          }
        }
      }

    } /* end of loop on fields */

  }

  /* Smooth and rough walls
     ---------------------- */

  /* Velocity and turbulent quantities have icodcl = 5
     Turbulent fluxes of scalars has 0 Dirichlet if scalars
     have a Dirichlet, otherwise treated in clptur.
     Other quantities are treated afterwards (Homogeneous Neumann) */

  for (int wall_bc_type = CS_SMOOTHWALL;
       wall_bc_type <= CS_ROUGHWALL;
       wall_bc_type++) {

    if (_n_bc_faces(wall_bc_type, bc_type_idx) <= 0)
      continue;

    const int wall_bc_code = (wall_bc_type == CS_ROUGHWALL) ? 6 : 5;

    const cs_lnum_t s_id = bc_type_idx[wall_bc_type-1];
    const cs_lnum_t e_id = bc_type_idx[wall_bc_type];

    for (int field_id = 0; field_id < n_fields; field_id++) {

      cs_field_t *f = cs_field_by_id(field_id);

      if (!(f->type & CS_FIELD_VARIABLE))
        continue;
      if (f->type & CS_FIELD_CDO)
        continue;

      int *icodcl = f->bc_coeffs->icodcl;
      cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1;
      cs_real_t *rcodcl2 = f->bc_coeffs->rcodcl2;
      cs_real_t *rcodcl3 = f->bc_coeffs->rcodcl3;

      if (f == CS_F_(vel)) {
        for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
          const cs_lnum_t f_id = itrifb[ii];
          if (icodcl[f_id] == 0) {
            icodcl[f_id] = wall_bc_code;

            for (cs_lnum_t k = 0; k < 3; k++) {
              /* rcodcl1_vel[f_id] = User */
              rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
              /* rcodcl3_vel[f_id] unused value or user */
            }
          }
        }
      }

      /* Turbulent variables */
      else if (   f == CS_F_(eps)  || f == CS_F_(rij)   || f == CS_F_(phi)
               || f == CS_F_(omg)  || f == CS_F_(f_bar) || f == CS_F_(k)
               || f == CS_F_(nusa) || f == CS_F_(alp_bl)) {

        if (f->dim == 1) {
          for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
            const cs_lnum_t f_id = itrifb[ii];
            if (icodcl[f_id] == 0) {
              icodcl[f_id] = wall_bc_code;
              rcodcl1[f_id] = 0.;
              rcodcl2[f_id] = cs_math_infinite_r;
              rcodcl3[f_id] = 0.;
            }
          }
        }
        else {
          for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
            const cs_lnum_t f_id = itrifb[ii];
            if (icodcl[f_id] == 0) {
              icodcl[f_id] = wall_bc_code ;
              for (cs_lnum_t k = 0; k < f->dim; k++) {
                rcodcl1[k*n_b_faces+f_id] = 0.;
                rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
                rcodcl3[k*n_b_faces+f_id] = 0.;
              }
            }
          }
        }

      }

      /* Homogeneous Neumann on the pressure */
      else if (f == CS_F_(p)) {
        for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
          const cs_lnum_t f_id = itrifb[ii];
          if (icodcl[f_id] == 0) {
            icodcl[f_id] = 3;
            rcodcl1[f_id] = 0.;
            rcodcl2[f_id] = cs_math_infinite_r;
            rcodcl3[f_id] = 0.;
          }
        }
      }

      /* Turbulent fluxes */

      else if (cs_field_get_key_int(f, keysca) <= 0)
        continue;

      /* Get the turbulent flux model for the scalar */
      const int turb_flux_model = cs_field_get_key_int(f, kturt);
      const int turb_flux_model_type = turb_flux_model / 10;

      if (turb_flux_model_type == 3) {

        cs_field_t *f_tf
          = cs_field_by_composite_name(f->name, "turbulent_flux");

        int *icodcl_tf = f_tf->bc_coeffs->icodcl;
        cs_real_t *rcodcl1_tf = f_tf->bc_coeffs->rcodcl1;

        for (cs_lnum_t jj = s_id; jj < e_id; jj++) {
          const cs_lnum_t f_id = itrifb[jj];
          if (icodcl_tf[f_id] == 0) {
            icodcl_tf[f_id] = wall_bc_code ;
            for (cs_lnum_t k = 0; k < f->dim; k++)
              rcodcl1_tf[k*n_b_faces+f_id] = 0.;
          }
        }

      }

      /* EB-GGDH/AFM/DFM alpha boundary conditions */

      if (turb_flux_model%10 == 1) {

        cs_field_t *f_al = cs_field_by_composite_name(f->name, "alpha");

        int *icodcl_al = f_al->bc_coeffs->icodcl;
        cs_real_t *rcodcl1_al = f_al->bc_coeffs->rcodcl1;

        for (cs_lnum_t jj = s_id; jj < e_id; jj++) {
          const cs_lnum_t f_id = itrifb[jj];
          if (icodcl_al[f_id] == 0) {
            icodcl_al[f_id] = wall_bc_code ;
            rcodcl1_al[f_id] = 0.;
          }
        }

      }

    } /* End of loop on fields */

  }

  /*  When called before time loop, some values are not yet available. */

  if (init) {
    BFT_FREE(bc_type_idx);
    return;
  }

  /* Convert to rcodcl and icodcl (if not already entered by user)
     =============================================================

     Now handle variables for which there is no special processing
     (outside of pressure, velocity ...) */

  /* Inlet and convective inlet bis
     ------------------------------ */

  /* Pressure is already treated (with a Neumann BC)
     Dirichlet BC for the velocity.
     Dirichlet BC for scalars if the user give a value, otherwise
     homogeneous Neumann if the mass flux is outgoing. */

  {
    int err_flags[2] = {0, -1};  /* error indicator, last matching field */

    cs_field_t *f_yplus = cs_field_by_name_try("yplus");
    cs_field_t *f_zground = cs_field_by_name_try("z_ground");

    int inlet_types[2] = {CS_INLET, CS_CONVECTIVE_INLET};
    int inlet_codes[2] = {1, 13};

    for (int field_id = 0; field_id < n_fields; field_id++) {

      cs_field_t *f = cs_field_by_id(field_id);

      if (!(f->type & CS_FIELD_VARIABLE))
        continue;
      if (f->type & CS_FIELD_CDO)
        continue;

      int *icodcl = f->bc_coeffs->icodcl;
      cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1;
      cs_real_t *rcodcl2 = f->bc_coeffs->rcodcl2;
      cs_real_t *rcodcl3 = f->bc_coeffs->rcodcl3;
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);

      /* Convective mass flux of the variable */
      const cs_real_t *b_massflux
        = cs_field_by_id(cs_field_get_key_int(f, kbmasf))->val;

      for (int type_idx = 0; type_idx < 2; type_idx++) {

        int inlet_bc_type = inlet_types[type_idx];
        int inlet_code = inlet_codes[type_idx];

        const cs_lnum_t s_id = bc_type_idx[inlet_bc_type-1];
        const cs_lnum_t e_id = bc_type_idx[inlet_bc_type];

        for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
          const cs_lnum_t f_id = itrifb[ii];

          if (icodcl[f_id] == 0) {

            /* For non convected variables,
               if nothing is defined: homogeneous Neumann */
            if (eqp->iconv == 0 && rcodcl1[f_id] > 0.5*cs_math_infinite_r) {
              icodcl[f_id] = 3;
              for (cs_lnum_t k = 0; k < f->dim; k++) {
                rcodcl1[k*n_b_faces+f_id] = 0.;
                rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
                rcodcl3[k*n_b_faces+f_id] = 0.;
              }
            }

            else if (f == CS_F_(vel)) {
              bool rcodcl1_set = true;
              for (cs_lnum_t k = 0; k < 3; k++) {
                if (rcodcl1[k*n_b_faces + f_id] > 0.5*cs_math_infinite_r)
                  rcodcl1_set = false;
              }
              if (rcodcl1_set == false) {
                bc_type[f_id] = - abs(bc_type[f_id]);
                if (err_flags[0] == 0 || err_flags[0] == 2)
                  err_flags[0] = err_flags[0] + 1;
              }
              else {
                icodcl[f_id] = inlet_code;
                for (cs_lnum_t k = 0; k < 3; k++) {
                  /* rcodcl1[f_id] = given by the user */
                  rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
                  rcodcl3[k*n_b_faces+f_id] = 0.;
                }
              }
            }

            else if (rcodcl1[f_id] > 0.5*cs_math_infinite_r) {
              const cs_real_t flumbf = b_massflux[f_id];
              /* Outgoing flux or yplus or z_ground */
              if (   flumbf >= - cs_math_epzero
                  || f == f_yplus
                  || f == f_zground) {
                icodcl[f_id] = 3;
                for (cs_lnum_t k = 0; k < f->dim; k++) {
                  rcodcl1[k*n_b_faces+f_id] = 0.;
                  rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
                  rcodcl3[k*n_b_faces+f_id] = 0.;
                }
              }
              else {
                bc_type[f_id] = - abs(bc_type[f_id]);
                if (err_flags[0] < 2) {
                  err_flags[0] = err_flags[0] + 2;
                  /* last field id */
                  err_flags[1] = CS_MAX(err_flags[1], f->id);
                }
              }
            }

            else {
              icodcl[f_id] = inlet_code;
              for (cs_lnum_t k = 0; k < f->dim; k++) {
                /* rcodcl1[f_id] = given by the user */
                rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
                rcodcl3[k*n_b_faces+f_id] = 0.;
              }
            }
          }
        }

      } /* End of loop on inlet types (int bc_type_idx) */


    } /* End of loop on fields */

    cs_parall_max(2, CS_INT_TYPE, err_flags);

    if (err_flags[0] > 0) {
      if (err_flags[0] == 1 || err_flags[0] == 3)
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("Error: incorrect or incomplete boundary conditions\n"
             "======\n\n"
             "At least one open boundary face is declared as inlet (or outlet) "
             "with prescribed velocity\n"
             "but its velocity value has not been set for all components.\n"));

      if (err_flags[0] == 2 || err_flags[0] == 3)
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("Error: incorrect or incomplete boundary conditions\n"
             "======\n\n"
             "At least one open boundary face declared as inlet (or outlet)\n"
             "with prescribed velocity for which the value of \"%s\"\n"
             "has not been specified (Dirichlet condition).\n"),
           cs_field_by_id(err_flags[1])->name);

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n"
           "The calculation will not be run.\n\n"
           "Verify the boundary condition definitions.\n"));

      cs_boundary_conditions_error(bc_type, NULL);
    }
  }

  /* Outlet (free inlet outlet) + Free inlet (outlet) Bernoulli + Free surface
     ------------------------------------------------------------------------- */

  /* Pressure has a Dirichlet code, velocities with code 9 have been handled
     earlier.
     The rest is Dirichlet if the user provides a value (inflow or outflow
     mass flux).
     If there is no user-provided values, use a homogeneous  Neumann. */

  /* Free inlet (Bernoulli, a dirichlet is needed on the other
     variables than velocity and pressure, already treated)
     Free outlet (homogeneous Neumann) */

  /* Outlet + Free Inlet + Free surface */

  {
    int f_bc_types[3] = {CS_OUTLET, CS_FREE_INLET, CS_FREE_SURFACE};

    for (int field_id = 0; field_id < n_fields; field_id++) {

      cs_field_t *f = cs_field_by_id(field_id);

      if (!(f->type & CS_FIELD_VARIABLE))
        continue;
      if (f->type & CS_FIELD_CDO)
        continue;

      int *icodcl = f->bc_coeffs->icodcl;
      cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1;
      cs_real_t *rcodcl2 = f->bc_coeffs->rcodcl2;
      cs_real_t *rcodcl3 = f->bc_coeffs->rcodcl3;

      for (int type_idx = 0; type_idx < 3; type_idx++) {

        int f_bc_type = f_bc_types[type_idx];

        const cs_lnum_t s_id = bc_type_idx[f_bc_type-1];
        const cs_lnum_t e_id = bc_type_idx[f_bc_type];

        for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
          const cs_lnum_t f_id = itrifb[ii];

          if (icodcl[f_id] == 0) {

            if (rcodcl1[f_id] > 0.5*cs_math_infinite_r) {
              icodcl[f_id] = 3;
              for (cs_lnum_t k = 0; k < f->dim; k++) {
                rcodcl1[k*n_b_faces+f_id] = 0.;
                rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
                rcodcl3[k*n_b_faces+f_id] = 0.;
              }
            }
            else {
              icodcl[f_id] = 1;
              for (cs_lnum_t k = 0; k < f->dim; k++) {
                /* rcodcl1[k*n_b_faces+f_id] is given by the user */
                rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
                rcodcl3[k*n_b_faces+f_id] = 0.;
              }
            }

          }
        } /* End of loop in faces */

      } /* End of loop on bc_type_idx */
    }
  }

  /* Smooth and rough walls bis
     -------------------------- */

  /* Velocity and turbulence variables have code 5 or 6, handled earlier,
     and the rest is Neumann. */

  for (int wall_bc_type = CS_SMOOTHWALL;
       wall_bc_type <= CS_ROUGHWALL;
       wall_bc_type++) {

    if (_n_bc_faces(wall_bc_type, bc_type_idx) <= 0)
      continue;

    const cs_lnum_t s_id = bc_type_idx[wall_bc_type-1];
    const cs_lnum_t e_id = bc_type_idx[wall_bc_type];

    for (int field_id = 0; field_id < n_fields; field_id++) {

      cs_field_t *f = cs_field_by_id(field_id);

      if (!(f->type & CS_FIELD_VARIABLE))
        continue;
      if (f->type & CS_FIELD_CDO)
        continue;

      int *icodcl = f->bc_coeffs->icodcl;
      cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1;
      cs_real_t *rcodcl2 = f->bc_coeffs->rcodcl2;
      cs_real_t *rcodcl3 = f->bc_coeffs->rcodcl3;

      for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
        const cs_lnum_t f_id = itrifb[ii];

        if (icodcl[f_id] == 0) {
          icodcl[f_id] = 3;
          for (cs_lnum_t k = 0; k < f->dim; k++) {
            rcodcl1[k*n_b_faces+f_id] = 0.;
            rcodcl2[k*n_b_faces+f_id] = cs_math_infinite_r;
            rcodcl3[k*n_b_faces+f_id] = 0.;
          }
        }

      }

    }

  }

  /* Automatic treatment for some variables
     -------------------------------------- */

  /* Put homogeneous Neumann on hydrostatic pressure for iphydr=1
     if not modified by the user */

  if (vp_param->icalhy > 0) {
    cs_field_t *f_hp = cs_field_by_name("hydrostatic_pressure");

    if (   (f_hp->type & CS_FIELD_VARIABLE)
        && (!(f_hp->type & CS_FIELD_CDO))) {

      int *icodcl_hp = f_hp->bc_coeffs->icodcl;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        if (icodcl_hp[f_id] == 0)
          icodcl_hp[f_id] = 3;
      }
    }
  }

  /* Ensure that for all scalars without diffusion
     wall values ignore diffusion. */

  for (int field_id = 0; field_id < n_fields; field_id++) {

    cs_field_t *f = cs_field_by_id(field_id);

    if (!(f->type & CS_FIELD_VARIABLE))
      continue;
    if (f->type & CS_FIELD_CDO)
      continue;

    /* Also ensure values of icodcl are the same for all dimensions */

    if (f->dim > 1) {
      int *icodcl = f->bc_coeffs->icodcl;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        for (cs_lnum_t k = 1; k < f->dim; k++)
          icodcl[k*n_b_faces+f_id] = icodcl[f_id];
      }
    }

    /* Wall values for scalars without diffusion */

    if (cs_field_get_key_int(f, keysca) <= 0)
      continue;

    cs_equation_param_t *eqp = cs_field_get_equation_param(f);

    if (eqp->idiff == 0) {
      int *icodcl = f->bc_coeffs->icodcl;
      cs_real_t *rcodcl3 = f->bc_coeffs->rcodcl3;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        if (icodcl[f_id] == 5 || icodcl[f_id] == 6) {
          icodcl[f_id] = 3;
          for (cs_lnum_t k = 0; k < f->dim; k++)
            rcodcl3[k*n_b_faces+f_id] = 0.;
        }
        else if (icodcl[f_id] == 3) {
          for (cs_lnum_t k = 0; k < f->dim; k++)
            rcodcl3[k*n_b_faces+f_id] = 0.;
        }
      }
    }

  } /* End of loop on fields for wall values ignoring diffusion */

  /* Put homogeneous Neumann on wall distance
     if not modified by the user */

  const cs_field_t *f_wall_dist = cs_field_by_name_try("wall_distance");
  if (f_wall_dist != NULL) {

    if (   (f_wall_dist->type & CS_FIELD_VARIABLE)
        && (!(f_wall_dist->type & CS_FIELD_CDO))) {
      int *icodcl = f_wall_dist->bc_coeffs->icodcl;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        if (icodcl[f_id] == 0)
          icodcl[f_id] = 3;
      }
    }

  }

  /* Matrix diagonal reinforcement if no dirichlet points
     ==================================================== */

  /* We reinforce if istat=0 and if the option is active (idircl=1).
     If one of these conditions is false, we force ndircl to be
     at least 1 to avoid shifting the diagonal. */

  {
    int **ndircl_p;
    cs_gnum_t *ndircl;
    BFT_MALLOC(ndircl_p, n_fields, int *);
    BFT_MALLOC(ndircl, n_fields, cs_gnum_t);

    int n_vars = 0;

    for (int field_id = 0; field_id < n_fields; field_id++) {

      cs_field_t *f = cs_field_by_id(field_id);

      if (!(f->type & CS_FIELD_VARIABLE))
        continue;
      if (f->type & CS_FIELD_CDO)
        continue;

      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      eqp->ndircl = 0;

      if (eqp->istat > 0 || eqp->idircl == 0)
        eqp->ndircl = 1;

      int *icodcl = f->bc_coeffs->icodcl;
      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        if (icodcl[f_id] == 1 || icodcl[f_id] == 5 || icodcl[f_id] == 15)
          eqp->ndircl += 1;
      }

      ndircl_p[n_vars] = &(eqp->ndircl);
      ndircl[n_vars] = eqp->ndircl;

      n_vars++;
    }

    cs_parall_counter(ndircl, n_vars);

    for (int i = 0; i < n_vars; i++) {
      *(ndircl_p[i]) = ndircl[i];
    }

    BFT_FREE(ndircl_p);
    BFT_FREE(ndircl);
  }

  /* Compute and log the mass flow at the different types of faces.
     ==============================================================

     It would be useful to do the logging in cs_log_iteration, but be careful,
     we log the mass flux of the previous time step, while in cs_log_iteration,
     we log at the end of the time step, so an incoherency would be possible.

     It could also be useful to log other values (balances, ...). */

  {
    /* Outputs: mass flux if verbosity or when log is on,
       and at the first two iterations and the two last iterations.
       Only the first iteration on cs_solve_navier_stokes. */

    /* Always print 2 first iterations and the last 2 iterations */
    cs_lnum_t modntl = 1;
    if (nt_cur - nt_prev < 2 || (nt_cur >= nt_max - 1) || eqp_vel->iwarni >= 1)
      modntl = 0;
    else if (cs_log_default_is_active())
      modntl = 0;

    if (modntl == 0) {

      /* Header */
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n"
           "   ** BOUNDARY MASS FLOW INFORMATION\n"
           "      ------------------------------\n\n"));

      cs_real_t *flumty = NULL;
      BFT_MALLOC(flumty, CS_MAX_BC_TYPE, cs_real_t);
      for (int i = 0; i < CS_MAX_BC_TYPE; i++)
        flumty[i] = 0;

      /* Convective mass flux of the velocity */
      const cs_real_t *b_massflux
        = cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf))->val;

      for (cs_lnum_t ii = 0; ii < CS_MAX_BC_TYPE; ii++) {
        const cs_lnum_t s_id = bc_type_idx[ii];
        const cs_lnum_t e_id = bc_type_idx[ii+1];

        for (int jj = s_id; jj < e_id; jj++) {
          const cs_lnum_t f_id = itrifb[jj];
          flumty[ii] = flumty[ii] + b_massflux[f_id];
        }
      }

      cs_log_separator(CS_LOG_DEFAULT);
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("Boundary type          Code    Nb faces           Mass flow\n"));
      cs_log_separator(CS_LOG_DEFAULT);

      cs_gnum_t inb[CS_MAX_BC_TYPE];
      char is_user_type[CS_MAX_BC_TYPE];

      for (int ii = 0; ii < CS_MAX_BC_TYPE; ii++) {
        inb[ii] = bc_type_idx[ii+1] - bc_type_idx[ii];
        is_user_type[ii] = 1;
      }

      cs_parall_counter(inb, CS_MAX_BC_TYPE);
      cs_parall_sum(CS_MAX_BC_TYPE, CS_REAL_TYPE, flumty);

      if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {

        /* Common type ids.
           We will log counts based on these type ids first, in that
           order, so as to make comparisons with prior code versions
           easier. */

        int log_always[11] = {1, 1, 1,
                              1, 1, 1,
                              1,
                              0, 0,
                              1, 1};

        if (cs_sat_coupling_n_couplings() >= 1) {
          if (cs_glob_sat_coupling_face_interpolation_type == 0)
            log_always[7] = 1;
          else
            log_always[8] = 1;
        }

        for (int jj = 0; jj < 11; jj++) {
          const int ii = type_id[jj] - 1;
          if (log_always[jj] || inb[ii] > 0)
            cs_log_printf(CS_LOG_DEFAULT, _("%-17s  %8d%12llu      %18.9e\n"),
                          info[jj], ii+1, (unsigned long long)inb[ii], flumty[ii]);
          is_user_type[ii] = 0;
        }

      }
      else {

        for (int jj = 0; jj < 9; jj++) {
          const int ii = type_id_c[jj] - 1;
          cs_log_printf(CS_LOG_DEFAULT, _("%-17s  %8d%12llu      %18.9e\n"),
                        info_c[jj], ii+1, (unsigned long long)inb[ii],
                        flumty[ii]);
          is_user_type[ii] = 0;
        }

      }

      /* User type */
      for (cs_lnum_t ii = 0; ii < CS_MAX_BC_TYPE; ii++) {
        if (is_user_type[ii] != 1)
          continue;
        cs_gnum_t inb_user = bc_type_idx[ii+1] - bc_type_idx[ii];
        if (inb_user > 0)
          cs_log_printf(CS_LOG_DEFAULT, _("%-17s %8d %12llu %e18.9\n"),
                        "User type", ii+1,
                        (unsigned long long)inb_user, flumty[ii]);
      }

      BFT_FREE(flumty);
      cs_log_separator(CS_LOG_DEFAULT);

    } /* Test on logging (modntl) */
  }

  BFT_FREE(bc_type_idx);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
