/*============================================================================
 * Coal combustion model boundary conditions.
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

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "base/cs_boundary.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_boundary_zone.h"
#include "comb/cs_coal.h"
#include "comb/cs_coal_ht_convert.h"
#include "pprt/cs_combustion_model.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parameters.h"
#include "base/cs_parameters_check.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "comb/cs_coal_boundary_conditions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal_boundary_conditions.cpp

  \brief Coal combustion model boundary conditions.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*! \brief Inlet definition for pulverized coal combustion */

struct _cs_coal_bc_inlet_state_t {

  cs_real_t  x20[CS_COMBUSTION_COAL_MAX_CLASSES];

};

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destructor for coal boundary conditions inlet structure.
 *
 * \param[ci_p]  pointer to associated structure
 *
 * \return: null pointer
 */
/*----------------------------------------------------------------------------*/

static void *
_destroy_inlet(void  *ci_p)
{
  if (ci_p != NULL) {

    cs_coal_bc_inlet_t *ci = (cs_coal_bc_inlet_t *)ci_p;

    CS_FREE(ci->state);
    CS_FREE(ci);

  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_eval_at_location_t function to compute the flow rate at
 *        boundary faces based on given air and coal flow rates.
 *
 * For the calling function, elt_ids is optional. If non-null, array(s) should
 * be accessed with an indirection.
 *
 * \param[in]       location_id   base associated mesh location id
 * \param[in]       n_elts        number of elements to consider
 * \param[in]       elt_ids       list of elements ids
 * \param[in]       input         pointer to cs_boundary_conditions_inlet_t
 * \param[in, out]  vals          resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_qm_from_oxidant(int               location_id,
                 cs_lnum_t         n_elts,
                 const cs_lnum_t  *elt_ids,
                 void             *input,
                 void             *vals_p)
{
  assert(location_id == CS_MESH_LOCATION_NONE);

  cs_coal_bc_inlet_t *ci = (cs_coal_bc_inlet_t *)input;

  cs_real_t  *vals = (cs_real_t *)vals_p;

  cs_coal_model_t *coal = cs_glob_coal_model;

  if (ci->qm_air_func != NULL) {
    cs_real_t qm_air[1] = {0};
    ci->qm_air_func(location_id,
                    n_elts,
                    elt_ids,
                    ci->qm_air_input,
                    qm_air);
    ci->qm_air = qm_air[0];
  }

  vals[0] = ci->qm_air;

  for (int i = 0; i < coal->n_coals; i++)
    vals[0] += ci->qimpcp[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Recompute inlet state for restart.
 */
/*----------------------------------------------------------------------------*/

static void
_recompute_inlet_state(void)
{
  /* Initialization */

  cs_coal_model_t *cm = cs_glob_coal_model;

  const int n_coals = cm->n_coals;

  /* Loop on coal inlet boundaries */

  assert(cs_glob_boundaries != NULL);

  const cs_boundary_t *bdy = cs_glob_boundaries;

  for (int bdy_idx = 0; bdy_idx < bdy->n_boundaries; bdy_idx += 1) {

    if (! (bdy->types[bdy_idx] & CS_BOUNDARY_INLET))
      continue;

    const cs_zone_t *z = cs_boundary_zone_by_id(bdy->zone_ids[bdy_idx]);

    auto ci = reinterpret_cast<cs_coal_bc_inlet_t *>
                (cs_boundary_conditions_get_model_inlet(z));

    /* An input ientre must be of type ientat = 1 or ientcp = 1
       -------------------------------------------------------- */

    if (ci->ientat == 1 || ci->ientcp == 1) {

      cs_real_t qimpc = cs_boundary_conditions_open_get_mass_flow_rate(z);

      int i_shift = 0;

      for (int icha = 0; icha < n_coals; icha++) {

        for (int iclapc = 0; iclapc < cm->n_classes_per_coal[icha]; iclapc++) {

          int icla = iclapc + i_shift;

          /* Calculation of total X2 per zone
             Small correction in case of a closed inlet */

          if (cs::abs(qimpc) < cs_math_epzero)
            ci->state->x20[icla] = 0.;
          else
            ci->state->x20[icla] =   ci->qimpcp[icha] / qimpc
                                   * ci->distch[icha][iclapc] * 1.e-2;

        }

        i_shift += cm->n_classes_per_coal[icha];

      }

    }

  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to coal boundary conditions inlet structure.
 *
 * If no such structure was previously present, it is created and linked
 * to the matching open boundary condition inlet.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone
 */
/*----------------------------------------------------------------------------*/

cs_coal_bc_inlet_t *
cs_coal_boundary_conditions_get_inlet(const  cs_zone_t   *zone)
{
  auto ci = reinterpret_cast<cs_coal_bc_inlet_t *>
              (cs_boundary_conditions_get_model_inlet(zone));

  /* Add and initialize coal inlet if not present */

  if (ci == NULL) {
    CS_MALLOC(ci, 1, cs_coal_bc_inlet_t);

    ci->zone = zone;

    ci->ientat = 0;
    ci->ientcp = 0;
    ci->inmoxy = 0;

    for (int i = 0; i < CS_COMBUSTION_MAX_COALS; i++) {
      ci->qimpcp[i] = 0;
      ci->timpcp[i] = 0;
      for (int j = 0; j < CS_COMBUSTION_MAX_CLASSES_PER_COAL; j++)
        ci->distch[i][j] = 0;
    }

    ci->t_air = 0;

    ci->qm_air = nan("");  /* Set to usable value if used */
    ci->qm_air_func = NULL;
    ci->qm_air_input = NULL;

    /* Add state-tracking structure */

    CS_MALLOC(ci->state, 1, cs_coal_bc_inlet_state_t);
    for (int i = 0; i < CS_COMBUSTION_COAL_MAX_CLASSES; i++)
      ci->state->x20[i] = 0;

    /* Now register in high-level boundary conditions mechanism */

    cs_boundary_conditions_assign_model_inlet(zone,
                                              (void *)ci,
                                              (void *)_destroy_inlet);

    /* Link with boundary mechanism;
       This will allow simplifying later loops */

    assert(cs_glob_boundaries != NULL);

    cs_boundary_t *bdy = cs_glob_boundaries;
    if (cs_boundary_id_by_zone_id(bdy, zone->id) < 0)
      cs_boundary_add(bdy, CS_BOUNDARY_INLET, zone->name);

  }

  return (ci);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant air mass flow rate to an inlet.
 *
 * The total mass flow rate will also include that of the pulverized coals.
 *
 * This is otherwise similar to
 * \ref cs_boundary_conditions_open_set_mass_flow_rate_by_value.
 *
 * \param[in]  z  pointer to associated zone
 * \param[in]  q  associated constant mass flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_boundary_conditions_inlet_set_air_mass_flow_rate_by_value
  (const  cs_zone_t  *z,
   cs_real_t          q)
{
  cs_coal_bc_inlet_t *ci = cs_coal_boundary_conditions_get_inlet(z);

  ci->qm_air = q;

  cs_boundary_conditions_open_set_mass_flow_rate_by_func(z,
                                                         _qm_from_oxidant,
                                                         ci);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign an air mass flow rate to an inlet based on provided function.
 *
 * The total mass flow rate will also include that of the pulverized coals.
 *
 * This is otherwise similar to
 * \ref cs_boundary_conditions_open_set_mass_flow_rate_by_func.
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated scalar (mass flow rate) evaluation function
 * \param[in]  input  optional function evaluation input, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_boundary_conditions_inlet_set_air_mass_flow_rate_by_func
  (const  cs_zone_t       *z,
   cs_eval_at_location_t  *func,
   void                   *input)
{
  cs_coal_bc_inlet_t *ci = cs_coal_boundary_conditions_get_inlet(z);

  ci->qm_air_func = func;
  ci->qm_air_input = input;

  cs_boundary_conditions_open_set_mass_flow_rate_by_func(z,
                                                         _qm_from_oxidant,
                                                         ci);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic boundary condition for pulverized coal combustion.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_boundary_conditions(int  bc_type[])
{
  /* Initialization
   * ============== */

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_coal_model_t *cm = cs_glob_coal_model;

  const int n_coals = cm->n_coals;

  cs_field_t *f = NULL;

  cs_real_t *b_x1 = cs_field_by_name("b_x_c")->val;

  cs_real_t *rcodcl1_age = NULL;
  cs_real_t **age_rcodcl1 = NULL;
  {
    const cs_field_t *f_bulk_age = cs_field_by_name_try("age");

    if (f_bulk_age != NULL) {
      CS_MALLOC(age_rcodcl1, cm->nclacp, cs_real_t *);

      for (int icla = 0; icla < cm->nclacp; icla++) {
        char name[20];
        snprintf(name, 19, "n_p_age_%02d", icla+1);
        name[19] = '\0';
        cs_field_t *f_p_age = cs_field_by_name(name);
        age_rcodcl1[icla] = f_p_age->bc_coeffs->rcodcl1;
      }

      rcodcl1_age = f_bulk_age->bc_coeffs->rcodcl1;
    }
  }

  /* Boundary conditions for HM */
  cs_real_t *rcodcl1_h = CS_F_(h)->bc_coeffs->rcodcl1;

  /* Boundary condition for x1*h1 */
  cs_real_t *rcodcl1_x_c_h = cs_field_by_name("x_c_h")->bc_coeffs->rcodcl1;

  /* Boundary conditions for X1.F4M (oxyd 2) */
  cs_real_t *rcodcl1_f4m = NULL;
  f = cs_field_by_name_try("fr_oxyd2");
  if (f != NULL)
    rcodcl1_f4m = f->bc_coeffs->rcodcl1;

  /* Boundary conditions for X1.F5M (oxyd 3) */
  cs_real_t *rcodcl1_f5m = NULL;
  f = cs_field_by_name_try("fr_oxyd3");
  if (f != NULL)
    rcodcl1_f5m = f->bc_coeffs->rcodcl1;

  /* Boundary conditions for X1.F6M (water) */
  cs_real_t *rcodcl1_f6m = NULL;
  f = cs_field_by_name_try("fr_h2o");
  if (f != NULL)
    rcodcl1_f6m = f->bc_coeffs->rcodcl1;

  /* Boundary conditions for X1.F7M_O2 */
  cs_real_t *rcodcl1_f7m = cs_field_by_name("fr_het_o2")->bc_coeffs->rcodcl1;

  /* Boundary conditions for X1.F8M.CO2 */
  cs_real_t *rcodcl1_f8m = NULL;
  f = cs_field_by_name_try("fr_het_co2");
  if (f != NULL)
    rcodcl1_f8m = f->bc_coeffs->rcodcl1;

  /* Boundary conditions for X1.F9M.H2O */
  cs_real_t *rcodcl1_f9m = NULL;
  f = cs_field_by_name_try("fr_het_h2o");
  if (f != NULL)
    rcodcl1_f9m = f->bc_coeffs->rcodcl1;

  /* Boundary conditions for X1.Variance */
  cs_real_t *rcodcl1_fvp2m
    = cs_field_by_name("f1f2_variance")->bc_coeffs->rcodcl1;

  /* Boundary conditions for X1.YCO2 */
  cs_real_t *rcodcl1_iyco2 = NULL;
  if (cm->ieqco2 == 1)
    rcodcl1_iyco2 = cs_field_by_name("x_c_co2")->bc_coeffs->rcodcl1;

  /* Boundary conditions for X1.HCN, X1.NO, T air */
  cs_real_t *rcodcl1_iyhcn = NULL;
  cs_real_t *rcodcl1_iyno = NULL;
  cs_real_t *rcodcl1_iynh3 = NULL;
  cs_real_t *rcodcl1_ihox = NULL;
  if (cm->ieqnox == 1) {
    rcodcl1_iyhcn = cs_field_by_name("x_c_hcn")->bc_coeffs->rcodcl1;
    rcodcl1_iyno = cs_field_by_name("x_c_no")->bc_coeffs->rcodcl1;
    rcodcl1_iynh3 = cs_field_by_name("x_c_nh3")->bc_coeffs->rcodcl1;
    rcodcl1_ihox = cs_field_by_name("x_c_h_ox")->bc_coeffs->rcodcl1;
  }

  /* Loop on coal inlet boundaries
     ============================= */

  assert(cs_glob_boundaries != NULL);

  const cs_boundary_t *bdy = cs_glob_boundaries;

  for (int bdy_idx = 0; bdy_idx < bdy->n_boundaries; bdy_idx += 1) {

    if (! (bdy->types[bdy_idx] & CS_BOUNDARY_INLET))
      continue;

    const cs_zone_t *z = cs_boundary_zone_by_id(bdy->zone_ids[bdy_idx]);

    const cs_lnum_t n_elts = z->n_elts;
    const cs_lnum_t *elt_ids = z->elt_ids;

    auto ci = reinterpret_cast<cs_coal_bc_inlet_t *>
                (cs_boundary_conditions_get_model_inlet(z));

    cs_real_t xsolid[CS_COMBUSTION_COAL_MAX_SOLIDS];
    for (int isol = 0; isol < CS_COMBUSTION_COAL_MAX_SOLIDS; isol++)
      xsolid[isol] = 0.;

    cs_real_t h1 = 0, h2[CS_COMBUSTION_MAX_CLASSES_PER_COAL];
    cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

    const int ico2 = cm->ico2 - 1;
    const int ih2o = cm->ih2o - 1;
    const int in2  = cm->in2 - 1;
    const int io2  = cm->io2 - 1;

    cs_real_t *rcodcl1_u = CS_F_(vel)->bc_coeffs->rcodcl1;
    cs_real_t *rcodcl1_v = rcodcl1_u + n_b_faces;
    cs_real_t *rcodcl1_w = rcodcl1_v + n_b_faces;

    cs_lnum_t i_shift = 0;
    cs_real_t x20t = 0, x2h20t = 0;

    /* Verify that coal distribution sum = 100% for area with ientcp = 1
       ----------------------------------------------------------------- */

    if (ci->ientcp == 1) {
      for (int icha = 0; icha < n_coals; icha++) {
        cs_real_t totcp = 0;
        for (int iclapc = 0; iclapc < cm->n_classes_per_coal[icha]; iclapc++)
          totcp += ci->distch[icha][iclapc];
        if (fabs(totcp - 100) > cs_math_epzero) {
          cs_parameters_error_header(CS_ABORT_DELAYED,
                                     "in pulverized coal inlet definitions");
          cs_log_t log = CS_LOG_DEFAULT;
          cs_log_printf(log,
                        _("\n"
                          "  zone: %s\n"
                          "  coal: %d\n"
                          "  distribution (in percentage):\n"), z->name, icha);
          for (int iclapc = 0; iclapc < cm->n_classes_per_coal[icha]; iclapc++)
            cs_log_printf(log,
                          _("    class %02d: %e12.5\n\n"),
                          iclapc+1, ci->distch[icha][iclapc]);
         cs_log_printf(log,
                        _("  sum of distributions: %g\n (100 expected)\n"),
                        totcp);
          cs_parameters_error_footer(CS_ABORT_DELAYED);
        }
      }
      cs_parameters_error_barrier();
    }

    /* An input ientre must be of type ientat = 1 or ientcp = 1
       -------------------------------------------------------- */

    if (ci->ientat == 1 || ci->ientcp == 1) {

      cs_real_t qimpc = cs_boundary_conditions_open_get_mass_flow_rate(z);

      i_shift = 0;

      for (int icha = 0; icha < n_coals; icha++) {

        for (int iclapc = 0; iclapc < cm->n_classes_per_coal[icha]; iclapc++) {

          int icla = iclapc + i_shift;

          /* Calculation of total X2 per zone
             Small correction in case of a closed inlet */

          if (cs::abs(qimpc) < cs_math_epzero)
            ci->state->x20[icla] = 0.;
          else
            ci->state->x20[icla] =   ci->qimpcp[icha] / qimpc
                                   * ci->distch[icha][iclapc] * 1.e-2;

          x20t += ci->state->x20[icla];

          /* Compute H2 of class icla */

          cs_real_t t2;

          int ch_id = cm->ich[icha] - 1;
          int ck_id = cm->ick[icha] - 1;
          int ash_id = cm->iash[icha] - 1;
          int wat_id = cm->iwat[icha] - 1;

          if (ci->ientcp == 1) {
            t2 = ci->timpcp[icha];
            xsolid[ch_id] = 1.0 - cm->xashch[icha];
            xsolid[ck_id] = 0;
            xsolid[ash_id] = cm->xashch[icha];

            /* Taking into account humidity */
            if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] == 1) {
              xsolid[ch_id] -= cm->xwatch[icha];
              xsolid[wat_id] = cm->xwatch[icha];
            }
            else
              xsolid[wat_id] = 0.;
          }
          else {
            t2 = ci->t_air;

            xsolid[ch_id] = 1.0 - cm->xashch[icha] - cm->xwatch[icha];
            xsolid[ck_id] = 0;
            xsolid[ash_id] = cm->xashch[icha];
            xsolid[wat_id] = cm->xwatch[icha];
          }

          h2[icla] = cs_coal_ht_convert_t_to_h_particles_by_yi(t2, icla, xsolid);
          x2h20t += ci->state->x20[icla] * h2[icla];

        }

        i_shift += cm->n_classes_per_coal[icha];

      }

      /* Compute H1 */
      for (int ige = 0;
           ige < CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS;
           ige++)
        coefe[ige] = 0.;

      int ioxy = ci->inmoxy - 1;
      cs_real_t dmas =   cm->wmole[io2]  * cm->oxyo2[ioxy]
                       + cm->wmole[in2]  * cm->oxyn2[ioxy]
                       + cm->wmole[ih2o] * cm->oxyh2o[ioxy]
                       + cm->wmole[ico2] * cm->oxyco2[ioxy];

      coefe[io2]  = cm->wmole[io2]  * cm->oxyo2[ioxy ]/dmas;
      coefe[ih2o] = cm->wmole[ih2o] * cm->oxyh2o[ioxy]/dmas;
      coefe[ico2] = cm->wmole[ico2] * cm->oxyco2[ioxy]/dmas;
      coefe[in2]  = cm->wmole[in2]  * cm->oxyn2[ioxy ]/dmas;

      h1 = cs_coal_ht_convert_t_to_h_gas_by_yi(ci->t_air, coefe);

    }

    /* Set boundary conditions for inlet zone faces
       -------------------------------------------- */

    i_shift = 0;

    for (int icha = 0; icha < n_coals; icha++) {

      char f_name[32];

      for (int iclapc = 0; iclapc < cm->n_classes_per_coal[icha]; iclapc++) {

        int icla = iclapc + i_shift;

        /* Boundary conditions for Xch of class icla */

        snprintf(f_name, 31, "x_p_coal_%02d", icla+1); f_name[31] = '\0';
        cs_real_t *rcodcl1_xch = cs_field_by_name(f_name)->bc_coeffs->rcodcl1;
        cs_real_t xch_in =   ci->state->x20[icla]
                           * (1. - cm->xashch[icha]);

        /* Taking into account humidity */
        if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] == 1) {
          xch_in =   ci->state->x20[icla]
                   * (1. - cm->xashch[icha]- cm->xwatch[icha]);
        }

        /* Boundary conditions for Xck of class icla */

        snprintf(f_name, 31, "x_p_char_%02d", icla+1); f_name[31] = '\0';
        cs_real_t *rcodcl1_xck = cs_field_by_name(f_name)->bc_coeffs->rcodcl1;

        /* Boundary conditions for Np of class icla */

        snprintf(f_name, 31, "n_p_%02d", icla+1); f_name[31] = '\0';
        cs_real_t *rcodcl1_inp = cs_field_by_name(f_name)->bc_coeffs->rcodcl1;

        cs_real_t inp_in = ci->state->x20[icla] / cm->xmp0[icla];

        /* Boundary conditions for Xwater of class icla */

        cs_real_t xwt = 0;
        cs_real_t *rcodcl1_xwt =  NULL;
        if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] == 1) {
          snprintf(f_name, 31, "x_p_wt_%02d", icla+1); f_name[31] = '\0';
          rcodcl1_xwt = cs_field_by_name(f_name)->bc_coeffs->rcodcl1;
          xwt = ci->state->x20[icla] * cm->xwatch[icha];
        }

        /* Boundary conditions for H2 of class icla */

        snprintf(f_name, 31, "x_p_h_%02d", icla+1); f_name[31] = '\0';
        cs_real_t *rcodcl1_ih2 = cs_field_by_name(f_name)->bc_coeffs->rcodcl1;
        cs_real_t ih2_in = ci->state->x20[icla] * h2[icla];

        /* Boundary conditions for particle age */

        cs_real_t *rcodcl1_iage = NULL;
        cs_real_t *rcodcl1_v_p_x = NULL;
        cs_real_t *rcodcl1_v_p_y = NULL;
        cs_real_t *rcodcl1_v_p_z = NULL;

        if (cm->idrift >= 1) {
          rcodcl1_iage = age_rcodcl1[icla];

          if (cm->idrift == 1) {
            snprintf(f_name, 31, "v_x_p_%02d", icla+1); f_name[31] = '\0';
            rcodcl1_v_p_x = cs_field_by_name(f_name)->bc_coeffs->rcodcl1;
            snprintf(f_name, 31, "v_y_p_%02d", icla+1); f_name[31] = '\0';
            rcodcl1_v_p_y = cs_field_by_name(f_name)->bc_coeffs->rcodcl1;
            snprintf(f_name, 31, "v_z_p_%02d", icla+1); f_name[31] = '\0';
            rcodcl1_v_p_z = cs_field_by_name(f_name)->bc_coeffs->rcodcl1;
          }
        }

        /* Apply to zone inlet faces */

        for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
          cs_lnum_t elt_id = elt_ids[elt_idx];
          if (   bc_type[elt_id] == CS_INLET
              || bc_type[elt_id] == CS_CONVECTIVE_INLET) {

            rcodcl1_xch[elt_id] = xch_in;
            rcodcl1_xck[elt_id] = 0;
            rcodcl1_inp[elt_id] = inp_in;
            if (rcodcl1_xwt !=  NULL)
              rcodcl1_xwt[elt_id] = xwt;
            rcodcl1_ih2[elt_id] = ih2_in;
            if (rcodcl1_iage != NULL)
              rcodcl1_iage[elt_id] = 0;

            if (rcodcl1_v_p_x != NULL) {
              rcodcl1_v_p_x[elt_id] = rcodcl1_u[elt_id];
              rcodcl1_v_p_y[elt_id] = rcodcl1_v[elt_id];
              rcodcl1_v_p_z[elt_id] = rcodcl1_w[elt_id];
            }

          }
        }

      } /* loop on coal classes */

      i_shift += cm->n_classes_per_coal[icha];

      /* Boundary conditions for x1f1m and x1f2m from coal icha */

      snprintf(f_name, 31, "fr_mv1_%02d", icha+1); f_name[31] = '\0';
      cs_real_t *rcodcl1_fm1 = cs_field_by_name(f_name)->bc_coeffs->rcodcl1;

      snprintf(f_name, 31, "fr_mv2_%02d", icha+1); f_name[31] = '\0';
      cs_real_t *rcodcl1_fm2 = cs_field_by_name(f_name)->bc_coeffs->rcodcl1;

      for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
        cs_lnum_t elt_id = elt_ids[elt_idx];
        if (   bc_type[elt_id] == CS_INLET
            || bc_type[elt_id] == CS_CONVECTIVE_INLET) {
          rcodcl1_fm1[elt_id] = 0.;
          rcodcl1_fm2[elt_id] = 0.;
        }
      }

    } /* loop on coals */

    /* Boundary conditions for HM */
    cs_real_t h_in = (1.-x20t)* h1 + x2h20t;

    /* Boundary condition for x1*h1 */
    cs_real_t hgas_in = (1.-x20t)* h1;

    /* Boundary conditions for X1.F4M (oxyd 2) */

    cs_real_t f4m_in = 0;
    if (ci->inmoxy == 2)
      f4m_in = 1. - x20t;

    /* Boundary conditions for X1.F5M (oxyd 3) */

    cs_real_t f5m_in = 0;
    if (ci->inmoxy == 3)
      f5m_in = 1. - x20t;

    /* Boundary conditions for X1.YCO2 */

    cs_real_t yco2_in = 0;

    if (cm->ieqco2 == 1) {
      int ioxy = ci->inmoxy - 1;
      cs_real_t dmas =   cm->wmole[io2]  * cm->oxyo2[ioxy]
                       + cm->wmole[in2]  * cm->oxyn2[ioxy]
                       + cm->wmole[ih2o] * cm->oxyh2o[ioxy]
                       + cm->wmole[ico2] * cm->oxyco2[ioxy];

      cs_real_t xco2 = cm->oxyco2[ioxy] * cm->wmole[ico2] / dmas;
      yco2_in = xco2 * (1. - x20t);
    }

    cs_real_t hox_in = 0;
    if (cm->ieqnox == 1)
      hox_in = (1. - x20t) * h1;

    /* Apply to zone inlet faces */

    for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
      cs_lnum_t elt_id = elt_ids[elt_idx];
      if (   bc_type[elt_id] == CS_INLET
          || bc_type[elt_id] == CS_CONVECTIVE_INLET) {

        if (rcodcl1_age != NULL)        /* age */
          rcodcl1_age[elt_id] = 0.;

        rcodcl1_h[elt_id] = h_in;
        rcodcl1_x_c_h[elt_id] = hgas_in;

        /* Store the Boundary value of X1 */
        b_x1[elt_id] = 1. - x20t;

        if (rcodcl1_f4m != NULL)
          rcodcl1_f4m[elt_id] = f4m_in;
        if (rcodcl1_f5m != NULL)
          rcodcl1_f5m[elt_id] = f5m_in;
        if (rcodcl1_f6m != NULL)
          rcodcl1_f6m[elt_id] = 0.;     /* water */
        rcodcl1_f7m[elt_id] = 0.;
        if (rcodcl1_f8m != NULL)
          rcodcl1_f8m[elt_id] = 0.;     /* CO2 */
        if (rcodcl1_f9m != NULL)
          rcodcl1_f9m[elt_id] = 0.;     /* H2O */

        rcodcl1_fvp2m[elt_id] = 0.;

        if (cm->ieqco2 == 1)
          rcodcl1_iyco2[elt_id] = yco2_in;

        if (cm->ieqnox == 1) {
          rcodcl1_iyhcn[elt_id] = 0.;
          rcodcl1_iyno[elt_id] = 0.;
          rcodcl1_iynh3[elt_id] = 0.;
          rcodcl1_ihox[elt_id] = hox_in;
        }

      }
    } /* loop on zone faces */

  } /* loop on zones */

  /* Wall BCs on the particle velocity: zero Dirichlet
     ------------------------------------------------- */

  /* TODO: merge v_p_x, v_p_y, and v_p_z fields into a single vector field */

  if (cm->idrift == 1) {

    for (int icha = 0; icha < n_coals; icha++) {

      for (int iclapc = 0; iclapc < cm->n_classes_per_coal[icha]; iclapc++) {

        int *icodcl_v_p_x = NULL;
        int *icodcl_v_p_y = NULL;
        int *icodcl_v_p_z = NULL;
        cs_real_t *rcodcl1_v_p_x = NULL;
        cs_real_t *rcodcl1_v_p_y = NULL;
        cs_real_t *rcodcl1_v_p_z = NULL;

        char f_name[32];

        snprintf(f_name, 31, "v_x_p_%02d", iclapc+1); f_name[31] = '\0';
        f = cs_field_by_name(f_name);
        icodcl_v_p_x = f->bc_coeffs->icodcl;
        rcodcl1_v_p_x = f->bc_coeffs->rcodcl1;

        snprintf(f_name, 31, "v_y_p_%02d", iclapc+1); f_name[31] = '\0';
        f = cs_field_by_name(f_name);
        icodcl_v_p_y = f->bc_coeffs->icodcl;
        rcodcl1_v_p_y = f->bc_coeffs->rcodcl1;

        snprintf(f_name, 31, "v_z_p_%02d", iclapc+1); f_name[31] = '\0';
        f = cs_field_by_name(f_name);
        icodcl_v_p_z = f->bc_coeffs->icodcl;
        rcodcl1_v_p_z = f->bc_coeffs->rcodcl1;

        for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
          if (   bc_type[face_id] == CS_SMOOTHWALL
              || bc_type[face_id] == CS_ROUGHWALL) {

            icodcl_v_p_x[face_id] = 1;
            icodcl_v_p_y[face_id] = 1;
            icodcl_v_p_z[face_id] = 1;

            rcodcl1_v_p_x[face_id] = 0;
            rcodcl1_v_p_y[face_id] = 0;
            rcodcl1_v_p_z[face_id] = 0;

          }
        }

      } /* loop on classes */

    } /* loop on coals */

  } /* case with drift */

  /* Cleanup */

  CS_FREE(age_rcodcl1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute density at boundaries for pulverized coal combustion.
 *
 * This is based on boundary condition definitions, but is called at an
 * earlier stage in the time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_boundary_conditions_density(void)
{
  /* Initialization
   * ============== */

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

  cs_coal_model_t *cm = cs_glob_coal_model;

  const cs_real_t p0 = cs_glob_fluid_properties->p0;

  const cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *brom = CS_F_(rho_b)->val;

  /* Mass density on edges for all faces */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    cs_lnum_t cell_id = b_face_cells[face_id];
    brom[face_id] = crom[cell_id];
  }

  /* Recompute density at coal inlets
     (except before or at first time step where x20 is not known yet) */

  if (cs_glob_time_step->nt_cur <= 1)
    return;

  /* Recompute state in case of restart */

  else if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_prev+1)
    _recompute_inlet_state();

  assert(cs_glob_boundaries != NULL);

  const cs_boundary_t *bdy = cs_glob_boundaries;

  const int ico2 = cm->ico2 - 1;
  const int ih2o = cm->ih2o - 1;
  const int in2  = cm->in2 - 1;
  const int io2  = cm->io2 - 1;

  /* loop on boundary zones, ignore non-inlet zones */

  for (int bdy_idx = 0; bdy_idx < bdy->n_boundaries; bdy_idx += 1) {

    if (! (bdy->types[bdy_idx] & CS_BOUNDARY_INLET))
      continue;

    const cs_zone_t *z = cs_boundary_zone_by_id(bdy->zone_ids[bdy_idx]);

    auto ci = reinterpret_cast<cs_coal_bc_inlet_t *>
                (cs_boundary_conditions_get_model_inlet(z));

    /* An input ientre must be of type ientat = 1 or ientcp = 1
       -------------------------------------------------------- */

    if (ci->ientat != 1 && ci->ientcp != 1)
      continue;

    const cs_lnum_t n_faces = z->n_elts;
    const cs_lnum_t *face_ids = z->elt_ids;

    cs_real_t x20t = 0, x20drho20 = 0;

    for (int iclapc = 0; iclapc < cm->nclacp; iclapc++) {
      x20drho20 += ci->state->x20[iclapc] / cm->rho20[iclapc];
      x20t      += ci->state->x20[iclapc];
    }

    int ioxy = ci->inmoxy - 1;
    cs_real_t dmas =   cm->wmole[io2]  * cm->oxyo2[ioxy]
                     + cm->wmole[in2]  * cm->oxyn2[ioxy]
                     + cm->wmole[ih2o] * cm->oxyh2o[ioxy]
                     + cm->wmole[ico2] * cm->oxyco2[ioxy];

    cs_real_t wmolme = (  cm->oxyo2[ioxy]
                        + cm->oxyn2[ioxy]
                        + cm->oxyh2o[ioxy]
                        + cm->oxyco2[ioxy]) / dmas;

    cs_real_t unsro1 = (wmolme * cs_physical_constants_r * ci->t_air) / p0;
    cs_real_t x1sro1 = (1.0 - x20t) * unsro1;

    cs_real_t rho_b_in = 1.0 / (x1sro1 + x20drho20);

    for (cs_lnum_t i = 0; i < n_faces; i++) {
      cs_lnum_t face_id = face_ids[i];
      brom[face_id] = rho_b_in;
    }

  } /* loop on boundary zones */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
