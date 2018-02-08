/*============================================================================
 * Define scaling parameter for electric model
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_math.h"
#include "cs_mesh_quantities.h"
#include "cs_elec_model.h"
#include "bft_mem.h"
#include "bft_printf.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_parall.h"
#include "cs_physical_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_electric_scaling.c
 *
 * \brief Define scaling parameter for electric model.
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Rescale all electro-magnetic physical fields
 *        (electric potential, current density and Joule effect).
 *
 * \param[in] mesh pointer to a cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to a cs_mesh_quantities_t structure
 * \param[in] dt pointer to a \ref cs_real_t
 *
 * These options allow defining the time step synchronization policy,
 * as well as a time step multiplier.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_scaling_elec(const cs_mesh_t             *mesh,
                     const cs_mesh_quantities_t  *mesh_quantities,
                           cs_real_t             *dt)
{
  /*! [electric_scaling] */

  cs_lnum_t  ncel   = mesh->n_cells;
  cs_lnum_t  ncelet = mesh->n_cells_with_ghosts;
  cs_real_t *xyzcen =  mesh_quantities->cell_cen;
  cs_real_t *volume = mesh_quantities->cell_vol;
  cs_lnum_t  nfac   = mesh->n_i_faces;
  const cs_real_3_t *surfac = (const cs_real_3_t *) mesh_quantities->b_face_normal;
  const cs_real_3_t *cdgfac = (const cs_real_3_t *) mesh_quantities->i_face_cog;

  cs_elec_option_t *elec_opt = cs_get_glob_elec_option();
  const int kivisl = cs_field_key_id("scalar_diffusivity_id");

  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  /* example of a restrike arc */
  if (ielarc >= 1) {
    if (cs_glob_time_step->nt_cur <= 200)
      elec_opt->couimp = 200.;
    else if (cs_glob_time_step->nt_cur > 200 &&
             cs_glob_time_step->nt_cur <= 400)
      elec_opt->couimp = 200. +  2. * (cs_glob_time_step->nt_cur - 200);
    else
      elec_opt->couimp = 600.;
  }

  if (cs_glob_time_step->nt_cur < 400 ||
      cs_glob_time_step->nt_cur == cs_glob_time_step->nt_prev + 1)
    elec_opt->irestrike = 0;

  double econs = 1.5e5;
  double coepot = 0.;
  double coepoa = 1.;

  if (cs_glob_elec_option->irestrike) {
    double amex = 1.e30;
    double aiex = -1.e30;
    double emax = 0.;
    double *w1;
    BFT_MALLOC(w1, ncelet, double);
    int diff_id = cs_field_get_key_int(CS_FI_(curre, 0), kivisl);
    cs_field_t *c_prop = NULL;
    c_prop = cs_field_by_id(diff_id);

    for (int iel = 0; iel < ncel; iel++) {
      double xelec = CS_FI_(curre, 0)->val[iel] / c_prop->val[iel];
      double yelec = CS_FI_(curre, 1)->val[iel] / c_prop->val[iel];
      double zelec = CS_FI_(curre, 2)->val[iel] / c_prop->val[iel];

      w1[iel] = pow(xelec * xelec + yelec * yelec + zelec * zelec, 0.5);
      amex = CS_MIN(amex, w1[iel]);
      aiex = CS_MAX(amex, w1[iel]);
    }
    cs_parall_min(1, CS_DOUBLE, &amex);
    cs_parall_max(1, CS_DOUBLE, &aiex);

    bft_printf("min and max for E : %14.5E %15.4E\n", amex, aiex);

    if (aiex > econs) {
      elec_opt->irestrike  = 1;
      elec_opt->ntdcla = cs_glob_time_step->nt_cur;

      /* initialize restrike point coordinates */
      elec_opt->restrike_point[0] = 1.e-8;
      elec_opt->restrike_point[1] = 1.e-8;
      elec_opt->restrike_point[2] = 1.e-8;
      double diff = 0.;
      double xyzmax[3] = {-1.e10, -1.e10, -1.e10};

      for (int iel = 0; iel < ncel; iel++) {
        diff = aiex - w1[iel];

        if (diff < 1.e-6) {
          emax = w1[iel];
          xyzmax[0] = xyzcen[3 * iel    ];
          xyzmax[1] = xyzcen[3 * iel + 1];
          xyzmax[2] = xyzcen[3 * iel + 2];
        }
      }
      /* we can only have a single restrike point */
      cs_parall_max_loc_vals(3, &emax, xyzmax);

      elec_opt->restrike_point[0] = xyzmax[0];
      elec_opt->restrike_point[1] = xyzmax[1];
      elec_opt->restrike_point[2] = xyzmax[2];

      bft_printf("restrike point : %14.5E %14.5E %14.5E\n",
                 elec_opt->restrike_point[0],
                 elec_opt->restrike_point[1],
                 elec_opt->restrike_point[2]);
    }

    BFT_FREE(w1);

    if (cs_glob_time_step->nt_cur <= elec_opt->ntdcla + 30) {
      double z1 = elec_opt->restrike_point[0] - 3.e-4;
      double z2 = elec_opt->restrike_point[0] + 3.e-4;
      if (z1 < 0.)
        z1 = 0.;
      if (z2 > 2.e-2)
        z2 = 2.e-2;
      for (int iel = 0; iel < ncel; iel++) {
        if (xyzcen[3 * iel + 2] > z1 && xyzcen[3 * iel + 2] < z2) {
          double rayo = elec_opt->restrike_point[0] * xyzcen[3 * iel    ]
                      - elec_opt->restrike_point[1] * xyzcen[3 * iel + 1];
          double denom = pow(elec_opt->restrike_point[0] * elec_opt->restrike_point[0]
                           + elec_opt->restrike_point[1] * elec_opt->restrike_point[1], 0.5);
          rayo /= denom;
          rayo += (xyzcen[3 * iel + 2] - elec_opt->restrike_point[2]
                 * xyzcen[3 * iel + 2] - elec_opt->restrike_point[2]);
          rayo = pow(rayo, 0.5);

          double posi = elec_opt->restrike_point[0] * xyzcen[3 * iel];

          if (rayo < 5.e-4 && posi <= 0.)
            CS_F_(h)->val[iel] = 8.e7;
        }
      }
    }
    else {
      elec_opt->irestrike = 0;
    }
    double somje = 0.;
    for (int iel = 0; iel < ncel; iel++) {
      somje += CS_F_(joulp)->val[iel] * volume[iel];
    }

    cs_parall_sum(1, CS_DOUBLE, &somje);

    if (fabs(somje) > 1.-20)
      coepot = cs_glob_elec_option->couimp * cs_glob_elec_option->pot_diff
              / CS_MAX(somje, cs_math_epzero);

    bft_printf("imposed current %14.5E, Dpot %14.5E, Somje %14.5E\n",
               cs_glob_elec_option->couimp,
               cs_glob_elec_option->pot_diff,
               somje);

    double elcou = 0.;
    for (int ifac = 0; ifac < nfac; ifac++) {
      if (fabs(surfac[ifac][0]) < 1.e-8 && fabs(surfac[ifac][1]) < 1.e-8 &&
               cdgfac[ifac][2] > 0.05e-2 &&     cdgfac[ifac][2] < 0.08e-2) {
        int iel = mesh->i_face_cells[ifac][0];
            elcou += CS_FI_(curre, 2)->val[iel] * surfac[ifac][2];
      }
    }

    cs_parall_sum(1, CS_DOUBLE, &elcou);

    if (fabs(elcou) > 1.e-6)
      elcou = fabs(elcou);
    else
      elcou = 0.;

    if (fabs(elcou) > 1.e20)
      coepoa = cs_glob_elec_option->couimp / elcou;

    coepot = coepoa;

    double dtj = 1.e15;
    double dtjm = dtj;
    double delhsh = 0.;
    double cdtj = 20.;

    for (int iel = 0; iel < ncel; iel++) {
      if (fabs(CS_F_(rho)->val[iel]) > 1.e-20)
        delhsh = CS_F_(joulp)->val[iel] * dt[iel]
               / CS_F_(rho)->val[iel];

      if (fabs(delhsh) > 1.e-20)
        dtjm = CS_F_(h)->val[iel] / delhsh;
      else
        dtjm = dtj;
      dtjm = fabs(dtjm);
      dtj = CS_MIN(dtj, dtjm);
    }
    cs_parall_min(1, CS_DOUBLE, &dtj);

    double cpmx = pow(cdtj * dtj, 0.5);
    coepot = cpmx;

    if (cs_glob_time_step->nt_cur > 3) {
      if (coepoa > 1.05)
        coepot = cpmx;
      else
        coepot = coepoa;
    }

    bft_printf(" Cpmx          = %14.5E\n", cpmx);
    bft_printf(" COEPOA        = %14.5E\n", coepoa);
    bft_printf(" COEPOT        = %14.5E\n", coepot);
    bft_printf(" Dpot rescaled = %14.5E\n",
               cs_glob_elec_option->pot_diff * coepot);

    /* scaling electric fields */
    elec_opt->pot_diff *= coepot;

    /* electric potential (for post treatment) */
    for (int iel = 0; iel < ncel; iel++)
      CS_F_(potr)->val[iel] *= coepot;

    /* current density */
    if (ielarc > 0)
      for (int i = 0; i < 3 ; i++)
        for (int iel = 0; iel < 3 ; iel++)
          CS_FI_(curre, i)->val[iel] *= coepot;

    /* joule effect */
    for (int iel = 0; iel < 3 ; iel++)
      CS_F_(joulp)->val[iel] *= coepot * coepot;
  }
  /*! [electric_scaling] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
