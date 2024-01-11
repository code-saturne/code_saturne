/*============================================================================
 * User function. Locally modify a given porosity to take into
 *  account erosion effect (for instance).
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"
#include "cs_ibm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_ibm_modify.c
 *
 * \brief User function. Locally modify a given porosity to take into
 *         account erosion effect (for instance).
 */
/*----------------------------------------------------------------------------*/

static void
_smoothe(const cs_mesh_t            *mesh,
         const cs_mesh_quantities_t *mesh_quantities,
         cs_real_t                   val[],
         int                         nloop)
{
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells     = mesh->n_cells;
  cs_lnum_t n_i_faces   = mesh->n_i_faces;

  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  cs_real_t *cell_vol  = mesh_quantities->cell_vol;

  cs_real_t *den, *val2;
  BFT_MALLOC(val2, n_cells_ext, cs_real_t);
  BFT_MALLOC(den, n_cells_ext, cs_real_t);

  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, val);

  for (int iloop = 1; iloop <= nloop; iloop++) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      val2[c_id] = val[c_id] * cell_vol[c_id];
      den[c_id] = cell_vol[c_id];
    }

    for (cs_lnum_t c_id = n_cells; c_id < n_cells_ext; c_id++) {
      val2[c_id] = 0.;
      den[c_id] = 0.;
    }

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id0 = i_face_cells[f_id][0];
      cs_lnum_t c_id1 = i_face_cells[f_id][1];

      val2[c_id0] += val[c_id1] * cell_vol[c_id1] * CS_F_(poro)->val[c_id1];
      val2[c_id1] += val[c_id0] * cell_vol[c_id0] * CS_F_(poro)->val[c_id0];
      den[c_id0] += cell_vol[c_id1] * CS_F_(poro)->val[c_id1];
      den[c_id1] += cell_vol[c_id0] * CS_F_(poro)->val[c_id0];
    }

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      val[c_id] = val2[c_id] / den[c_id];

    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, val);
  }

  BFT_FREE(val2);
  BFT_FREE(den);
}

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function. Locally modify a given porosity to take into
 *         account erosion effect (for instance)
 *
 *  This function is called for each time step.
 *  Porosity will be modified if
 *  cs_ibm->porosity_user_source_term_modification = true
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 *
 *----------------------------------------------------------------------------*/

void
cs_user_ibm_modify(const cs_mesh_t            *mesh,
                   const cs_mesh_quantities_t *mesh_quantities)
{

  /*!< [loc_var_def_init] */
  cs_lnum_t n_cells     = mesh->n_cells;
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_i_faces   = mesh->n_i_faces;
  cs_lnum_t n_b_faces   = mesh->n_b_faces;

  cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;

  cs_real_t *cell_vol = mesh_quantities->cell_vol;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)mesh_quantities->cell_cen;
  const cs_real_3_t *i_face_normal
    = (const cs_real_3_t *)mesh_quantities->i_face_normal;
  const cs_real_3_t *c_w_face_normal
    = (const cs_real_3_t *)mesh_quantities->c_w_face_normal;
  const cs_real_t *c_w_face_surf
    = (const cs_real_t *)mesh_quantities->c_w_face_surf;
  /*!< [loc_var_def_init] */

  /*!< [example_1] */

  /* Initialize porosity to previous time-step porosity in
   * a defined region (here 11 < x < 16) */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (cell_cen[c_id][0] > 10. && cell_cen[c_id][0] < 16.)
      CS_F_(poro)->val[c_id] = CS_F_(poro)->val_pre[c_id];

  cs_real_t *source_term, *por_init;
  cs_real_3_t *convective_term;

  BFT_MALLOC(source_term, n_cells_ext, cs_real_t);
  BFT_MALLOC(por_init, n_cells_ext, cs_real_t);
  BFT_MALLOC(convective_term, n_cells_ext, cs_real_3_t);

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    for (int i = 0; i < 3; i++)
      convective_term[c_id][i] = 0.;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    cs_real_t flux[3];
    for (int idim = 0; idim < 3; idim++)
      flux[idim] = 0.;

    cs_real_t roij = 0.5 * (CS_F_(rho)->val[c_id0] + CS_F_(rho)->val[c_id1]);
    cs_real_t alpi = 1.;
    cs_real_t alpj = 1.;

    cs_real_t fluij = 0.;
    for (int idim = 0; idim < 3; idim++)
      fluij += 0.5 * (  alpi * CS_F_(vel)->val_pre[3 * c_id0 + idim]
                      + alpj * CS_F_(vel)->val_pre[3 * c_id1 + idim])
                   * i_face_normal[f_id][idim];
    fluij *= roij;

    cs_real_t  alpij = 1.;

    cs_real_t uij[3];
    for (int idim = 0; idim < 3; idim++)
      uij[idim] = CS_F_(vel)->val_pre[3 * c_id0 + idim];

    if (fluij < 0.)
      for (int idim = 0; idim < 3; idim++)
        uij[idim] = CS_F_(vel)->val_pre[3 * c_id1 + idim];

    for (int idim = 0; idim < 3; idim++)
      flux[idim] += alpij * fluij * uij[idim];

    for (int idim = 0; idim < 3; idim++) {
      convective_term[c_id0][idim] += flux[idim];
      convective_term[c_id1][idim] -= flux[idim];
    }
  }

  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  cs_real_t *restrict b_massflux =
    cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf))->val;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    cs_real_3_t flux;
    for (int idim = 0; idim < 3; idim++)
      flux[idim] = 0.;

    cs_real_t fluij = b_massflux[f_id];
    cs_real_t alpij = 1.;

    for (int idim = 0; idim < 3; idim++)
      flux[idim] += alpij * fluij * CS_F_(vel)->val_pre[3 * c_id + idim];

    for (int idim = 0; idim < 3; idim++)
      convective_term[c_id][idim] += flux[idim];
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_3_t np;
    for (int idim = 0; idim < 3; idim++)
      np[idim] = c_w_face_normal[c_id][idim];

    cs_real_t nn = c_w_face_surf[c_id];
    nn = cs_math_fmax(nn, 1.e-20);

    for (int idim = 0; idim < 3; idim++)
      np[idim] /= nn;

    cs_real_t cc = cs_math_3_dot_product(convective_term[c_id], np);
    cc = cs_math_fabs(cc);
    source_term[c_id] = cc;
  }

  int nloop = 2;
  _smoothe(mesh, mesh_quantities, source_term, nloop);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (cell_cen[c_id][0] < 11 || cell_cen[c_id][0] > 16)
      source_term[c_id] = 0.;

  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, source_term);

  /* Mass transfer term in kg/m^2/s */
  /* Solid density */
  cs_real_t ros = 2500.;

  /* Multiplier coeff. */
  cs_real_t coeff = 1000.;

  cs_real_t dt = cs_glob_time_step->dt[0];

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    por_init[c_id] = CS_F_(poro)->val[c_id];

  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, por_init);

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    cs_real_t pori = por_init[c_id0];
    cs_real_t porj = por_init[c_id1];
    if (pori > 0.9999 && porj < 0.0001) {
      cs_real_t surfi = c_w_face_surf[c_id0];
      cs_real_t ci = source_term[c_id0];
      cs_real_t gama = coeff * ci / cs_math_fmax(porj, 0.5);
      CS_F_(poro)->val[c_id1] += gama * surfi / (ros * cell_vol[c_id1]) * dt;
    } else if (porj > 0.9999 && pori < 0.0001) {
      cs_real_t surfj = c_w_face_surf[c_id1];
      cs_real_t cj = source_term[c_id1];
      cs_real_t gama = coeff * cj / cs_math_fmax(pori, 0.5);
      CS_F_(poro)->val[c_id0] += gama * surfj / (ros * cell_vol[c_id0]) * dt;
    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t surf = c_w_face_surf[c_id];
    cs_real_t cc = source_term[c_id];
    cs_real_t por = por_init[c_id];
    cs_real_t gama = coeff * cc / cs_math_fmax(por, 0.5);
    CS_F_(poro)->val[c_id] += gama * surf / (ros * cell_vol[c_id]) * dt;
  }

  /* Clipping of the porosity */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t porosi = CS_F_(poro)->val[c_id];
    if (porosi > 1.) {
      CS_F_(poro)->val[c_id] = 1.;
    } else if (porosi < 1.e-5) {
      CS_F_(poro)->val[c_id] = 0.;
    }
  }

  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, CS_F_(poro)->val);

  BFT_FREE(por_init);
  BFT_FREE(convective_term);
  BFT_FREE(source_term);

  /*!< [example_1] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
