/*============================================================================
 * Functions associated to VOF model
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_cdo_quantities.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_main.h"
#include "cs_convection_diffusion.h"
#include "cs_domain.h"
#include "cs_domain_setup.h"
#include "cs_equation.h"
#include "cs_equation_iterative_solve.h"
#include "cs_face_viscosity.h"
#include "cs_divergence.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_interface.h"
#include "cs_log.h"
#include "cs_physical_constants.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_bad_cells.h"
#include "cs_parall.h"
#include "cs_time_step.h"
#include "cs_rotation.h"
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

  @}

  @}

*/
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*!
  @addtogroup vof
  @{

  VOF model
  - -1: model disabled
  -  0: model enabled

  Note that void fraction variable tracks fluid 2.
*/

int cs_glob_vof_model = -1;

/*
  @}
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static variables
 *============================================================================*/

static cs_vof_parameters_t  _vof_parameters =
{
  .rho1 =  1.e3,
  .rho2 =    1.,
  .mu1  = 1.e-3,
  .mu2  = 1.e-5
};

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_vof_get_pointers(int **ivofmt,
                      double **rho1,
                      double **rho2,
                      double **mu1,
                      double **mu2);

void
cs_f_vof_update_phys_prop(void);

void
cs_f_vof_log_mass_budget(void);

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
 *   ivofmt --> pointer to cs_glob_vof_model
 *   rho1   --> pointer to cs_glob_vof_parameters->rho1
 *   rho2   --> pointer to cs_glob_vof_parameters->rho2
 *   mu1    --> pointer to cs_glob_vof_parameters->mu1
 *   mu2    --> pointer to cs_glob_vof_parameters->mu2
 *----------------------------------------------------------------------------*/

void
cs_f_vof_get_pointers(int    **ivofmt,
                      double **rho1,
                      double **rho2,
                      double **mu1,
                      double **mu2)
{
  *ivofmt = &cs_glob_vof_model;
  *rho1   = &(_vof_parameters.rho1);
  *rho2   = &(_vof_parameters.rho2);
  *mu1    = &(_vof_parameters.mu1);
  *mu2    = &(_vof_parameters.mu2);
}

/*----------------------------------------------------------------------------
 * wrappers to vof functions, intended for use by Fortran wrappers only.
 *----------------------------------------------------------------------------*/

void
cs_f_vof_update_phys_prop(void)
{
  cs_vof_update_phys_prop(cs_glob_domain);
}

void
cs_f_vof_log_mass_budget(void)
{
  cs_vof_log_mass_budget(cs_glob_domain);
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
 * \brief  Compute the mixture density, mixture dynamic viscosity and mixture
 *        mass flux given the volumetric flux, the volume fraction and the
 *        reference density and dynamic viscosity \f$ \rho_l, \mu_l \f$
 *        (liquid), \f$ \rho_v, \mu_v \f$ (gas) as follows:
 *
 * \f[
 * \rho_\celli = \alpha_\celli \rho_v + (1-\alpha_\celli) \rho_l,
 * \f]
 * \f[
 * \mu_\celli = \alpha_\celli \mu_v + (1-\alpha_\celli) \mu_l,
 * \f]
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
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_update_phys_prop(const cs_domain_t *domain)
{
  const cs_mesh_t *m = domain->mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_t *b_face_cells = m->b_face_cells;

  cs_halo_type_t halo_type = m->halo_type;

  cs_real_t *cvar_voidf = CS_F_(void_f)->val;
  cs_real_t *a_voidf = CS_F_(void_f)->bc_coeffs->a;
  cs_real_t *b_voidf = CS_F_(void_f)->bc_coeffs->b;

  cs_real_t *cpro_rom = CS_F_(rho)->val;
  cs_real_t *bpro_rom = CS_F_(rho_b)->val;

  cs_real_t *cpro_viscl = CS_F_(mu)->val;

  const cs_real_t rho1 =_vof_parameters.rho1;
  const cs_real_t rho2 = _vof_parameters.rho2;
  const cs_real_t mu1 = _vof_parameters.mu1;
  const cs_real_t mu2 = _vof_parameters.mu2;

  /*  Update mixture density and viscocity on cells */

# pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t vf = cvar_voidf[c_id];
    cpro_rom[c_id]   = rho2*vf + rho1*(1. - vf);
    cpro_viscl[c_id] =  mu2*vf +  mu1*(1. - vf);
  }

  cs_field_synchronize(CS_F_(rho), halo_type);
  cs_field_synchronize(CS_F_(mu), halo_type);

  /* Update mixture density on boundary faces */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];
    cs_real_t vf = a_voidf[f_id] + b_voidf[f_id]*cvar_voidf[c_id];

    bpro_rom[f_id]   = rho2*vf + rho1*(1. - vf);
  }

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
 * \sum_i\left(
 * |\Omega_i|\dfrac{\alpha_i^n - \alpha_i^{n-1}}{\Delta t} +
 * \sum_{j\in\Face{\celli}}\left(\rho\vect{u}\vect{S}\right)_{ij}^n
 * \right).
 * \f]
 */
/*----------------------------------------------------------------------------*/

void
cs_vof_log_mass_budget(const cs_domain_t *domain)
{
  const cs_mesh_t *m = domain->mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_with_ghosts = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const cs_mesh_quantities_t *mq = domain->mesh_quantities;
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
    BFT_MALLOC(i_massflux_abs, n_i_faces, cs_real_t);
    BFT_MALLOC(b_massflux_abs, n_b_faces, cs_real_t);

    const int *cell_rotor_num = cs_turbomachinery_get_cell_rotor_num();

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      i_massflux_abs[f_id] = i_massflux[f_id];

      cs_lnum_t  c_id_i = i_face_cells[f_id][0];
      cs_lnum_t  c_id_j = i_face_cells[f_id][1];
      int rot_ce_i = cell_rotor_num[c_id_i];
      int rot_ce_j = cell_rotor_num[c_id_j];

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
      int rot_ce_i = cell_rotor_num[c_id];

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

  bft_printf(_("   ** VOF MODEL, MASS BALANCE at iteration %6i: %12.4e\n\n"),
             cs_glob_time_step->nt_cur, glob_m_budget);

  BFT_FREE(divro);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
