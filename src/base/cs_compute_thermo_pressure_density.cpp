/*============================================================================
 * Thermodynamic pressure and density.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_boundary_conditions.h"
#include "cs_domain.h"
#include "cs_equation.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_velocity_pressure.h"
#include "cs_volume_mass_injection.h"
#include "cs_wall_condensation.h"
#include "cs_wall_condensation_1d_thermal.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_compute_thermo_pressure_density.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_compute_thermo_pressure_density.c

  Compute the thermodynamic pressure and dansity

*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the thermodynamic pressure for a perfect gas.
 *
 * \param[in]  n_cells      number of cellules
 * \param[in]  n_b_faces    number of boundary faces
 * \param[in]  verbosity    verbosity level
 * \param[in]  crom         density in cellules
 * \param[in]  cell_vol     volume of cellules
 * \param[out] new_pther    thermodynamic pressure p_ther^(n+1)
 */
/*----------------------------------------------------------------------------*/

static void
_compute_thermodynamic_pressure_perfect_gas(const cs_lnum_t n_cells,
                                            const cs_lnum_t n_b_faces,
                                            const int       verbosity,
                                            const cs_real_t *crom,
                                            const cs_real_t *cell_vol,
                                            cs_real_t       *new_pther)
{
  const cs_real_t *b_face_surf = cs_glob_mesh_quantities->b_face_surf;

  const int *bc_type = cs_glob_bc_type;

  const cs_lnum_t ncmast = cs_glob_wall_condensation->ncmast;
  const cs_lnum_t nfbpcd = cs_glob_wall_condensation->nfbpcd;
  const cs_lnum_t *ifbpcd = cs_glob_wall_condensation->ifbpcd;

  const int var_id_key = cs_field_key_id("variable_id");
  const int ipr = cs_field_get_key_int(CS_F_(p), var_id_key) - 1;
  const cs_real_t *spcond = cs_glob_wall_condensation->spcond + ipr*nfbpcd;
  const cs_real_t *svcond = cs_glob_wall_condensation->svcond + ipr*ncmast;

  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmab =  cs_field_get_key_int(CS_F_(p), kbmasf);
  const cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  cs_real_t _new_pther = 0.0;

  cs_real_t *cromo = CS_F_(rho)->val;
  if (CS_F_(rho)->n_time_vals > 1)
    cromo = CS_F_(rho)->val_pre;

  cs_real_t debin  = 0.0;
  cs_real_t debout = 0.0;
  cs_real_t debtot = 0.0;

  // Computation of mass flux imposed on the boundary faces

# pragma omp parallel for reduction(-:debin, debout) if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    if (   (bc_type[face_id] == CS_INLET)
        || (bc_type[face_id] == CS_CONVECTIVE_INLET)    )
      // the inlet mass flux integrated in space
      debin -= bmasfl[face_id];
    else if (   (bc_type[face_id] == CS_OUTLET)
             || (bc_type[face_id] == CS_FREE_INLET)    )
       debout -= bmasfl[face_id];
  }

  debtot = debin + debout;

  // Computation of the inlet mass flux imposed on the cells volume

  cs_lnum_t ncetsm = 0;
  const cs_lnum_t *icetsm = nullptr;
  cs_real_t *smacel;
  cs_volume_mass_injection_get_arrays(CS_F_(p),
                                      &ncetsm,
                                      &icetsm,
                                      nullptr,
                                      &smacel,
                                      nullptr);

  for (cs_lnum_t ii = 0; ii < ncetsm; ii++) {
    const cs_lnum_t c_id = icetsm[ii];
    debtot += smacel[ii] * cell_vol[c_id];
  }

  /* Flow rate computation associated to the condensation phenomena
     -------------------------------------------------------------- */

  /* Sink source term associated to
     the surface condensation modelling */
# pragma omp parallel for reduction(+:debtot) if (nfbpcd > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < nfbpcd; ii++) {
    const cs_lnum_t face_id = ifbpcd[ii];
    debtot += b_face_surf[face_id] * spcond[ii];
  }

  /* Sink source term associated to
     the metal structures condensation modelling */
  if (cs_glob_wall_condensation->icondv == 0) {
    const cs_lnum_t *ltmast = cs_glob_wall_condensation->ltmast;
    const cs_lnum_t *izmast = cs_glob_wall_condensation->izmast;
    const cs_real_t *volume_surf = cs_glob_wall_cond_0d_thermal->volume_surf;
    const cs_real_t *volume_measure
      = cs_glob_wall_cond_0d_thermal->volume_measure;

    for (cs_lnum_t ii = 0; ii < ncmast; ii++) {
      const cs_lnum_t c_id = ltmast[ii];
      const cs_lnum_t v_id = izmast[ii];
      const cs_real_t surfbm
        = volume_surf[v_id]*cell_vol[c_id]/volume_measure[v_id];
      debtot += surfbm*svcond[ii];
    }
  }

  cs_parall_sum(1, CS_REAL_TYPE, &debtot);

  /* Global leak
     ----------- */
  const cs_fluid_properties_t *fp = cs_glob_fluid_properties;

  cs_real_t rho = fp->roref;
  cs_real_t dp = fp->pther - fp->p0;
  if (dp > 0.0)
    rho = fp->ro0;

  // for sign(1.0, dp) in fortran
  if (fabs(dp) <= 0.0)
    dp = 1.0;
  else
    dp = dp/fabs(dp);

  debtot -= dp * fp->sleak * sqrt(2.0*rho/fp->kleak*fabs(dp));

  // for the first time step : rho^(n-1) = rho^(n)
  if ((cs_restart_present() == 0) && (cs_glob_time_step->nt_cur == 1))
    cs_array_real_copy(n_cells, crom, cromo);

  // initialization
  cs_real_t romoy  = 0.0;
  cs_real_t roamoy = 0.0;
# pragma omp parallel reduction(+:romoy, roamoy) if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      romoy += crom[c_id]*cell_vol[c_id];
      roamoy += cromo[c_id]*cell_vol[c_id];
  }

  cs_parall_sum(1, CS_REAL_TYPE, &romoy);
  cs_parall_sum(1, CS_REAL_TYPE, &roamoy);

  // Compute the thermodynamic pressure p_ther^(n+1)
  const cs_real_t dt = CS_F_(dt)->val[0];
  _new_pther = fp->pther*(roamoy/romoy + dt*debtot/romoy);

  // pthermodynamic pressure clipping in case of user venting
  if (fp->pthermax > 0.0)
    _new_pther = cs_math_fmin(_new_pther, fp->pthermax);

  if (   (verbosity >= 1)
      && (cs_glob_velocity_pressure_model->idilat == 3)   )
    if (cs_log_default_is_active())
      bft_printf("** Perfect gas computation of average td_pressure:\n"
                 "   -----------------------------------------------\n"
                 "-------------------------------------------------------\n"
                 "time       rho(n-1)/rho(n)     dt.debtot/ro(n)"
                 "-------------------------------------------------------\n"
                 "%10.16lf       %10.16lf     %10.16lf\n",cs_glob_time_step->t_cur,
                 roamoy/romoy, (dt*debtot)/romoy);

  *new_pther = _new_pther;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the density \f$ \rho^{n+1}\f$ with the
 *        \f$ \rho^{n-\frac{1}{2}} \f$ density with the state law
 *        and a thermodynamic pressure \f$ p_{ther}^{n+1} \f$
 *        estimated from the integral over the total
 *        fluid domain of the mass conservation equation.
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_thermo_pressure_density(void)
{
  const cs_mesh_t *m = cs_glob_mesh;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *brom = CS_F_(rho_b)->val;

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(CS_F_(p));

  /* Flow rate computation for the inlet and oulet conditions
     -------------------------------------------------------- */

  /* Update the thermodynamic pressure
   * for the previous time step
   * ---------------------------------- */

  fp->pthera = fp->pther;
  cs_real_t new_pther = 0.0;
  _compute_thermodynamic_pressure_perfect_gas(n_cells,
                                              n_b_faces,
                                              eqp->verbosity,
                                              crom,
                                              cell_vol,
                                              &new_pther);

  cs_user_physical_properties_td_pressure(&new_pther);

  fp->pther = new_pther;

  /* Thermodynamic pressure and density computation
     ---------------------------------------------- */
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    crom[c_id] = fp->pther/fp->pthera*crom[c_id];

  if (m->halo != nullptr)
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, crom);

  /* Update the density at the boundary face
   * with cell value for severe accident low-Mach algorithm */
  if (cs_glob_velocity_pressure_model->idilat == 3)
#   pragma omp parallel for if (n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      const cs_lnum_t c_id = b_face_cells[face_id];
      brom[face_id] = crom[c_id];
    }
  // else with boundary values
  else
#   pragma omp parallel for if (n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
      brom[face_id] = fp->pther/fp->pthera*brom[face_id];

  /* Change the reference variable rho0
     ---------------------------------- */
  cs_real_t ro0moy = 0.0;
# pragma omp parallel for reduction(+:ro0moy) if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    ro0moy += crom[c_id]*cell_vol[c_id];

  cs_parall_sum(1, CS_REAL_TYPE, &ro0moy);

  fp->ro0 = ro0moy/cs_glob_mesh_quantities->tot_vol;

  if (eqp->verbosity >= 1)
    if (cs_log_default_is_active())
      bft_printf("** Thermodynamic pressure computation:\n"
                 "   -----------------------------------------------\n"
                 "-------------------------------------------------------\n"
                 "time       pther(n+1)     pther(n)     Dpther/Dt     ro0"
                 "-------------------------------------------------------\n"
                 "%10.16lf       %10.16lf     %10.16lf     %10.16lf     %10.16lf\n",
                 cs_glob_time_step->t_cur, fp->pther, fp->pthera,
                 (fp->pther-fp->pthera)/CS_F_(dt)->val[0], fp->ro0);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
