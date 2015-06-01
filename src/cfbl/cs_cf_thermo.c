/*============================================================================
 * Thermodynamic laws for the compressible module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_thermal_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cf_thermo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cf_thermo.c
        Define thermodynamic laws for the compressible flow module.
        Only the perfect gas law is available for now. The molar mass has to be
        provided either in the GUI or in cs_user_parameters.f90.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set variability of isobaric specific heat and isochoric specific heat
 * according to the chosen thermodynamic law.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_set_thermo_options(void)
{
  if (cs_glob_fluid_properties->ieos == 1) {
    /* Calculation options: constant Cp and Cv (perfect gas)
       specific heat Cv0 is calculated in a subsequent section (from Cp0) */
    cs_glob_fluid_properties->icp = 0;
    cs_glob_fluid_properties->icv = 0;
  }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve molar mass.
 *
 * \param[out]   xmasml    molar mass
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_get_molar_mass(cs_real_t *xmasml)
{
  if (cs_glob_fluid_properties->ieos == 1) {
    /*  Molar mass of the gas (kg/mol) */
    *xmasml = cs_glob_fluid_properties->xmasmr;
  }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute \f$\gamma\f$.
 *
 * \param[out]   gamma    ratio of specific heat
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_gamma(cs_real_t *gamma)
{
  /*  Local variables */
  cs_real_t r_pg = cs_glob_physical_constants->r;
  cs_real_t cp0 = cs_glob_fluid_properties->cp0;
  cs_real_t gamagp = 0;
  cs_real_t xmasml;
  cs_cf_get_molar_mass(&xmasml);

  /*  Gamagp is supposed to be superior or equal to 1.
      It is computed at each call, even if this may seem costly,
      to be coherent with the "constant gamma" case for which this
      constant is not saved. A ''save'' instruction and a test would
      be sufficient to avoid computing gamagp at each call if necessary. */
  if (cs_glob_fluid_properties->ieos == 1) {
    gamagp = 1. + r_pg/(xmasml*cp0-r_pg);

    if (gamagp < 1.)
      bft_error(__FILE__, __LINE__, 0,
                _("Error in thermodynamics computations for "
                  "compressible flows:\n"
                  "the thermodynamic law is perfect gas with constant Gamma\n"
                  "with Gamma=%12.4e <= 1.\n"
                  "Gamma must be a real number greater or equal to 1."),
                gamagp);
  }

  *gamma = gamagp;
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize density, total energy and isochoric specific heat
 * according to the chosen thermodynamic law using the default parameters.
 *
 * \param[in]     isuite    indicator of resumed calculation
 * \param[in]     l_size    local size
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_default_init(int isuite, cs_lnum_t l_size)
{
  /* Local variables */
  cs_real_t xmasml;
  cs_real_t  r_pg  = cs_glob_physical_constants->r;
  cs_real_t  p0  = cs_glob_fluid_properties->p0;
  cs_real_t  t0  = cs_glob_fluid_properties->t0;
  cs_real_t  cp0 = cs_glob_fluid_properties->cp0;
  cs_real_t *cv0 = &cs_glob_fluid_properties->cv0;

  /* Default initializations
     t0 is positive (this assumption has been checked in verini) */
  cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *cvar_en = CS_F_(energy)->val;

  cs_cf_get_molar_mass(&xmasml);

  if (cs_glob_fluid_properties->ieos == 1) {
    *cv0 = cp0 - r_pg/xmasml;

    if (isuite == 0) {
      for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++) {
        crom[cell_id] = p0 * xmasml/(r_pg*t0);
        cvar_en[cell_id] = *cv0 * t0;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

// TODO: the check function should be generalized (pass the name as argument).

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the positivity of the pressure.
 *
 * \param[in]     pres    array of pressure values
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_check_pressure(cs_real_t *pres,
                     cs_lnum_t l_size)
{
  /* Local variables */
  cs_gnum_t ierr;

  /* If the pressure is lower or equal to zero, stop the computation.
     Indeed, if this is the case, the thermodynamic computations will most
     probably fail. This call is done at the end of the density calculation */
  ierr = 0;
  for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++)
    if (pres[cell_id] <= 0.)
      ierr = ierr + 1;

  if (cs_glob_rank_id >= 0) cs_parall_counter(&ierr, 1);

  if (ierr > 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error in thermodynamics computations for compressible flows\n"
                "(perfect gas with constant Gamma):\n"
                "Negative values of the pressure were encountered in %lu"
                " cells.\n"), ierr);
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the positivity of the internal energy.

 * \param[in]     ener    array of total energy values
 * \param[in]     l_size  l_size of the array
 * \param[in]     vel     array of velocity values
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_check_internal_energy(cs_real_t   *ener,
                            cs_lnum_t    l_size,
                            cs_real_3_t *vel)
{
  /* Local variables */
  cs_gnum_t ierr;
  cs_real_t enint;

  /* If the internal energy <= zero: stop the computation.
     Indeed, if this is the case, the thermodynamic computations will
     most probably fail. */
  ierr = 0;
  for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++) {
    cs_real_t v2 = cs_math_3_square_norm(vel[cell_id]);
    enint = ener[cell_id] - 0.5*v2;

    if (enint <= 0.)
      ierr++;
  }

  if (cs_glob_rank_id >= 0)
    cs_parall_counter(&ierr, 1);

  if (ierr > 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error in thermodynamics computations for compressible flows\n"
                "(perfect gas with constant Gamma):\n"
                "Negative values of the internal energy were encountered in %lu"
                " cells.\n"), ierr);
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the positivity of the density given by the user.
 *
 * \param[in]     dens    array of density values
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_check_density(cs_real_t *dens,
                    cs_lnum_t l_size)
{
  /* Local variables */
  cs_gnum_t ierr;

  /* Verification of the values of the density
     Stop if a negative value is detected (since the density has been
     provided by the user, one potential cause is a wrong user
     initialization). */
  ierr = 0;
  for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++)
    if (dens[cell_id] <= 0.)
      ierr = ierr + 1;

  if (cs_glob_rank_id >= 0) cs_parall_counter(&ierr, 1);

  if (ierr > 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error in thermodynamics computations for compressible flows\n"
                "(perfect gas with constant Gamma):\n"
                "Negative values of the density were encountered in %lu"
                " cells.\n"), ierr);
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check strict positivity of temperature (Celsius) given by the user.

 * \param[in]     temp    array of temperature values
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_check_temperature(cs_real_t *temp,
                     cs_lnum_t l_size)
{
  /* Local variables */
  cs_gnum_t ierr;

  /* Verification of the values of the temperature
     Stop if a negative value is detected (since the temperature has been
     provided by the user, one potential cause is a wrong user
     initialization). */
  ierr = 0;
  for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++)
    if (temp[cell_id] <= 0.)
      ierr++;

  if (cs_glob_rank_id >= 0) cs_parall_counter(&ierr, 1);

  if (ierr > 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error in thermodynamics computations for compressible flows\n"
                "(perfect gas with constant Gamma):\n"
                "Negative values of the temperature were encountered in %lu"
                " cells.\n"), ierr);
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute temperature and total energy from density and pressure.
 *
 * \param[in]     pres    array of pressure values
 * \param[in]     dens    array of density values
 * \param[out]    temp    array of temperature values
 * \param[out]    ener    array of total energy values
 * \param[in]     vel     array of velocity component values
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_te_from_dp(cs_real_t   *pres,
                        cs_real_t   *dens,
                        cs_real_t   *temp,
                        cs_real_t   *ener,
                        cs_real_3_t *vel,
                        cs_lnum_t    l_size)
{
  /* local variables */
  cs_real_t xmasml;
  cs_real_t r_pg = cs_glob_physical_constants->r;

  /*  calculation of temperature and energy from pressure and density */
  cs_cf_get_molar_mass(&xmasml);

  if (cs_glob_fluid_properties->ieos == 1) {
    cs_real_t cv0 = cs_glob_fluid_properties->cv0;

    for (cs_lnum_t ii = 0; ii < l_size; ii++) {
      /*  temperature */
      temp[ii] = xmasml * pres[ii] / (r_pg*dens[ii]);
      /*  total energy */
      cs_real_t v2 = cs_math_3_square_norm(vel[ii]);
      ener[ii] =  cv0*temp[ii] + 0.5*v2;
    }
  }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute density and total energy from pressure and temperature
 *
 * \param[in]     pres    array of pressure values
 * \param[in]     temp    array of temperature values
 * \param[out]    dens    array of density values
 * \param[out]    ener    array of total energy values
 * \param[in]     vel     array of velocity component values
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_de_from_pt(cs_real_t   *pres,
                        cs_real_t   *temp,
                        cs_real_t   *dens,
                        cs_real_t   *ener,
                        cs_real_3_t *vel,
                        cs_lnum_t    l_size)
{
  /* Local variables */
  cs_real_t xmasml;
  cs_real_t r_pg = cs_glob_physical_constants->r;

  /*  Calculation of density and energy from pressure and temperature */
  cs_cf_get_molar_mass(&xmasml);

  if (cs_glob_fluid_properties->ieos == 1) {
    cs_real_t cv0 = cs_glob_fluid_properties->cv0;

    for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++) {
      /*  Temperature */
      dens[cell_id] = xmasml * pres[cell_id] / (r_pg*temp[cell_id]);
      /*  Total energy */
      cs_real_t v2 = cs_math_3_square_norm(vel[cell_id]);
      ener[cell_id] =  cv0*temp[cell_id] + 0.5*v2;
    }
  }
}

/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute density and temperature from pressure and total energy;
 *        (interleaved version).

 * \param[in]     pres    array of pressure values
 * \param[in]     ener    array of total energy values
 * \param[out]    dens    array of density values
 * \param[out]    temp    array of temperature values
 * \param[in]     vel     array of velocity component values
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_dt_from_pe(cs_real_t   *pres,
                        cs_real_t   *ener,
                        cs_real_t   *dens,
                        cs_real_t   *temp,
                        cs_real_3_t *vel,
                        cs_lnum_t    l_size)
{
  /* Local variables */
  cs_real_t enint, gamagp, xmasml;
  cs_real_t r_pg = cs_glob_physical_constants->r;

  /*  Calculation of density and temperature from pressure and energy */
  cs_cf_get_molar_mass(&xmasml);

  if (cs_glob_fluid_properties->ieos == 1) {
    cs_cf_thermo_gamma(&gamagp);

    for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++) {
      /*  Internal energy (to avoid the need to divide by the temperature
          to compute density) */
     cs_real_t v2 = cs_math_3_square_norm(vel[cell_id]);
      enint =  ener[cell_id] - 0.5*v2;

      /*  Density */
      dens[cell_id] = pres[cell_id] / ( (gamagp-1.) * enint );
      /*  Temperature */
      temp[cell_id] = xmasml * (gamagp-1.) * enint / r_pg;
    }
  }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute pressure and total energy from density and temperature

 * \param[in]     dens    array of density values
 * \param[in]     temp    array of temperature values
 * \param[out]    pres    array of pressure values
 * \param[out]    ener    array of total energy values
 * \param[in]     vel     array of velocity component values
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_pe_from_dt(cs_real_t   *dens,
                        cs_real_t   *temp,
                        cs_real_t   *pres,
                        cs_real_t   *ener,
                        cs_real_3_t *vel,
                        cs_lnum_t    l_size)
{
  /* Local variables */
  cs_real_t xmasml;
  cs_real_t r_pg = cs_glob_physical_constants->r;

  /* Calculation of pressure and energy from density and temperature */
  cs_cf_get_molar_mass(&xmasml);

  if (cs_glob_fluid_properties->ieos == 1) {
    cs_real_t cv0 = cs_glob_fluid_properties->cv0;

    for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++) {
      /*  Pressure */
      pres[cell_id] = dens[cell_id]*temp[cell_id]*r_pg/xmasml;
      /*  Total energy */
      cs_real_t v2 = cs_math_3_square_norm(vel[cell_id]);
      ener[cell_id] =  cv0*temp[cell_id] + 0.5*v2;
    }
  }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute pressure and temperature from density and total energy.

 * \param[in]     dens    array of density values
 * \param[in]     ener    array of total energy values
 * \param[out]    pres    array of pressure values
 * \param[out]    temp    array of temperature values
 * \param[in]     vel     array of velocity component values
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_pt_from_de(cs_real_t   *dens,
                        cs_real_t   *ener,
                        cs_real_t   *pres,
                        cs_real_t   *temp,
                        cs_real_3_t *vel,
                        cs_lnum_t    l_size)
{
  /*  Local variables */
  cs_real_t enint, gamagp, xmasml;
  cs_real_t r_pg = cs_glob_physical_constants->r;

  /*  Calculation of pressure and temperature from density and energy */
  cs_cf_get_molar_mass(&xmasml);

  if (cs_glob_fluid_properties->ieos == 1) {
    cs_cf_thermo_gamma(&gamagp);

    for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++) {
      /*  Internal energy (to avoid the need to divide by the temperature
          to compute density) */
      cs_real_t v2 = cs_math_3_square_norm(vel[cell_id]);
      enint =  ener[cell_id] - 0.5*v2;

      /*  Pressure */
      pres[cell_id] = (gamagp-1.) * dens[cell_id] * enint;
      /*  Temperature */
      temp[cell_id] = xmasml * (gamagp-1.) * enint / r_pg;
    }
  }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute square of sound velocity:
 *
 * \f[c^2  = \left(\frac{\partial p}{\partial \rho}\right)_s\f];
 *
 * for perfect gas, this expression simply writes:
 *
 * \f[c^2  = \gamma \frac{p}{\rho}\f]
 *
 * \param[in]    pres    array of pressure values
 * \param[in]    dens    array of density values
 * \param[out]   c2      array of the values of the square of sound velocity
 * \param[in]    l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_c_square(cs_real_t *pres,
                      cs_real_t *dens,
                      cs_real_t *c2,
                      cs_lnum_t l_size)
{
  /*  Local variables */
  cs_real_t gamagp;

  if (cs_glob_fluid_properties->ieos == 1) {
    cs_cf_thermo_gamma(&gamagp);

    for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++)
      c2[cell_id] = gamagp * pres[cell_id] / dens[cell_id];
  }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the thermal expansion coefficient:
 *
 * \f[ \beta = \left(\frac{\partial p}{\partial s}\right)_\rho \f]
 *
 * for a perfect gas, the explicit formula is:
 *
 * \f[ \beta = \rho^\gamma \f]
 *
 * \param[in]    dens    array of density values
 * \param[out]   beta    array of beta values
 * \param[in]    l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_beta(cs_real_t *dens,
                  cs_real_t *beta,
                  cs_lnum_t  l_size)
{
  /*  Local variables */
  cs_real_t gamagp;

  if (cs_glob_fluid_properties->ieos == 1) {
    cs_cf_thermo_gamma(&gamagp);

    for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++)
      beta[cell_id] = pow(dens[cell_id],gamagp);
  }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the isochoric specific heat:
 *
 * \f[C_v = \left(\frac{\partial e}{\partial T}\right)_\rho\f]
 *
 * \param[in]     cp      array of isobaric specific heat values
 * \param[out]    cv      array of isochoric specific heat values
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_cv(cs_real_t *cp,
                cs_real_t *cv,
                cs_lnum_t  l_size)
{
  /*  Constant quantity if perfect gas chosen,
      nothing to be done. */
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute entropy from pressure and density:
 *
 * \f[s = \frac{p}{\rho^\gamma}\f]
 *
 * \param[in]     dens    array of density values
 * \param[in]     pres    array of pressure values
 * \param[out]    entr    array of total energy values
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_s_from_dp(cs_real_t *dens,
                       cs_real_t *pres,
                       cs_real_t *entr,
                       cs_lnum_t  l_size)
{
  /*  Local variables */
  cs_real_t gamagp;

  if (cs_glob_fluid_properties->ieos == 1) {
    cs_cf_check_density(dens, l_size);

    cs_cf_thermo_gamma(&gamagp);

    for (cs_lnum_t cell_id = 0; cell_id < l_size; cell_id++)
      entr[cell_id] = pres[cell_id] / pow(dens[cell_id],gamagp);
  }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute wall boundary condition values.
 *
 * \param[in,out] wbfb    output work array
 * \param[in]     face_id boundary face index
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_wall_bc(cs_real_t *wbfb,
                     cs_lnum_t  face_id)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;

  /*  Map field arrays */
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(u)->val;
  cs_real_t *cvar_pr = CS_F_(p)->val;
  cs_real_t *crom = CS_F_(rho)->val;

  cs_real_t gamagp;

  if (cs_glob_fluid_properties->ieos == 1) {

    /*  Calculation of the boundary conditions on the wall face face_id */

    cs_cf_thermo_gamma(&gamagp);

    cs_lnum_t cell_id = b_face_cells[face_id];

    /*  Calculation of the Mach number at the boundary face, using the
        cell center velocity projected on the vector normal to the boundary */

    cs_real_t uni = cs_math_3_dot_product(vel[cell_id],b_face_normal[face_id])
                   /b_face_surf[face_id];
    cs_real_t xmach = uni / sqrt( gamagp * cvar_pr[cell_id] / crom[cell_id] );

    /*  Pressure */

    /* A Neumann boundary condition is used. This does not allow to use
       the Rusanov scheme, but some stabilization effect is expected.
       A test based on the value of coefb at the previous time step
       is implemented to avoid oscillating between a rarefaction
       situation and a shock configuration from one time step to the
       next. */

    /*  Rarefaction !FIXME with the new cofaf cofbf */
    if (xmach < 0. && wbfb[face_id] <= 1.) {

      if (xmach > 2./(1.-gamagp)) {

        wbfb[face_id] = pow(1. + (gamagp-1.)/2. * xmach, 2.*gamagp/(gamagp-1.));
      }
      else {
        /* In case the rarefaction is too strong, a zero Dirichlet value
           is used for pressure (the value of wbfb is used here as an
           indicator) */
        wbfb[face_id] = 1.e30;
      }

    }
    /*  Shock */
    else if (xmach > 0. && wbfb[face_id] >= 1.) {

      wbfb[face_id] = 1. + gamagp*xmach
                          *((gamagp+1.)/4.*xmach
                            + sqrt(1. + pow(gamagp+1.,2/16.*xmach*xmach)));
    }
    /*  Oscillation between rarefaction and shock or zero Mach number */
    else {
      wbfb[face_id] = 1.;
    }

  }

}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute subsonic outlet boundary conditions values.

 * \param[in,out] bc_en   total energy values at boundary faces
 * \param[in,out] bc_pr   pressure values at boundary faces
 * \param[in,out] bc_vel  velocity values at boundary faces
 * \param[in]     face_id    boundary face index
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_subsonic_outlet_bc(cs_real_t   *bc_en,
                                cs_real_t   *bc_pr,
                                cs_real_3_t *bc_vel,
                                cs_lnum_t    face_id)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;

  cs_real_t gamagp;
  cs_real_t roi, ro1, pri, uni, un1, uns;
  cs_real_t ci, c1, mi, a, b, sigma1, pinf;

  /*  Map field arrays */
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(u)->val;
  cs_real_t *cvar_pr = CS_F_(p)->val;
  cs_real_t *cvar_en = CS_F_(energy)->val;
  cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *brom = CS_F_(rho_b)->val;

  /*  Calculation of the boundary conditions on the subsonic outlet face face_id */

  if (cs_glob_fluid_properties->ieos == 1) {

    cs_cf_thermo_gamma(&gamagp);

    cs_lnum_t cell_id = b_face_cells[face_id];

    pinf = bc_pr[face_id];
    pri  = cvar_pr[cell_id];
    roi  = crom[cell_id];

    ci   = sqrt(gamagp * pri / roi);
    uni = cs_math_3_dot_product(vel[cell_id],b_face_normal[face_id])
         /b_face_surf[face_id];

    /*  Rarefaction case */
    if (pinf <= pri) {

      /* Computation of the velocity in state 1 using Riemann invariants
         of the 1-rarefaction */
      a = 2 * ci / (gamagp - 1.)
          * (1. -  pow(pinf / pri, (gamagp - 1.) /(2. * gamagp)) );
      un1 = uni + a;

      /* Computation of the density in state 1 using Rieman invariants
         of the 1-rarefaction */
      ro1 = roi * pow(pinf/pri, 1./gamagp);

      /* Subsonic inlet - state 2 should be imposed but too few information
         is available to compute it
         for want of anything better, state 1 is imposed */
      if (un1 < 0.) {

        /*  Density */
        brom[face_id] = ro1;
        /*  Velocity */
        bc_vel[face_id][0] =  vel[cell_id][0]
                       + a * b_face_normal[face_id][0] / b_face_surf[face_id];
        bc_vel[face_id][1] =  vel[cell_id][1]
                       + a * b_face_normal[face_id][1] / b_face_surf[face_id];
        bc_vel[face_id][2] =  vel[cell_id][2]
                       + a * b_face_normal[face_id][2] / b_face_surf[face_id];
        /*  Total energy */
        bc_en[face_id] =  pinf / ((gamagp - 1.) * ro1)
                        + 0.5 * cs_math_3_square_norm(bc_vel[face_id]);

      }
      /*  Outlet */
      else {

        /*  Computation of the sound speed in state 1 */
        c1 = sqrt(gamagp * pinf / ro1);

        /*  Subsonic outlet - state 1 is imposed */
        if ((un1-c1) < 0.) {

          /*  Density */
          brom[face_id] = ro1;
          /*  Velocity */
          bc_vel[face_id][0] =  vel[cell_id][0]
                         + a * b_face_normal[face_id][0] / b_face_surf[face_id];
          bc_vel[face_id][1] =  vel[cell_id][1]
                         + a * b_face_normal[face_id][1] / b_face_surf[face_id];
          bc_vel[face_id][2] =  vel[cell_id][2]
                         + a * b_face_normal[face_id][2] / b_face_surf[face_id];
          /*  Total energy */
          bc_en[face_id] =  pinf / ((gamagp - 1.) * ro1)
                          + 0.5 * cs_math_3_square_norm(bc_vel[face_id]);

        }
        /*  Sonic outlet */
        else if ((uni-ci) < 0.) {

          /*  Mach number in the domain */
          mi = uni / ci;

          b = (gamagp - 1.) / (gamagp + 1.) * (mi + 2. / (gamagp - 1));

          /*  Sonic state pressure */
          bc_pr[face_id] = pri * pow(b, 2. * gamagp / (gamagp - 1.));
          /*  Sonic state density */
          brom[face_id] = roi * pow(b, 2. / (gamagp - 1.));
          /*  Sonic state velocity */
          uns = b * ci;
          bc_vel[face_id][0] = uns * b_face_normal[face_id][0] / b_face_surf[face_id];
          bc_vel[face_id][1] = uns * b_face_normal[face_id][1] / b_face_surf[face_id];
          bc_vel[face_id][2] = uns * b_face_normal[face_id][2] / b_face_surf[face_id];
          /*  Sonic state energy */
          bc_en[face_id] =  bc_pr[face_id]/((gamagp - 1.) * brom[face_id])
                          + 0.5 * uns*uns;

        }
        /*  Supersonic outlet */
        else {

          /*  pb = pri */
          bc_pr[face_id] = pri;
          /*  ub = uni */
          bc_vel[face_id][0] = vel[cell_id][0];
          bc_vel[face_id][1] = vel[cell_id][1];
          bc_vel[face_id][2] = vel[cell_id][2];
          /*  rob = roi */
          brom[face_id] = roi;
          /*  eb = ei */
          bc_en[face_id] = cvar_en[cell_id];

        }


      }

    }
    /*  Shock case */
    else {

      /*  Computation of the density in state 1 with Rankine-Hugoniot relations */
      ro1 = roi * ((gamagp - 1.) * pri  + (gamagp + 1.) * pinf)
                / ((gamagp - 1.) * pinf + (gamagp + 1.) * pri );

      /* Computation of the velocity in state 1 with Rankine-Hugoniot relations
         un1 = un2 */
      a = sqrt( (pinf - pri) * (1./roi - 1./ro1) );
      un1 = uni - a;

      /* Subsonic inlet - state 2 should be imposed but too few information
         is available to compute it
         for want of anything better, state 1 is imposed */
      if (un1 <= 0.) {

        /*  Density */
        brom[face_id] = ro1;
        /*  Velocity */
        bc_vel[face_id][0] =  vel[cell_id][0]
                       - a * b_face_normal[face_id][0] / b_face_surf[face_id];
        bc_vel[face_id][1] =  vel[cell_id][1]
                       - a * b_face_normal[face_id][1] / b_face_surf[face_id];
        bc_vel[face_id][2] =  vel[cell_id][2]
                       - a * b_face_normal[face_id][2] / b_face_surf[face_id];
        /*  Total energy */
        bc_en[face_id] =  pinf / ((gamagp-1.) * brom[face_id])
                        + 0.5 * cs_math_3_square_norm(bc_vel[face_id]);

      }
      /*  Outlet */
      else {

        /*  Computation of the shock velocity */
        sigma1 = (roi * uni - ro1 * un1) / (roi - ro1);

        /*  Subsonic outlet - state 1 is imposed */
        if (sigma1 <= 0.) {

          /*  Density */
          brom[face_id] = ro1;
          /*  Velocity */
          bc_vel[face_id][0] =  vel[cell_id][0]
                              - a * b_face_normal[face_id][0] / b_face_surf[face_id];
          bc_vel[face_id][1] =  vel[cell_id][1]
                              - a * b_face_normal[face_id][1] / b_face_surf[face_id];
          bc_vel[face_id][2] =  vel[cell_id][2]
                              - a * b_face_normal[face_id][2] / b_face_surf[face_id];
        /*  Total energy */
        bc_en[face_id] =  pinf / ((gamagp-1.) * brom[face_id])
                        + 0.5 * cs_math_3_square_norm(bc_vel[face_id]);

        }
        /*  Supersonic outlet */
        else {

          /*  pb = pri */
          bc_pr[face_id] = pri;
          /*  unb = uni */
          bc_vel[face_id][0] = vel[cell_id][0];
          bc_vel[face_id][1] = vel[cell_id][1];
          bc_vel[face_id][2] = vel[cell_id][2];
          /*  rob = roi */
          brom[face_id] = roi;
          /*  eb = ei */
          bc_en[face_id] = cvar_en[cell_id];

        }
        /*  test on shock speed sign */

      }
      /*  test on state 1 velocity sign */

    }
    /*  test on pinf-pri sign */

  }

}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute inlet boundary condition with total pressure and total
 * enthalpy imposed.
 *
 * \param[in,out] bc_en   total energy values at boundary faces
 * \param[in,out] bc_pr   pressure values at boundary faces
 * \param[in,out] bc_vel  velocity values at boundary faces
 * \param[in]     face_id    boundary face index
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_ph_inlet_bc(cs_real_t   *bc_en,
                         cs_real_t   *bc_pr,
                         cs_real_3_t *bc_vel,
                         cs_lnum_t    face_id)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;

  int niter, nitermax;
  cs_real_t gamagp, bMach, eps, pstat, old_pstat, ptot, res, rhotot;
  cs_real_t roi, ro1, pri, ei, uni, un1, y, uns, bc, cosalp, norm;
  cs_real_t ci, c1, mi, a, sigma1, utxi, utyi, utzi;
  cs_real_3_t dir;

  /*  Map field arrays */
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(u)->val;
  cs_real_t *cvar_pr = CS_F_(p)->val;
  cs_real_t *cvar_en = CS_F_(energy)->val;
  cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *brom = CS_F_(rho_b)->val;

  if (cs_glob_fluid_properties->ieos == 1) {

    /*  Calculation of the boundary conditions on the inlet face face_id */

    cs_cf_thermo_gamma(&gamagp);

    cs_lnum_t cell_id = b_face_cells[face_id];

    niter = 0;

    roi  = crom[cell_id];
    pri  = cvar_pr[cell_id];

    /*  Normalize the direction vector given by the user */
    norm = sqrt(cs_math_3_square_norm(bc_vel[face_id]));
    if (norm < cs_defs_epzero)
      bft_error(__FILE__, __LINE__, 0,
                _("Error in thermodynamics computations for compressible "
                  "flows\n"
                  "(perfect gas with constant Gamma):\n"
                  "The computation of the subsonic inlet boundary condition\n"
                  "with imposed total pressure and total enthalpy failed at\n"
                  "boundary face %i. The direction vector given by the user\n"
                  "can't be null."),
              face_id);

    dir[0] = bc_vel[face_id][0] / norm;
    dir[1] = bc_vel[face_id][1] / norm;
    dir[2] = bc_vel[face_id][2] / norm;

    /*  Angle between the imposed direction and the inlet normal */
    cosalp =  cs_math_3_dot_product(dir,b_face_normal[face_id])
      / b_face_surf[face_id];

    /*  If direction vector is outward, warn the user */
    if (cosalp > cs_defs_epzero)
      bft_printf("Warning in thermodynamics computations for compressible"
                   "flows\n"
                   "(perfect gas with constant Gamma):\n"
                   "The computation of the subsonic inlet boundary condition\n"
                   "with imposed total pressure and total enthalpy failed at\n"
                   "boundary face %i. The direction vector given by the user\n"
                   "points outward the fluid domain.\n",
                 face_id);

    /*  Computation of the sound speed inside the domain */
    ci = sqrt(gamagp * pri / roi);

    uni =  cs_math_3_dot_product(vel[cell_id],b_face_normal[face_id])
      / b_face_surf[face_id];

    bMach = uni / ci;

    utxi = vel[cell_id][0] - uni * b_face_normal[face_id][0]
      * b_face_surf[face_id];
    utyi = vel[cell_id][0] - uni * b_face_normal[face_id][1]
      * b_face_surf[face_id];
    utzi = vel[cell_id][0] - uni * b_face_normal[face_id][2]
      * b_face_surf[face_id];

    cs_real_t v2 = cs_math_3_square_norm(vel[cell_id]);
    ei   = cvar_en[cell_id] - 0.5 * v2;

    ptot = bc_pr[face_id];
    rhotot = gamagp / (gamagp - 1.) * ptot / bc_en[face_id];
    old_pstat = ptot;

    int key_cal_opt_id = cs_field_key_id("var_cal_opt");
    cs_var_cal_opt_t var_cal_opt;
    cs_field_get_key_struct(CS_F_(p), key_cal_opt_id, &var_cal_opt);

    eps = var_cal_opt.epsrsm;
    nitermax = 100;
    res = 1.;

    while (niter <= nitermax && res > eps) {

      pstat =  ptot*pow(1.+(gamagp - 1.)*0.5*bMach*bMach,gamagp/(1.-gamagp));
      y = pri / pstat;

      /*  1-shock */
      if (y < 1.) {

        /* Computation of the density in state 1 with Rankine-Hugoniot relations */
        ro1 = roi * ((gamagp - 1.) * pri   + (gamagp + 1.) * pstat)
          / ((gamagp - 1.) * pstat + (gamagp + 1.) * pri);

        /* Computation of the velocity in state 1 with Rankine-Hugoniot relations
           un1 = un2 */
        un1 = uni - sqrt( (pstat - pri) * (1./roi - 1./ro1) );

        /*  Subsonic inlet */
        if (un1 <= 0.) {

          /*  unb = u2 */
          bc_vel[face_id][0] = un1 / cosalp * dir[0];
          bc_vel[face_id][1] = un1 / cosalp * dir[1];
          bc_vel[face_id][2] = un1 / cosalp * dir[2];
          /*  rob = ro2 */
          brom[face_id] = pow(pstat/ptot,1./gamagp) * rhotot;
          /*  eb = e2 */
          bc_en[face_id] = pstat / ((gamagp - 1.) * brom[face_id])
                        + 0.5 * cs_math_3_square_norm(bc_vel[face_id]);
        }
        /*  Outlet */
        else {
          /*  Computation of the shock velocity */
          sigma1 = (roi * uni - ro1 * un1) / (roi - ro1);

          /*  subsonic outlet */
          if (sigma1 <= 0.) {

            /*  unb = u1 */
            bc_vel[face_id][0] = utxi + un1 * b_face_normal[face_id][0]
              / b_face_surf[face_id];
            bc_vel[face_id][1] = utyi + un1 * b_face_normal[face_id][1]
              / b_face_surf[face_id];
            bc_vel[face_id][2] = utzi + un1 * b_face_normal[face_id][2]
              / b_face_surf[face_id];
            /*  rob = ro1 */
            brom[face_id] = ro1;
            /*  eb = e1 */
            bc_en[face_id] =  ei
              - 0.5 * (pstat + pri) * (1. / ro1 - 1. / roi)
              + 0.5 * (un1*un1 + utxi*utxi + utyi*utyi + utzi*utzi);

          }
          /*  supersonic outlet */
          else {

            /*  pb = pri */
            pstat = pri;
            /*  unb = uni */
            bc_vel[face_id][0] = vel[cell_id][0];
            bc_vel[face_id][1] = vel[cell_id][1];
            bc_vel[face_id][2] = vel[cell_id][2];
            /*  rob = roi */
            brom[face_id] = roi;
            /*  eb = ei */
            bc_en[face_id] = cvar_en[cell_id];

          }


        }

      }
      /*  1-rarefaction */
      else {

        /* Computation of the velocity in state 1 using Riemann invariants
           of the 1-rarefaction */
        un1 =  uni +  2 * ci / (gamagp - 1.)
          * (1. - pow(pstat / pri,  (gamagp - 1.) / (2. * gamagp)));

        /* Computation of the density in state 1 using Riemann invariants
           of the 1-rarefaction */
        ro1 = pow(pstat / pri, 1. / gamagp) * roi;

        /*  Subsonic inlet */
        if (un1 <= 0.) {

          /*  unb = u2 */
          bc_vel[face_id][0] = un1 / cosalp * dir[0];
          bc_vel[face_id][1] = un1 / cosalp * dir[1];
          bc_vel[face_id][2] = un1 / cosalp * dir[2];
          /*  rob = ro2 */
          brom[face_id] = pow(pstat / ptot,1./gamagp) * rhotot;
          /*  eb = e2 */
          bc_en[face_id] = pstat / ((gamagp - 1.) * brom[face_id])
                        + 0.5 * cs_math_3_square_norm(bc_vel[face_id]);
        }
        /*  Outlet */
        else {

          /*  Computation of the sound speed in state 1 */
          c1 = sqrt(gamagp * pstat / ro1);

          /*  Subsonic outlet */
          if ((un1 - c1) < 0.) {

            /*  unb = u1 */
            bc_vel[face_id][0] = utxi + un1 * b_face_normal[face_id][0]
              / b_face_surf[face_id];
            bc_vel[face_id][1] = utyi + un1 * b_face_normal[face_id][1]
              / b_face_surf[face_id];
            bc_vel[face_id][2] = utzi + un1 * b_face_normal[face_id][2]
              / b_face_surf[face_id];
            /*  rob = ro1 */
            brom[face_id] = ro1;
            /*  eb = e1 */
            bc_en[face_id] =  pstat / (ro1 * (gamagp - 1.))
              + 0.5 * (un1*un1 + utxi*utxi + utyi*utyi + utzi*utzi);

          }
          /*  Supersonic outlet */
          else if ((uni - ci) >= 0.) {

            /*  pb = pri */
            pstat = pri;
            /*  ub = uni */
            bc_vel[face_id][0] = vel[cell_id][0];
            bc_vel[face_id][1] = vel[cell_id][1];
            bc_vel[face_id][2] = vel[cell_id][2];
            /*  rob = roi */
            brom[face_id] = roi;
            /*  eb = ei */
            bc_en[face_id] = cvar_en[cell_id];

          }
          /*  Outlet in sonic state */
          else {

            /*  Mach number in the domain */
            mi = uni / ci;

            a = (gamagp - 1.) / (gamagp + 1.) * (mi + 2. / (gamagp - 1));

            /*  Sonic state pressure */
            pstat = pri * pow(a,2.*gamagp/(gamagp-1.));
            /*  Sonic state density */
            brom[face_id] = roi * pow(a,2./(gamagp-1.));
            /*  Sonic state velocity */
            uns = a * ci;
            bc_vel[face_id][0] = uns * b_face_normal[face_id][0]
              / b_face_surf[face_id];
            bc_vel[face_id][1] = uns * b_face_normal[face_id][1]
              / b_face_surf[face_id];
            bc_vel[face_id][2] = uns * b_face_normal[face_id][2]
              / b_face_surf[face_id];
            /*  Sonic state energy */
            bc_en[face_id] =  pstat / ((gamagp - 1.) * brom[face_id])
              + 0.5 * uns*uns;

          }

        }

      }


      bc = sqrt(gamagp * pstat / brom[face_id]);
      bMach = cs_math_3_dot_product(bc_vel[face_id],b_face_normal[face_id])
              / b_face_surf[face_id] / bc;

      bc_pr[face_id] = pstat;

      /*  Pressure residual */
      res = CS_ABS((pstat - old_pstat) / ptot);

      /*  Prepare next iteration */
      old_pstat = pstat;
      niter++;
    }

    /*  Warn the user if fixed point algorithm did not converge */
    if (niter > nitermax)
      bft_printf("Warning in thermodynamics computations for compressible"
                   "flows\n"
                   "(perfect gas with constant Gamma):\n"
                   "Fixed point algorithm did not converge when computing\n"
                   "the subsonic inlet boundary condition with total\n"
                   "pressure and total enthalpy imposed.\n"
                   "At boundary face %i,\n"
                   "boundary Mach number residual = %12.4e,\n"
                   "maximum number of iterations (%i) was reached.\n",
                 face_id, res, nitermax);
  }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute epsilon sup:
 *
 * \f[\epsilon_{\textrm{sup}} = e - C_v T\f]
 *
 * for perfect gas: \f[\epsilon_{\textrm{sup}} = 0\f]

 * \param[out]    eps_sup epsilon sup array
 * \param[in]     l_size  l_size of the array
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo_eps_sup(cs_real_t *eps_sup,
                     cs_lnum_t  l_size)
{
  if (cs_glob_fluid_properties->ieos == 1) {
    /*  It is zero for a perfect gas */
    for (cs_lnum_t ii = 0; ii < l_size; ii++)
      eps_sup[ii] = 0.;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This subroutine is a driver allowing to call the appropriate
 * thermodynamicalfunctions depending on the quantities provided by the user.
 * Hence it is only used during the initialization step and at the boundaries
 * of type supersonic inlet. It is described in the following how to select the
 * quantity to be returned.
 *
 * When calling the user subroutine, the integer 'iccfth' specifies which
 * calculation has to be performed (and which quantity has to be returned).
 * The values for 'iccfth' for each case are provided below.
 *
 *   The variables are referred to using a different index i:
 *
 *     - pressure: 2
 *     - density: 3
 *     - temperature: 5
 *     - internal energy: 7
 *     - entropy: 13
 *
 *   iccfth is as follows, depending on which quantity needs to be computed:
 *     - variables at cell centers from variable i and variable j (i<j):
 *           iccfth = i*j*10000
 *     - variables at boundary faces from variable i and variable j (i<j):
 *           iccfth = i*j*10000+900
 *
 * Detailed values of iccfth and corresponding computations:
 *
 *   Values at the cell centers:
 *
 *     - temperature and energy from pressure and density: iccfth =  60000
 *     - density and energy from pressure and temperature: iccfth =  100000
 *     - density and temperature from pressure and energy: iccfth =  140000
 *     - pressure and energy from density and temperature: iccfth =  150000
 *     - pressure and temperature from density and energy: iccfth =  210000
 *
 *   Values at the faces for boundary conditions:
 *     - temperature and energy from pressure and density: iccfth = 60900
 *     - density and energy from pressure and temperature: iccfth = 100900
 *     - density and temperature from pressure and energy: iccfth = 140900
 *     - pressure and energy from density and temperature: iccfth = 150900
 *     - pressure and temperature from density and energy: iccfth = 210900
 *
 * \param[in]     iccfth        id of computation
 * \param[in]     face_id       face index if the computation is for a B.C.
 * \param[in,out] bc_en         total energy values at boundary faces
 * \param[in,out] bc_pr         pressure values at boundary faces
 * \param[in,out] bc_tk         temperature values at boundary faces
 * \param[in,out] bc_vel        velocity values at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_thermo(const int    iccfth,
             cs_lnum_t    face_id,
             cs_real_t   *bc_en,
             cs_real_t   *bc_pr,
             cs_real_t   *bc_tk,
             cs_real_3_t *bc_vel)
{
  /* Local variables */
  cs_lnum_t l_size;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  /*  Map field arrays */
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(u)->val;
  cs_real_t *cvar_pr = (cs_real_t *)CS_F_(p)->val;
  cs_real_t *crom = (cs_real_t *)CS_F_(rho)->val;
  cs_real_t *brom = (cs_real_t *)CS_F_(rho_b)->val;
  cs_real_t *cvar_tk = (cs_real_t *)CS_F_(t)->val;
  cs_real_t *cvar_en = (cs_real_t *)CS_F_(energy)->val;

  /*  0. Initialization. */

  /*  Calculation of temperature and energy from pressure and density */
  if (iccfth == 60000) {
    l_size = n_cells;
    cs_cf_check_density(crom, l_size);
    cs_cf_thermo_te_from_dp(cvar_pr, crom, cvar_tk, cvar_en, vel, l_size);
  }
  /*  Calculation of density and energy from pressure and temperature: */
  else if (iccfth == 100000) {
    l_size = n_cells;
    cs_cf_check_temperature(cvar_tk, l_size);
    cs_cf_thermo_de_from_pt(cvar_pr, cvar_tk, crom, cvar_en, vel, l_size);
  }
  /*  Calculation of density and temperature from pressure and energy */
  else if (iccfth == 140000) {
    l_size = n_cells;
    cs_cf_thermo_dt_from_pe(cvar_pr, cvar_en, crom, cvar_tk, vel, l_size);
  }
  /*  Calculation of pressure and energy from density and temperature */
  else if (iccfth == 150000) {
    l_size = n_cells;
    cs_cf_thermo_pe_from_dt(crom, cvar_tk, cvar_pr, cvar_en, vel, l_size);
  }
  /*  Calculation of pressure and temperature from density and energy */
  else if (iccfth == 210000) {
    l_size = n_cells;
    cs_cf_thermo_pt_from_de(crom, cvar_en, cvar_pr, cvar_tk, vel, l_size);
  }
  /*  Calculation of temperature and energy from pressure and density
      (it is postulated that the pressure and density values are strictly
      positive) */
  else if (iccfth == 60900) {
    l_size = 1;
    cs_cf_thermo_te_from_dp(&bc_pr[face_id], &brom[face_id], &bc_tk[face_id],
                            &bc_en[face_id], &bc_vel[face_id], l_size);
  }
  /*  Calculation of density and energy from pressure and temperature */
  else if (iccfth == 100900) {
    l_size = 1;
    cs_cf_thermo_de_from_pt(&bc_pr[face_id], &bc_tk[face_id], &brom[face_id],
                            &bc_en[face_id], &bc_vel[face_id], l_size);
  }
  /*  Calculation of density and temperature from pressure and total energy */
  else if (iccfth == 140900) {
    l_size = 1;
    cs_cf_thermo_dt_from_pe(&bc_pr[face_id], &bc_en[face_id], &brom[face_id],
                            &bc_tk[face_id], &bc_vel[face_id], l_size);
  }
  /*  Calculation of pressure and energy from density and temperature */
  else if (iccfth == 150900) {
    l_size = 1;
    cs_cf_thermo_pe_from_dt(&brom[face_id], &bc_tk[face_id], &bc_pr[face_id],
                            &bc_en[face_id], &bc_vel[face_id], l_size);
  }
  /*  Calculation of pressure and temperature from density and energy */
  else if (iccfth == 210900) {
    l_size = 1;
    cs_cf_thermo_pt_from_de(&brom[face_id], &bc_en[face_id], &bc_pr[face_id],
                            &bc_tk[face_id], &bc_vel[face_id], l_size);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
