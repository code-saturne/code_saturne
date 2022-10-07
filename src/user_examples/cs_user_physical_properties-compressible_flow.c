/*============================================================================
 * User definition of physical properties.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_physical_properties.c
 *
 * \brief User definition of physical properties.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties(cs_domain_t   *domain)
{
  /* List of examples

   * Ex. 1: molecular viscosity varying with temperature
   * Ex. 2: molecular volumetric viscosity varying with temperature
   * Ex. 3: isobaric specific heat varying with temperature
   * Ex. 4: molecular thermal conductivity varying with temperature
   * Ex. 5: molecular diffusivity of user-defined scalars varying
   *        with temperature */

  /*! [compressible_properties_init] */
  const cs_lnum_t n_cells = domain->mesh->n_cells;

  cs_real_t *cpro_cp = NULL;
  cs_real_t *cpro_cv = NULL;
  cs_real_t *cpro_vtmpk = NULL;
  cs_real_t *cpro_viscl = CS_F_(mu)->val;
  cs_real_t *cvar_t = CS_F_(t_kelvin)->val;
  cs_real_t *mix_mol_mas = cs_field_by_name("mix_mol_mas")->val;

  /* Molecular volumetric viscosity */
  cs_real_t *cpro_viscv = cs_field_by_name_try("volume_viscosity")->val;

  const int kivisl = cs_field_key_id("diffusivity_id");
  int ifcvsl = cs_field_get_key_int(CS_F_(t_kelvin), kivisl);
  if (ifcvsl > -1)
    cpro_vtmpk = cs_field_by_id(ifcvsl)->val;

  if (CS_F_(cp) != NULL)
    cpro_cp = CS_F_(cp)->val;

  if (CS_F_(cv) != NULL)
      cpro_cv = CS_F_(cv)->val;
  /*! [compressible_properties_init] */

  /* Ex. 1: molecular viscosity varying with temperature
   * =====
   */

  {
    /*! [example_1] */

    /* Molecular dynamic viscosity 'cpro_viscl' */

    /*  User-defined coefficients for the selected law.
     *  The values hereafter are provided as a mere example. They
     *  are physically meaningless. */

    cs_real_t varam = -3.4016e-9;
    cs_real_t varbm =  6.2332e-7;
    cs_real_t varcm = -4.5577e-5;
    cs_real_t vardm =  1.6935e-3;

    /* Molecular dynamic viscosity mu at the cell centers, kg/(m s)
     * In this example, mu is provided as a function of the temperature T:
     *   mu(T)           = T  *( T  *( am  * T +  bm  )+ cm  )+ dm
     * that is:
     *   cpro_viscl(iel) = xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t xvart = cvar_t[c_id];
      cpro_viscl[c_id] = xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm;
    }
    /*! [example_1] */
  }

  /* Ex. 2: molecular volumetric viscosity varying with temperature
   * =====
   */

  {
    /*! [example_2] */

    /* Stop if the volumetric viscosity has not been defined as variable */
    if (cpro_viscv == NULL)
      bft_error(__FILE__, __LINE__, 0,
                "%s: cpro_viscv not available.", __func__);

    /* User-defined coefficients for the selected law.
     * The values provided hereafter are provided as a mere example.
     * They are physically meaningless. */

    const cs_real_t varam = -3.4016e-9;
    const cs_real_t varbm =  6.2332e-7;
    const cs_real_t varcm = -4.5577e-5;
    const cs_real_t vardm =  1.6935e-3;

    /* Molecular dynamic volumetric viscosity kappa at the cell centers, kg/(m s)
     * In this example, kappa is provided as a function of the temperature T:
     *   kappa(T)        = T  *( T  *( am  * T +  bm  )+ cm  )+ dm
     * that is:
     *   cpro_viscv(iel) = xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm */

    for (cs_lnum_t c_id =0; c_id < n_cells; c_id++) {
      const cs_real_t xvart = cvar_t[c_id];
      cpro_viscv[c_id] = xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm;
    }
    /*! [example_2] */
  }

  /* Ex. 3: isobaric specific heat varying with temperature
   * =====
   */

  {
    /*! [example_3] */

    /* Stop if the isobaric or isochoric specific heat (cpro_cp or cpro_cv)
     * has not been defined as variable */

    if ((cpro_cp == NULL) || (cpro_cv == NULL))
      bft_error(__FILE__, __LINE__, 0,
                "%s: cpro_cp or cpro_cv not available.", __func__);

    /* User-defined coefficients for the selected law.
     * The values provided hereafter are provided as a mere example.
     * They are physically meaningless.*/

    const cs_real_t varac = 0.00001;
    const cs_real_t varbc = 1000.0;

    /* Isobaric specific heat cpro_cp at the cell centers, J/(kg degree)
     * In this example, cpro_cp is provided as a function of the temperature T:
     *   cpro_cp(T)   = ac * T  + ab
     * that is:
     *   cpro_cp(iel) = varac*xvart+varbc */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t  xvart = cvar_t[c_id];
      cpro_cp[c_id] = varac*xvart + varbc;
    }

    /* The isochoric specific heat is deduced from the isobaric specific heat */
    cs_cf_thermo_cv(cpro_cp, mix_mol_mas, cpro_cv, n_cells);

    /*! [example_3] */
  }

  /* Ex. 4: molecular thermal conductivity varying with temperature
   * =====
   */

  {
    /*! [example_4] */
    /* Stop if the molecular thermal conductivity has not
     * been defined as variable */

    if (cpro_vtmpk == NULL)
      bft_error(__FILE__, __LINE__, 0,
                "%s: cpro_vtmpk not available.", __func__);

    /* User-defined coefficients for the selected law.
     * The values provided hereafter are provided as a mere example.
     * They are physically meaningless. */

    cs_real_t varal = -3.3283e-7;
    cs_real_t varbl =  3.6021e-5;
    cs_real_t varcl =  1.2527e-4;
    cs_real_t vardl =  0.589230;

    /* Molecular thermal conductivity lambda at the cell centers, W/(m degree)
     * In this example, lambda is provided as a function of the temperature T:
     *   lambda(T)          =    T  *( T  *( al  * T +  bl  )+ cl  )+ dl
     * that is:
     *   cpro_vtmpk(iel) =   xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t xvart = cvar_t[c_id];
      cpro_vtmpk[c_id] = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl);
    }
    /*!  [example_4] */
  }

  /* Ex. 5: molecular diffusivity of user-defined scalars varying
   * =====
   *        with temperature
   * The molecular diffusivity can be set for all the user-defined scalars
   * ** except **:
   *  - temperature and enthalpy (already dealt with above: for these
   *    variables, the 'diffusivity' is the thermal conductivity)
   *  - variances of the fluctuations of another scalar variable (the
   *    diffusivity is assumed to be equal to that of the associated
   *    scalar)
   * The values of the molecular diffusivity are provided as a function
   * of the temperature. All variables are evaluated at the cell centers. */

  {
    /*! [example_5] */
    const int n_fields = cs_field_n_fields();
    const int keysca = cs_field_key_id("scalar_id");
    const int kscavr = cs_field_key_id("first_moment_id");

    /* Loop on the scalars fields */
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *fld = cs_field_by_id(f_id);

      /* If the scalar is the temperature, it will be skipped. */
      if (fld == CS_F_(t_kelvin))
        continue;

      /* Here we only handle user or model scalar-type variables
         which are not fluctuations */

      int sc_id = -1;
      if (fld->type & CS_FIELD_VARIABLE)
        sc_id = cs_field_get_key_int(fld, keysca) - 1;
      if (sc_id < 0)
        continue;

      int variance_id = cs_field_get_key_int(fld, kscavr);
      int diffusivity_id = cs_field_get_key_int(fld, kivisl);

      if (variance_id > -1 || diffusivity_id < 0)
        continue;

      cs_real_t *cpro_vscal = cs_field_by_id(diffusivity_id)->val;

      /* User-defined coefficients for the selected law.
       * The values provided hereafter are provided as a mere example.
       * They are physically meaningless. */

      const cs_real_t varal = -3.3283e-7;
      const cs_real_t varbl =  3.6021e-5;
      const cs_real_t varcl =  1.2527e-4;
      const cs_real_t vardl =  0.5892300;

      /* Molecular diffusivity lambda at the cell centers, kg/(m s)
       * In this example, lambda is provided as a function of the temperature T:
       *   lambda(T)       = T  *( T  *( al  * T +  bl  )+ cl  )+ dl
       * that is:
       *   cpro_vscal(iel) = xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl */

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t xvart = cvar_t[c_id];
        cpro_vscal[c_id] = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl);
      }

    } /* End of the loop on the scalars */
    /*! [example_5] */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
