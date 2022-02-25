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
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

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
  CS_NO_WARN_IF_UNUSED(domain);

  /* Check fields exists */
  if (CS_F_(lambda) == NULL)
    bft_error(__FILE__, __LINE__, 0,_("error lambda not variable\n"));
  if (CS_F_(rho) == NULL)
    bft_error(__FILE__, __LINE__, 0,_("error rho not variable\n"));
  if (CS_F_(cp) == NULL)
    bft_error(__FILE__, __LINE__, 0,_("error cp not variable\n"));

  cs_real_t *cpro_lambda = CS_F_(lambda)->val;
  cs_real_t *cpro_rho = CS_F_(rho)->val;
  cs_real_t *cpro_cp = CS_F_(cp)->val;

  /* Impose thermal conductivity, density and specific heat for solid zones */
  {
    /* Volume zone "CM" must be defined in the GUI or in cs_user_zones.c */
    const cs_zone_t *z = cs_volume_zone_by_name("CM");

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {
      cs_lnum_t cell_id = z->elt_ids[i];
      cpro_lambda[cell_id] = 3.3;
      cpro_rho[cell_id] = 7800.;
      cpro_cp[cell_id] = 444.;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of enthalpy to temperature conversion.
 *
 * This allows overwriting the solver defaults if necessary.
 *
 * This function may be called on a per-zone basis, so as to allow different
 * conversion relations in zones representing solids or different fluids.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       z        zone (volume or boundary) applying to current call
 * \param[in]       z_local  if true, h and t arrays are defined in a compact
 *                           (contiguous) manner for this zone only;
 *                           if false, h and t are defined on the zone's parent
 *                           location (usually all cells or boundary faces)
 * \param[in]       h        enthalpy values
 * \param[in, out]  t        temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_h_to_t(cs_domain_t      *domain,
                                   const cs_zone_t  *z,
                                   bool              z_local,
                                   const cs_real_t   h[],
                                   cs_real_t         t[])
{
  CS_NO_WARN_IF_UNUSED(domain);
  CS_NO_WARN_IF_UNUSED(z);

  /* Tabulated values */

  /*! [tabulation] */
  static const int n_tv = 5;
  static const cs_real_t ht[5] = {100000.0, 200000.0, 300000.0,
                                  400000.0, 00000.0};
  static const cs_real_t th[5] = {100.0, 200.0, 300.0, 400.0, 500.0};
  /*! [tabulation] */

  /* Conversion:
     Note that z->name or z->location_id can be used as a filter
     if "per-zone" properties are needed (such as with solid zones) */

  /*! [z_h_to_t] */
  for (cs_lnum_t i_l = 0; i_l < z->n_elts; i_l++) {

    cs_lnum_t i = (z_local) ? i_l : z->elt_ids[i_l];

    cs_real_t temperature = 0;  /* Default initialization */

    /* If H is outside the tabulated value range, use range
       start or end value. */

    if (h[i] < ht[0])
      temperature = th[0];
    else if (h[i] > ht[n_tv - 1])
      temperature = th[n_tv - 1];

    /* Otherwise, use piecewise linear interpolation */

    else
      for (int j = 1; j < n_tv; j++) {
        if (h[j] < ht[j]) {
          temperature = th[j-1] +   (h[i]-ht[j-1])*(th[j]-th[j-1])
                                  / (ht[j]-ht[j-1]);
          break;
        }
      }

    t[i] = temperature;

  }
  /*! [z_h_to_t] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of temperature to enthalpy conversion.
 *
 * This allows overwriting the solver defaults if necessary.
 *
 * This function may be called on a per-zone basis, so as to allow different
 * conversion relations in zones representing solids or different fluids.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       z        zone (volume or boundary) applying to current call
 * \param[in]       z_local  if true, h and t arrays are defined in a compact
 *                           (contiguous) manner for this zone only;
 *                           if false, h and t are defined on the zone's parent
 *                           location (usually all cells or boundary faces)
 * \param[in]       h        temperature values
 * \param[in, out]  t        enthalpy values
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_t_to_h(cs_domain_t      *domain,
                                   const cs_zone_t  *z,
                                   bool              z_local,
                                   const cs_real_t   t[],
                                   cs_real_t         h[])
{
  CS_NO_WARN_IF_UNUSED(domain);
  CS_NO_WARN_IF_UNUSED(z);
  CS_NO_WARN_IF_UNUSED(domain);
  CS_NO_WARN_IF_UNUSED(z);

  /* Tabulated values */

  static const int n_tv = 5;
  static const cs_real_t ht[5] = {100000.0, 200000.0, 300000.0,
                                  400000.0, 00000.0};
  static const cs_real_t th[5] = {100.0, 200.0, 300.0, 400.0, 500.0};

  /* Conversion:
     Note that z->name or z->location_id can be used as a filter
     if "per-zone" properties are needed (such as with solid zones) */

  /*! [z_t_to_h] */
  for (cs_lnum_t i_l = 0; i_l < z->n_elts; i_l++) {

    cs_lnum_t i = (z_local) ? i_l : z->elt_ids[i_l];

    cs_real_t enthalpy = 0;  /* Default initialization */

    /* If H is outside the tabulated value range, use range
       start or end value. */

    if (t[i] < th[0])
      enthalpy = ht[0];
    else if (t[i] > th[n_tv - 1])
      enthalpy = ht[n_tv - 1];

    /* Otherwise, use piecewise linear interpolation */

    else
      for (int j = 1; j < n_tv; j++) {
        if (t[j] < th[j]) {
          enthalpy = ht[j-1] +   (t[i]-th[j-1])*(ht[j]-ht[j-1])
                               / (th[j]-th[j-1]);
          break;
        }
      }

    h[i] = enthalpy;

  }
  /*! [z_t_to_h] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User modification of the Smagorinsky constant for the
 *        dynamic Smagorinsky model.
 *
 * CS = Mij.Lij / Mij.Mij
 *
 * The local averages of the numerator and denominator are done before calling
 * this function, so
 *
 * CS = < Mij.Lij > / < Mij.Mij >
 *
 * In this subroutine, Mij.Lij and Mij.Mij are passed as arguments before
 * the local average.
 *
 * \param[in, out]   domain      pointer to a cs_domain_t structure
 * \param[in]        mijlij      mij.lij before the local averaging
 * \param[in]        mijmij      mij.mij before the local averaging
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_smagorinsky_c(cs_domain_t      *domain,
                                          const cs_real_t   mijlij[],
                                          const cs_real_t   mijmij[])
{
  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_real_t tot_vol = domain->mesh_quantities->tot_vol;
  const cs_real_t *cell_vol = domain->mesh_quantities->cell_vol;

  cs_real_t *cpro_smago
    = cs_field_by_name("smagorinsky_constant^2")->val;

  cs_real_t mijmijmoy = 0, mijlijmoy = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    mijlijmoy += mijlij[c_id]*cell_vol[c_id];
    mijmijmoy += mijmij[c_id]*cell_vol[c_id];
  }

  cs_parall_sum(1, CS_REAL_TYPE, &mijlijmoy);
  cs_parall_sum(1, CS_REAL_TYPE, &mijmijmoy);

  mijmijmoy /= tot_vol;
  mijlijmoy /= tot_vol;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    cpro_smago[c_id] = mijlijmoy/mijmijmoy;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
