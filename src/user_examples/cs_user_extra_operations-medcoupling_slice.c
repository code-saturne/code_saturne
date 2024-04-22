/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
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
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 *  User Header
 *----------------------------------------------------------------------------*/

#include "cs_user_profile.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function)
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize variables.
 *
 * This function is called at beginning of the computation
 * (restart or not) before the time step loop.
 *
 * This is intended to initialize or modify (when restarted)
 * variable and time step values.

 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations_initialize(cs_domain_t *domain)
{
  CS_UNUSED(domain);

  {
    /*![medcpl_slice_init_1]*/
    /* Add a slice plane intersecting the origin (0,0,0) with a normal
     * along the Z axis (0,0,1), and length of 1m along both axis.
     *
     * If code_saturne is compiled without MEDCoupling, this call will
     * raise an error during runtime.
     */
    cs_real_t o[3] = {0.};
    cs_real_t n[3] = {0., 0., 1.};
    cs_medcoupling_postprocess_add_plane_slice("slice_OZ",  /* Name of the slice */
                                               "all[]",     /* Selection criteria for cells to intersect */
                                               o,           /* Slice origin point */
                                               n,           /* Slice normal vector */
                                               1.,          /* Plane length along first axis */
                                               1.);         /* Plane length along second axis */
    /*![medcpl_slice_init_1]*/
  }

  {
    /*!(medcpl_slice_activate_postprocessing]*/
    /* Activate intersection flag and surface postprocessing functions.
     */
    cs_medcoupling_slice_activate_postprocess("slice_OZ");
    /*!(medcpl_slice_activate_postprocessing]*/
  }

  {
    /*![medcpl_slice_init_2]*/
    /* Add a disc slice intersecting the origin (0,0,0) with a normal
     * along the Z axis (0,0,1), and a radius of 0.5m.
     *
     * If code_saturne is compiled without MEDCoupling, this call will
     * raise an error during runtime.
     */
    cs_real_t o[3] = {0.};
    cs_real_t n[3] = {0., 0., 1.};
    cs_medcoupling_postprocess_add_disc_slice("disc1",  /* Name of the slice */
                                              "all[]",     /* Selection criteria for cells to intersect */
                                               o,           /* Slice origin point */
                                               n,           /* Slice normal vector */
                                               0.5,         /* Disc radius */
                                               -1);         /* Number of sectors. if < 0 default value of 36 is used. */
    /*![medcpl_slice_init_2]*/
  }

  {
    /*![medcpl_slice_init_3]*/
    /* Add a slice plane intersecting the origin (0,0,0) with a normal
     * along the Z axis (0,0,1), and length of 1m along both axis.
     *
     * If code_saturne is compiled without MEDCoupling, this call will
     * raise an error during runtime.
     */
    cs_real_t o[3] = {0.};
    cs_real_t n[3] = {0., 0., 1.};
    cs_medcoupling_postprocess_add_annulus_slice("annulus_slice", /* Name of the slice */
                                                 "all[]",     /* Selection criteria for cells to intersect */
                                                 o,           /* Slice origin point */
                                                 n,           /* Slice normal vector */
                                                 0.2,         /* Inner radius (hole)  */
                                                 0.5,         /* Outer radius */
                                                 72);         /* Number of sectors. if < 0 default value of 36 is used */
    /*![medcpl_slice_init_3]*/
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t *domain)
{

  {
    /*![medcpl_slice_mean]*/
    /* Compute mean of a scalar over the slice "slice_OZ" */
    cs_real_t *cpro_scal = cs_field_by_name("my_scalar")->val;

    cs_real_t mean_val = cs_medcoupling_slice_scalar_mean("slice_OZ", cpro_scal);
    /*![medcpl_slice_mean]*/
  }

  {
    /*![medcpl_slice_integral]*/
    /* Compute an integral of a given scalar over the slice "annulus_slice" */
    cs_real_t *cpro_scal = cs_field_by_name("my_scalar")->val;

    cs_real_t int_val =
      cs_medcoupling_slice_scalar_integral("annulus_slice", /* Name of the slice */
                                           cpro_scal);      /* Pointer to scalar values */
    /*![medcpl_slice_integral]*/
  }

  {
    /*![medcpl_slice_integral_weighted]*/
    /* Compute Tbulk of a disc slice.
     * Tbulk = Int(T * cp * rho * <u|n>) / Int(cp * rho * <u|n>)
     */

    const cs_lnum_t n_cells = domain->mesh->n_cells;
    cs_real_t *rho_cp = NULL;
    BFT_MALLOC(rho_cp, n_cells, cs_real_t);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rho_cp[c_id] = CS_F_(cp)->val[c_id] * CS_F_(rho)->val[c_id];

    cs_real_3_t *cvar_vel = (cs_real_3_t *)CS_F_(vel)->val;

    cs_real_t Tbulk = cs_medcoupling_slice_scalar_mean_weighted("disc1",
                                                                CS_F_(t)->val,
                                                                rho_cp,
                                                                cvar_vel);

    BFT_FREE(rho_cp);
    /*![medcpl_slice_integral_weighted]*/
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of the calculation.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.

 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations_finalize(cs_domain_t *domain)
{
  CS_UNUSED(domain);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
