/*============================================================================
 * Checkpoint/restart handling for default application.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_restart_default.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_restart_default.c
        Checkpoint/restart handling for default application.
*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

const char *_coeff_name[] = {"bc_coeffs::a", "bc_coeffs::b",
                             "bc_coeffs::af", "bc_coeffs::bf",
                             "bc_coeffs::ad", "bc_coeffs::bd",
                             "bc_coeffs::ac", "bc_coeffs::bc"};

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */


/*! \endcond (end ignore by Doxygen) */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read boundary condition coefficients for all fields from checkpoint.
 *
 * \param[in, out]  restart  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_read_bc_coeffs(cs_restart_t  *restart)
{
  int c_id, f_id;

  int errcount = 0;
  const int coupled_key_id = cs_field_key_id_try("coupled");
  const int n_fields = cs_field_n_fields();

  /* Loop on all fields, to search for those defined on all cells
     and with BC coefficients */

  for (f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    if (   f->location_id == CS_MESH_LOCATION_CELLS
        && f->bc_coeffs != NULL) {

      /* Check for presence of coefficients */

      int coupled = 0;
      int n_loc_vals = 1;

      int32_t coeff_p[] = {0, 0, 0, 0, 0, 0, 0, 0};

      cs_real_t *p[] = {f->bc_coeffs->a,
                        f->bc_coeffs->b,
                        f->bc_coeffs->af,
                        f->bc_coeffs->bf,
                        f->bc_coeffs->ad,
                        f->bc_coeffs->bd,
                        f->bc_coeffs->ac,
                        f->bc_coeffs->bc};

      for (c_id = 0; c_id < 8; c_id++) {
        if (p[c_id] != NULL) {
          coeff_p[c_id] = 1;
          /* avoid double reads/writes in case of aliasing */
          for (int i = 0; i < c_id; i++) {
            if (p[i] == p[c_id])
              coeff_p[c_id] = 0;
          }
        }
      }

      cs_parall_max(8, CS_INT32, coeff_p);

      if (f->dim > 1 && coupled_key_id > -1)
        coupled = cs_field_get_key_int(f, coupled_key_id);

      for (c_id = 0; c_id < 8; c_id++) {

        int retval;
        char *sec_name = NULL;
        cs_real_t *c = p[c_id];

        if (coeff_p[c_id] == 0)
          continue;

        if (coupled) {
          if (c_id %2 == 0)
            n_loc_vals = f->dim;
          else
            n_loc_vals = f->dim * f->dim;
        }
        else { /* uncoupled */
          n_loc_vals = f->dim;
          if (f->dim > 1 && !f->interleaved) { /* Interleave values if not done yet */
            const cs_lnum_t *n_elts
              = cs_mesh_location_get_n_elts(CS_MESH_LOCATION_BOUNDARY_FACES);
            BFT_MALLOC(c, f->dim*n_elts[0], cs_real_t);
          }
        }

        BFT_MALLOC(sec_name,
                   strlen(f->name) + strlen(_coeff_name[c_id]) + 3,
                   char);
        sprintf(sec_name, "%s::%s", f->name, _coeff_name[c_id]);

        retval = cs_restart_read_section(restart,
                                         sec_name,
                                         3, /* location_id */
                                         n_loc_vals,
                                         CS_TYPE_cs_real_t,
                                         c);

        if (retval != CS_RESTART_SUCCESS)
          errcount += 1;

        BFT_FREE(sec_name);

        if (f->dim > 1 && !f->interleaved && coupled == 0) {

          /* De-interleave values (obsolete case) */

          const cs_lnum_t *n_elts
            = cs_mesh_location_get_n_elts(CS_MESH_LOCATION_BOUNDARY_FACES);
          cs_lnum_t _n_elts = n_elts[0];
          BFT_MALLOC(c, f->dim*_n_elts, cs_real_t);
          for (cs_lnum_t j = 0; j < _n_elts; j++) {
            for (int k = 0; k < f->dim; k++)
              p[c_id][j + k*n_elts[2]] = c[j*f->dim + k];
          }
          BFT_FREE(c);

        }

      } /* End of loop in i (coeff type) */

    } /* End for field with BC coeffs */

  } /* End of loop on fields */

  if (errcount > 0) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("\nSome boundary condition coefficients "
                 "could not be read from a restart file;\n"
                 "they will be initialized with default values.\n\n"));
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write boundary condition coefficients for all fields to checkpoint.
 *
 * \param[in, out]  restart  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_write_bc_coeffs(cs_restart_t  *restart)
{
  int c_id, f_id;

  const int coupled_key_id = cs_field_key_id_try("coupled");
  const int n_fields = cs_field_n_fields();

  /* Loop on all fields, to search for those defined on all cells
     and with BC coefficients */

  for (f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    if (   f->location_id == CS_MESH_LOCATION_CELLS
        && f->bc_coeffs != NULL) {

      /* Check for presence of coefficients */

      int coupled = 0;
      int n_loc_vals = 1;

      int32_t coeff_p[] = {0, 0, 0, 0, 0, 0, 0, 0};

      cs_real_t *p[] = {f->bc_coeffs->a,
                        f->bc_coeffs->b,
                        f->bc_coeffs->af,
                        f->bc_coeffs->bf,
                        f->bc_coeffs->ad,
                        f->bc_coeffs->bd,
                        f->bc_coeffs->ac,
                        f->bc_coeffs->bc};

      for (c_id = 0; c_id < 8; c_id++) {
        if (p[c_id] != NULL) {
          coeff_p[c_id] = 1;
          /* avoid double reads/writes in case of aliasing */
          for (int i = 0; i < c_id; i++) {
            if (p[i] == p[c_id])
              coeff_p[c_id] = 0;
          }
        }
      }

      cs_parall_max(8, CS_INT32, coeff_p);

      if (f->dim > 1 && coupled_key_id > -1)
        coupled = cs_field_get_key_int(f, coupled_key_id);

      for (c_id = 0; c_id < 8; c_id++) {

        char *sec_name = NULL;

        cs_real_t *c = p[c_id];

        if (coeff_p[c_id] == 0)
          continue;

        if (coupled) {
          if (c_id %2 == 0)
            n_loc_vals = f->dim;
          else
            n_loc_vals = f->dim * f->dim;
        }
        else { /* uncoupled */
          n_loc_vals = f->dim;
          if (f->dim > 1 && !f->interleaved) { /* Interleave values if not done yet */
            const cs_lnum_t *n_elts
              = cs_mesh_location_get_n_elts(CS_MESH_LOCATION_BOUNDARY_FACES);
            cs_lnum_t _n_elts = n_elts[0];
            BFT_MALLOC(c, f->dim*_n_elts, cs_real_t);
            for (cs_lnum_t j = 0; j < _n_elts; j++) {
              for (int k = 0; k < f->dim; k++)
                c[j*f->dim + k] = p[c_id][j + k*n_elts[2]];
            }
          }

        }

        BFT_MALLOC(sec_name,
                   strlen(f->name) + strlen(_coeff_name[c_id]) + 3,
                   char);
        sprintf(sec_name, "%s::%s", f->name, _coeff_name[c_id]);

        cs_restart_write_section(restart,
                                 sec_name,
                                 3, /* location_id */
                                 n_loc_vals,
                                 CS_TYPE_cs_real_t,
                                 c);

        BFT_FREE(sec_name);

        if (c != p[c_id])
          BFT_FREE(c);

      } /* End of loop in i (coeff type) */

    } /* End for field with BC coeffs */

  } /* End of loop on fields */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
