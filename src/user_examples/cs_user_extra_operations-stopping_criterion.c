/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
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
#include <stdio.h>

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
 * \file cs_user_extra_operations-stopping_criterion.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function).
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Example of extra operations allowing to properly stop a computation
 * when the L2 time residuals (displayed in the run_solver.log file) of all
 * solved variables have decreased below a value of 1e-3.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  /*! [extra_stopping_criterion] */

  /* get total number of fields */
  int n_fields = cs_field_n_fields();

  /* get time step structure */
  cs_time_step_t *ts = cs_get_glob_time_step();

  /* declare a C structure holding solving information values */
  cs_solving_info_t sinfo;

  bool cved = true;

  /* criterion is set here to 1e-3 */
  cs_real_t epsres = 1e-3;

  /* loop over all fields */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    /* filter fields of type variable (i.e. the solved ones) */
    if (f->type & CS_FIELD_VARIABLE) {
      /* get solving info structure for current field */
      cs_field_get_key_struct(f, cs_field_key_id("solving_info"), &sinfo);

      /* has the value of the L2 time residual decreased below criterion ? */
      cved = cved && (sinfo.l2residual  <= epsres);
    }
  }

  /* if current iteration is the second to last one, nothing to do */
  cved = cved && (ts->nt_cur < ts->nt_max);

  /* otherwise if converged */
  if (cved) {
    /* set maximum number of iteration to the next one */
    ts->nt_max = ts->nt_cur+1;

    /* Warning: L2 time residual is computed after extra operations */
    bft_printf("Converged at %d\n",ts->nt_max);
  }

  /*! [extra_stopping_criterion] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
