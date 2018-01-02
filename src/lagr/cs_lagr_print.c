/*============================================================================
 * Methods for particle info printing
 *============================================================================*/

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

/*============================================================================
 * Functions dealing with file printing
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_parall.h"
#include "cs_time_step.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_tracking.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_print.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

static FILE *flal = NULL; /* associated file */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write Lagrangian particle info file
 *
 * This file logs, for each time step:
 *  - number of particles in the domain
 *  - number of entered particles
 *  - number of exited particles
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_print(cs_real_t ttcabs)
{
  static cs_lnum_t _ipass = 0;  /* number of calls */

  _ipass++;

  const cs_lagr_model_t *lagr_model = cs_glob_lagr_model;

  /* Parallelism management    */

  const cs_lagr_particle_counter_t *pc = cs_lagr_update_particle_counter();

  /* Open file at first pass */

  if (   cs_glob_rank_id <= 0
      && flal == NULL && _ipass == 1)
    flal = fopen("listla","w");

  /* ====================================================================   */
  /* 2. OUVERTURE DU FICHIER DE STOCKAGE */
  /* ====================================================================   */
  /*     Seul le premier processeur ecrit les informations   */

  if (flal != NULL) {

    if (cs_glob_rank_id <= 0) {

      if (_ipass == 1) {

        fprintf(flal,
                "# ** Information on Lagrangian computation\n"
                "#    --------------------------------------\n"
                "#\n"
                "# column  1: time step number\n"
                "# column  2: physical time\n"
                "# column  3: inst. number of particles\n"
                "# column  4: inst. number of particles (weighted)\n"
                "# column  5: inst. number of injected particles\n"
                "# column  6: inst. number of injected particles (weighted)\n"
                "# column  7: inst. number of exited, or deposited and removed particles\n"
                "# column  8: inst. number of exited, or deposited and removed particles (weighted)\n"
                "# column  9: inst. number of deposited particles\n"
                "# column 10: inst. number of deposited particles (weighted)\n");

        if (   lagr_model->physical_model == 2
            && lagr_model->fouling == 1)
          fprintf(flal,
                  "# column 11: inst. number of fouled particles (coal)\n"
                  "# column 12: inst. number of fouled particles (coal, weighted)\n"
                  "# column 13: inst. number of lost particles\n"
                  "#\n");

        else if (lagr_model->resuspension > 0)
          fprintf(flal,
                  "# column 11: inst. number of resuspended particles\n"
                  "# column 12: inst. number of resuspended particles (weighted)\n"
                  "# column 13: inst. number of lost particles\n"
                  "#\n");

        else
          fprintf(flal,
                  "# column 11: inst. number of lost particles\n"
                  "#\n");

      }

      /* ====================================================================
       * 2 - Ecriture des INFORMATIONS
       * ====================================================================   */

      if (   lagr_model->physical_model == 2
          && lagr_model->fouling == 1)
        fprintf(flal, " %8d %11.4E %8llu %11.4E %8llu %11.4E %8llu %11.4E %8llu %11.4E %8llu %11.4E %8llu\n",
                cs_glob_time_step->nt_cur,
                ttcabs,
                (unsigned long long)(pc->n_g_total), pc->w_total,
                (unsigned long long)(pc->n_g_new), pc->w_new,
                (unsigned long long)(pc->n_g_exit - pc->n_g_fouling),
                pc->w_exit - pc->w_fouling,
                (unsigned long long)(pc->n_g_deposited), pc->w_deposited,
                (unsigned long long)(pc->n_g_fouling), pc->w_fouling,
                (unsigned long long)(pc->n_g_failed));

      else if (lagr_model->resuspension > 0)
        fprintf(flal, " %8d %11.4E %8llu %11.4E %8llu %11.4E %8llu %11.4E %8llu %11.4E %8llu %11.4E %8llu\n",
                cs_glob_time_step->nt_cur,
                ttcabs,
                (unsigned long long)(pc->n_g_total), pc->w_total,
                (unsigned long long)(pc->n_g_new), pc->w_new,
                (unsigned long long)(pc->n_g_exit), pc->w_exit,
                (unsigned long long)(pc->n_g_deposited), pc->w_deposited,
                (unsigned long long)(pc->n_g_resuspended), pc->w_resuspended,
                (unsigned long long)(pc->n_g_failed));

      else
        fprintf(flal, " %8d %11.4E %8llu %11.4E %8llu %11.4E %8llu %11.4E %8llu %11.4E %8llu\n",
                cs_glob_time_step->nt_cur,
                ttcabs,
                (unsigned long long)(pc->n_g_total), pc->w_total,
                (unsigned long long)(pc->n_g_new), pc->w_new,
                (unsigned long long)(pc->n_g_exit),
                pc->w_exit,
                (unsigned long long)(pc->n_g_deposited), pc->w_deposited,
                (unsigned long long)(pc->n_g_failed));

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Close Lagrangian particle info file
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_print_finalize(void)
{
  if (flal != NULL) {
    fclose(flal);
    flal = NULL;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
