/*============================================================================
 * User functions for input of calculation parameters.
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
 * \file cs_user_parameters-ctwr.c
 *
 * \brief Cooling towers parameters example
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate cooling tower model.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{
  /* Activate cooling tower model */

  /*! [ctwr_user_model_1] */
  cs_glob_physical_model_flag[CS_COOLING_TOWERS] = 0;
  /*! [ctwr_user_model_1] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /*
   * We define a cooling tower zone
   */

  /*! [ctwr_user_1] */
  {

    cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();

    /* Evaporation model:
       CS_CTWR_NONE None,
       CS_CTWR_POPPE Poppe,
       CS_CTWR_MERKEL Merkel*/
    ct_opt->evap_model = CS_CTWR_POPPE;

    cs_real_t surface = 0.48 * 6540.; /* 48% of the total disc */
    cs_real_t qw = surface *  2.64; /* Water flow rate (kg/s) */

    cs_ctwr_define(
        "2 or 3", /* selction criterion */
        CS_CTWR_COUNTER_CURRENT, /*Type:
                                   CS_CTWR_COUNTER_CURRENT counter current,
                                   CS_CTWR_CROSS_CURRENT cross,
                                   CS_CTWR_RAIN rain zone */
        -1., /* Imposed delta temperature if positive */
        0.1, /* Associated relaxation time */
        36., /* Liquid injected water temperature */
        qw,
        0.2, /* Evaportaion law constant A */
        0.5, /* Evaportaion law constant n */
        surface,
        -1.); /* Leaking factor, not taken into account if negative */

  }
  /*! [ctwr_user_1] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
