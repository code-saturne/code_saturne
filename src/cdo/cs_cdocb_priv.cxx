/*============================================================================
 * Structure and functions common to all CDO cell-based schemes but not public
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_param_saddle.h"
#include "cs_param_sles.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdocb_priv.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

BEGIN_C_DECLS

/*!
  \file cs_cdocb_priv.c

  \brief Common functions for CDO cell-based schemes which are not visible by
         default from user-defined functions

*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the default settings for a scalar-valued CDO cell-based scheme
 *
 * \param[in, out] eqp  set of equation parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_init_default_param(cs_equation_param_t  *eqp)
{
  if (eqp == nullptr)
    return;

  /* Discretization parameters */
  /* ------------------------- */

  eqp->space_scheme = CS_SPACE_SCHEME_CDOCB;
  eqp->space_poly_degree = 0;

  eqp->diffusion_hodgep.inv_pty = true;
  eqp->diffusion_hodgep.type = CS_HODGE_TYPE_FPED;
  eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_COST;

  eqp->time_hodgep.type = CS_HODGE_TYPE_VDCP;
  eqp->time_hodgep.algo = CS_HODGE_ALGO_VORONOI;

  eqp->reaction_hodgep.type = CS_HODGE_TYPE_VDCP;
  eqp->reaction_hodgep.algo = CS_HODGE_ALGO_VORONOI;

  /* Linear algebra
   * -------------- *
   *
   * This is a saddle-point system. Define a default settings
   */

  /* Set the name of the saddle-point problem to solve */

  char *saddle_name = nullptr;
  int  len = strlen(eqp->name) + strlen("_saddle_point");

  BFT_MALLOC(saddle_name, len + 1, char);
  sprintf(saddle_name, "%s_saddle_point", eqp->name);

  cs_param_saddle_set_name(saddle_name, eqp->saddle_param);

  BFT_FREE(saddle_name);

  /* Associate the cs_param_sles_t structure related to the (1,1)-block */

  cs_param_saddle_t  *saddlep = cs_equation_param_get_saddle_param(eqp);
  cs_param_sles_t  *b11_slesp = cs_equation_param_get_sles_param(eqp);

  cs_param_saddle_set_block11_sles_param(saddlep, b11_slesp);

  /* Default saddle-point solver */

#if defined(HAVE_MUMPS)
  cs_equation_param_set(eqp, CS_EQKEY_SADDLE_SOLVER, "alu");
  cs_equation_param_set(eqp, CS_EQKEY_SADDLE_AUGMENT_SCALING, "100");
  cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "mumps");

  /* Advanced settings */

  cs_param_sles_mumps(b11_slesp,
                      false, /* is single-precision */
                      CS_PARAM_MUMPS_FACTO_LDLT_SYM);
#else
    cs_equation_param_set(eqp, CS_EQKEY_SADDLE_SOLVER, "gcr");
    cs_equation_param_set(eqp, CS_EQKEY_SADDLE_PRECOND, "sgs");
    cs_equation_param_set(eqp, CS_EQKEY_SADDLE_SCHUR_APPROX, "mass_scaled");
    cs_equation_param_set(eqp, CS_EQKEY_SADDLE_SOLVER_RESTART, "40");
    cs_equation_param_set(eqp, CS_EQKEY_SADDLE_MAX_ITER, "100");

    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER, "fcg");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_MAX_ITER, "50");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_RTOL, "1e-1");
#endif
}

/*----------------------------------------------------------------------------*/
END_C_DECLS
