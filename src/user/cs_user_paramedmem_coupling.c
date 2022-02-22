/*============================================================================
 * User functions for input of ParaMEDMEM coupling parameters
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
 * \file cs_user_paramedmem_coupling.c
 *
 * \brief User functions for input of ParaMEDMEM coupling parameters
 *
 * \brief User functions for input of calculation parameters.
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define ParaMEDMEM coupling(s)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_paramedmem_define_couplings
void
cs_user_paramedmem_define_couplings(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupled meshes
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_paramedmem_define_meshes
void
cs_user_paramedmem_define_meshes(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define fields to couple with ParaMEDMEM
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_paramedmem_define_fields
void
cs_user_paramedmem_define_fields(void)
{


}

END_C_DECLS
