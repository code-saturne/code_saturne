/*============================================================================
 * Routines to handle the definition and usage of material properties
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_base.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_param_cdo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_PARAM_CDO_DBG  0

/*============================================================================
 * Global variables
 *============================================================================*/

/* Separation lines: long, medium, short */
const char lsepline[80] =
  "# =======================================================================\n";
const char msepline[60] =
  "# =========================================\n";
const char ssepline[40] =
  "# =================\n";

/*============================================================================
 * Global static variables
 *============================================================================*/

static const char
cs_param_hodge_type_desc[CS_PARAM_N_HODGE_TYPES][CS_BASE_STRING_LEN] =
  { N_("VpCd"),
    N_("EpFd"),
    N_("FpEd"),
    N_("EdFp"),
    N_("CpVd")  };

static const char
cs_param_hodge_algo_desc[CS_PARAM_N_HODGE_ALGOS][CS_BASE_STRING_LEN] =
  { N_("Voronoi"),
    N_("Whitney on the Barycentric Subdivision (WBS)"),
    N_("COnsistency-STabilization splitting (COST)"),
    N_("Automatic switch") };

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of algorithm related to a discrete Hdoge operator
 *
 * \param[in] h_info     cs_param_hodge_t structure
 *
 * \return the name of the algorithm
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_hodge_get_algo_name(const cs_param_hodge_t   h_info)
{
  return cs_param_hodge_algo_desc[h_info.algo];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the type of discrete Hodge operator
 *
 * \param[in] h_info     cs_param_hodge_t structure
 *
 * \return the name of the type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_hodge_get_type_name(const cs_param_hodge_t   h_info)
{
  return cs_param_hodge_type_desc[h_info.type];
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
