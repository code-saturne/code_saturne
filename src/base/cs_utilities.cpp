/*============================================================================
 * Global utilities functions used for initialize and finalize calls.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_medcoupling_intersector.h"
#include "cs_medcoupling_postprocess.h"
#include "cs_medcoupling_remapper.h"
#include "cs_medcoupling_utils.h"
#include "cs_stl.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_utilities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy different intersectors/remapping structures related to STL or
 * MEDCoupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_utilities_destroy_all_remapping(void)
{
  cs_medcoupling_intersector_destroy_all();
  cs_medcoupling_remapper_destroy_all();
  cs_medcoupling_slice_destroy_all();

  cs_medcoupling_free_meshes();
  cs_stl_mesh_destroy_all();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
