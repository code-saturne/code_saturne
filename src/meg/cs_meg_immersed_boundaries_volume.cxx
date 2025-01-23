/*----------------------------------------------------------------------------*/

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

#include "base/cs_defs.h"

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_meg_immersed_boundaries_volume.cxx
 *
 * \brief This function is used to indicate to compute user defined values
 *        over a given immersed object. The mathematical expression is defined
 *        in the GUI.
 *
 * \param[in] object_name      Name of the immersed object
 * \param[in] var_name         Name of the variable
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_meg_ibm_volume_func_by_name
cs_ibm_volume_func_t *
cs_meg_ibm_volume_func_by_name(const char *object_name,
                               const char *gui_var_name)
{
  CS_NO_WARN_IF_UNUSED(object_name);
  CS_NO_WARN_IF_UNUSED(gui_var_name);

  cs_ibm_volume_func_t *f = nullptr;

  return f;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
