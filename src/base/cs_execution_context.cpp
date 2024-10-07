/*============================================================================
 * Class to handle different execution policies (MPI, OpenMP, CUDA, ...)
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_execution_context.h"

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_execution_context *_glob_context = nullptr;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* Get the current execution context, for the moment global is returned */
/*----------------------------------------------------------------------------*/

const cs_execution_context *
cs_execution_context_get(void)
{
  return _glob_context;
}

/*----------------------------------------------------------------------------*/
/* Get the global execution context. */
/*----------------------------------------------------------------------------*/

const cs_execution_context *
cs_execution_context_glob_get(void)
{
  return _glob_context;
}

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/* Initialize the global execution context. */
/*----------------------------------------------------------------------------*/

void
cs_execution_context_glob_init(void)
{
  _glob_context = new cs_execution_context();
#if defined(HAVE_MPI)
  _glob_context->set_comm(cs_glob_mpi_comm);
#endif
}

/*----------------------------------------------------------------------------*/
/* Free the global execution context pointer. */
/*----------------------------------------------------------------------------*/

void
cs_execution_context_glob_finalize(void)
{
  if (_glob_context != nullptr)
    delete _glob_context;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
