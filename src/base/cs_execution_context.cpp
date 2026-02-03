/*============================================================================
 * Class to handle different execution policies (MPI, OpenMP, CUDA, ...)
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
#include "base/cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_execution_context.h"

/*============================================================================
 * Static global variables
 *============================================================================*/


static cs::execution::environment *_default_env = nullptr;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* Get the current execution context, for the moment global is returned */
/*----------------------------------------------------------------------------*/

const cs::execution::environment *
cs::execution::default_env(void)
{
  return _default_env;
}

/*----------------------------------------------------------------------------*/
/* Get the global dispatch context. */
/*----------------------------------------------------------------------------*/

cs_dispatch_context&
cs::execution::default_context(void)
{
  return _default_env->g_ctx;
}

/*----------------------------------------------------------------------------*/
/* Get the global host context. */
/*----------------------------------------------------------------------------*/

cs_host_context&
cs::execution::default_h_context(void)
{
  cs_host_context h_ctx = static_cast<cs_host_context>(_default_env->g_ctx);
  return h_ctx;
}

cs::execution::mpi_wrapper&
cs::execution::default_mpi(void)
{
  return _default_env->mpi;
}
/*----------------------------------------------------------------------------*/
/* Initialize the global execution context. */
/*----------------------------------------------------------------------------*/

void
cs_execution_default_env_init(void)
{
  _default_env = new cs::execution::environment();
#if defined(HAVE_MPI)
  _default_env->mpi.set_comm(cs_glob_mpi_comm);
#endif
}

/*----------------------------------------------------------------------------*/
/* Free the global execution context pointer. */
/*----------------------------------------------------------------------------*/

void
cs_execution_default_env_finalize(void)
{
  if (_default_env != nullptr)
    delete _default_env;
}

/*----------------------------------------------------------------------------*/
