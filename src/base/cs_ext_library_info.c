/*============================================================================
 * External library information.
 *============================================================================*/

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
#include <string.h>
#include <time.h>

#if defined(__linux__)
#include <stdio.h>
#endif

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined HAVE_SYS_UTSNAME_H
#include <sys/utsname.h>
#endif

#if defined(HAVE_SYS_SYSINFO_H) && defined(HAVE_SYSINFO)
#include <sys/sysinfo.h>
#endif

#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)
#include <pwd.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "cs_log.h"

#if defined(HAVE_CUDA)
#include "cs_base_cuda.h"
#endif

#include "cs_partition.h"

#if defined(HAVE_PETSC)
#if 0
#include "cs_sles_petsc.h"
#else
/* Duplicate prototype here to avoid requiring PETSc headers */
void
cs_sles_petsc_library_info(cs_log_t  log_type);
#endif
#endif

#if defined(HAVE_HYPRE)
#include "cs_sles_hypre.h"
#endif

#if defined(HAVE_AMGX)
#include "cs_sles_amgx.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ext_library_info.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print external library info.
 *
 * This may be incomplete.
 *
 * \param[in]  log   if true, standard logging; otherwise, single output
 */
/*----------------------------------------------------------------------------*/

static void
_ext_library_version_info(bool  log)
{
  int  n_logs = (log) ? 2 : 1;
  cs_log_t logs[] = {CS_LOG_DEFAULT, CS_LOG_PERFORMANCE};

  int n_ext = 0;

#if defined(HAVE_PETSC)
  n_ext += 1;
#endif

#if defined(HAVE_HYPRE)
  n_ext += 1;
#endif

#if defined(HAVE_AMGX)
  n_ext += 1;
#endif

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)
  n_ext += 1;
#endif
#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)
  n_ext += 1;
#endif

  if (n_ext < 1)
    return;

  for (int log_id = 0; log_id < n_logs; log_id++) {

    cs_log_printf(logs[log_id],
                  "\n  External libraries:\n");

#if defined(HAVE_PETSC)
    cs_sles_petsc_library_info(log_id);
#endif
#if defined(HAVE_HYPRE)
    cs_sles_hypre_library_info(log_id);
#endif
#if defined(HAVE_AMGX)
    cs_sles_amgx_library_info(log_id);
#endif

#if    defined(HAVE_METIS)  || defined(HAVE_PARMETIS) \
    || defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)
    cs_partition_external_library_info(log_id);
#endif

  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print available system information.
 *
 * \param[in]  comm  associated MPI communicator
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_library_info(void)
{
  _ext_library_version_info(true);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print available system information, without additional logging
 *
 * \param[in]  comm  associated MPI communicator
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_library_info_no_log(void)
{
  _ext_library_version_info(false);
}

/*-----------------------------------------------------------------------------*/

END_C_DECLS
