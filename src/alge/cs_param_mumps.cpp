/*============================================================================
 * Routines and structure to handle the MUMPS settings
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

#include "cs_base.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_param_mumps.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize with the default settings a new structure
 *        storing a set of parameters used when calling MUMPS.
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_param_mumps_t *
cs_param_mumps_create(void)
{
  cs_param_mumps_t  *mumpsp = nullptr;

  BFT_MALLOC(mumpsp, 1, cs_param_mumps_t);

  mumpsp->is_single = false;
  mumpsp->facto_type = CS_PARAM_MUMPS_FACTO_LU;

  /* Advanced options */

  mumpsp->analysis_algo = CS_PARAM_MUMPS_ANALYSIS_AUTO;
  mumpsp->mem_usage = CS_PARAM_MUMPS_MEMORY_AUTO;

  mumpsp->advanced_optim = false;   /* No advanced MPI/OpenMP optimization */
  mumpsp->blr_threshold = 0;        /* No BLR */
  mumpsp->mem_coef = -1;            /* No additional memory range */
  mumpsp->block_analysis = 0;       /* No clustered analysis */
  mumpsp->ir_steps = 0;             /* No iterative refinement */

  return mumpsp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy into a new structure the given set of parameters used when
 *        calling MUMPS
 *
 * \param[in] mumpsp   set of mumps parameters
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_param_mumps_t *
cs_param_mumps_copy(const cs_param_mumps_t  *mumpsp)
{
  cs_param_mumps_t  *cpy = cs_param_mumps_create();

  cpy->analysis_algo = mumpsp->analysis_algo;
  cpy->facto_type = mumpsp->facto_type;
  cpy->mem_usage = mumpsp->mem_usage;

  cpy->is_single = mumpsp->is_single;
  cpy->advanced_optim = mumpsp->advanced_optim;
  cpy->blr_threshold = mumpsp->blr_threshold;
  cpy->mem_coef = mumpsp->mem_coef;
  cpy->block_analysis = mumpsp->block_analysis;
  cpy->ir_steps = mumpsp->ir_steps;

  return cpy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the structure storing the set of parameters used with MUMPS
 *
 * \param[in] name     name related to the current SLES
 * \param[in] mumpsp   set of mumps parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_param_mumps_log(const char              *name,
                   const cs_param_mumps_t  *mumpsp)
{
  if (mumpsp == nullptr)
    return;

  char type = (mumpsp->is_single) ? 's' : 'd';
  char tag[32];

  switch (mumpsp->facto_type) {

  case CS_PARAM_MUMPS_FACTO_LU:
    sprintf(tag, "%cmumps_lu", type);
    break;
  case CS_PARAM_MUMPS_FACTO_LDLT_SYM:
    sprintf(tag, "%cmumps_ldlt_sym", type);
    break;
  case CS_PARAM_MUMPS_FACTO_LDLT_SPD:
    sprintf(tag, "%cmumps_ldlt_spd", type);
    break;

  default:
    sprintf(tag, "undefined");
    break;

  }

  cs_log_printf(CS_LOG_SETUP, "  * %s | MUMPS_type:               %s\n",
                name, tag);

  /* Strategy for the memory usage */

  switch (mumpsp->mem_usage) {
  case CS_PARAM_MUMPS_MEMORY_CONSTRAINED:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Memory_usage:             %s\n",
                  name, "constrained");
    break;
  case CS_PARAM_MUMPS_MEMORY_AUTO:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Memory_usage:             %s\n",
                  name, "automatic");
    break;
  case CS_PARAM_MUMPS_MEMORY_CPU_DRIVEN:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Memory_usage:             %s\n",
                  name, "CPU-driven (efficiency first)");
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Memory_usage:             %s\n",
                  name, "Undefined");
    break;
  }

  /* Algorithm used for the analysis step */

  switch (mumpsp->analysis_algo) {
  case CS_PARAM_MUMPS_ANALYSIS_AMD:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:            %s\n",
                  name, "AMD");
    break;
  case CS_PARAM_MUMPS_ANALYSIS_QAMD:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:            %s\n",
                  name, "QAMD");
    break;
  case CS_PARAM_MUMPS_ANALYSIS_PORD:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:            %s\n",
                  name, "PORD");
    break;
  case CS_PARAM_MUMPS_ANALYSIS_SCOTCH:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:            %s\n",
                  name, "SCOTCH");
    break;
  case CS_PARAM_MUMPS_ANALYSIS_PTSCOTCH:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:            %s\n",
                  name, "PT-SCOTCH");
    break;
  case CS_PARAM_MUMPS_ANALYSIS_METIS:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:            %s\n",
                  name, "METIS");
    break;
  case CS_PARAM_MUMPS_ANALYSIS_PARMETIS:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:            %s\n",
                  name, "PARMETIS");
    break;
  case CS_PARAM_MUMPS_ANALYSIS_AUTO:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:            %s\n",
                  name, "automatic choice done by MUMPS");
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "  * %s | Analysis_algo:            %s\n",
                  name, "Undefined");
    break;
  }

  cs_log_printf(CS_LOG_SETUP, "  * %s | Advanced_Optim:           %s\n",
                name, cs_base_strtf(mumpsp->advanced_optim));

  if (mumpsp->block_analysis > 1)
    cs_log_printf(CS_LOG_SETUP, "  * %s | Block_Size in analysis:   %d\n",
                  name, mumpsp->block_analysis);

  if (mumpsp->ir_steps != 0)
    cs_log_printf(CS_LOG_SETUP, "  * %s | Iterative_Refinement:      %d\n",
                  name, CS_ABS(mumpsp->ir_steps));

  if (fabs(mumpsp->blr_threshold) > FLT_MIN)
    cs_log_printf(CS_LOG_SETUP, "  * %s | BLR_threshold:             %e\n",
                  name, mumpsp->blr_threshold);

  if (mumpsp->mem_coef > 0)
    cs_log_printf(CS_LOG_SETUP, "  * %s | Memory pct. increase:     %f\n",
                  name, mumpsp->mem_coef);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
