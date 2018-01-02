#ifndef __CS_RANDOM_H__
#define __CS_RANDOM_H__

/*============================================================================
 * Random number generation.
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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize random number generator.
 *
 * Generates initial seed buffer by linear congruential method.
 * Taken from Marsaglia, FSU report FSU-SCRI-87-50.
 *
 * \param[in]  seed  variable seed, with 0 < seed < 31328
 */
/*----------------------------------------------------------------------------*/

void
cs_random_seed(int  seed);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Uniform distribution random number generator.
 *
 * Portable lagged Fibonacci series uniform random number generator
 * with "lags" -273 und -607:
 * W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92
 *
 * \param[in]   n  number of values to compute
 * \param[out]  a  pseudo-random numbers following uniform distribution
 */
/*----------------------------------------------------------------------------*/

void
cs_random_uniform(cs_lnum_t  n,
                  cs_real_t  a[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Normal distribution random number generator.
 *
 * Box-Muller method for Gaussian random numbers.
 *
 * \param[in]   n  number of values to compute
 * \param[out]  x  pseudo-random numbers following normal distribution
 */
/*----------------------------------------------------------------------------*/

void
cs_random_normal(cs_lnum_t  n,
                 cs_real_t  x[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Poisson distribution random number generator.
 *
 * q(mu,p) = exp(-mu) mu**p/p!
 *
 * \param[in]   n   number of values to compute
 * \param[in]   mu  Poisson distribution parameter
 * \param[out]  p   pseudo-random numbers following Poisson distribution
 */
/*----------------------------------------------------------------------------*/

void
cs_random_poisson(cs_lnum_t  n,
                  cs_real_t  mu,
                  int        p[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Save static variables used by random number generator.
 *
 * \param[out]  save_block  saved state values
 */
/*----------------------------------------------------------------------------*/

void
cs_random_save(cs_real_t  save_block[1634]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Restore static variables used by random number generator.
 *
 * \param[out]  save_block  saved state values
 */
/*----------------------------------------------------------------------------*/

void
cs_random_restore(cs_real_t  save_block[1634]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RANDOM_H__ */
