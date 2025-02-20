#ifndef __CS_COMBUSTION_MODEL_H__
#define __CS_COMBUSTION_MODEL_H__

/*============================================================================
 * Combustion model parameters.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include <stdarg.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute Enthalpy and Cp based on the JANAF band.
 *
 * \param[in]   ncoel        number of elementary constituents
 * \param[in]   ngazem       number of elementary constituents
 * \param[in]   npo          number of interpolation points
 * \param[in]   nomcoel      names of elementary constituants
 * \param[out]  ehcoel       enthalpy for each elementary species
 *                           (for point i and species j, ehcoel[i*ngazem + j])
 * \param[out]  cpcoel       cp for each elementary species
 *                           (for point i and species j, cpcoel[i*ngazem + j])
 * \param[out]  coeff_therm  coefficients for the Burke-Scumann model,
 *                           or null otherwise
 * \param[in]   wmolce       molar mass of each species
 * \param[in]   th           temperature in K
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_enthalpy_and_cp_from_janaf(int           ncoel,
                                         int           ngazem,
                                         int           npo,
                                         const char    nomcoel[][13],
                                         double        ehcoel[],
                                         double        cpcoel[],
                                         double        coeff_therm[][2][5],
                                         const double  wmolce[],
                                         const double  th[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute rectangle-Dirac pdf parameters.
 *
 * From P. Plion & A. Escaich
 *
 * \param[in]       n_cells       number of cells
 * \param[out]      indpdf        indicator for pdf integration or mean value
 * \param[out]      tpdf          indicator for pdf shape:
 *                               - 0: Dirac at mean value
 *                               - 1: rectangle
 *                               - 2: Dirac's peak at \f$ f_{min} \f$
 *                               - 3: Dirac's peak at \f$ f_{max} \f$
 *                               - 4: rectangle and 2 Dirac's pics
 * \param[in]       fm            mean mixture fraction at cell centers
 * \param[in, out]  fp2m          mean mixture fraction variance at cell centers
 * \param[in, out]  fmini         mixture fraction low boundary
 * \param[in]       fmaxi         mixture fraction high boundary
 * \param[out]      dirmin        Dirac's peak value at \f$ f_{min} \f$
 * \param[out]      dirmax        Dirac's peak value at \f$ f_{max} \f$
 * \param[out]      fdeb          abscissa of rectangle low boundary
 * \param[out]      ffin          abscissa of rectangle high boundary
 * \param[out]      hrec          rectangle height
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_dirac_pdf(cs_lnum_t         n_cells,
                        int               indpdf[],
                        cs_real_t         tpdf[],
                        const cs_real_t   fm[],
                        cs_real_t         fp2m[],
                        const cs_real_t   fmini[],
                        const cs_real_t   fmaxi[],
                        cs_real_t         dirmin[],
                        cs_real_t         dirmax[],
                        cs_real_t         fdeb[],
                        cs_real_t         ffin[],
                        cs_real_t         hrec[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COMBUSTION_MODEL_H__ */
