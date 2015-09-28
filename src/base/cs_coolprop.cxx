/*============================================================================
 * Equation of state
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * CoolProp library headers
 *----------------------------------------------------------------------------*/

#include "CoolProp.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coolprop.hxx"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes physical properties in (P,h,Yi) for compressible flow.
 *
 * \param[in]  CoolPropMaterial type of material
 * \param[in]  thermo_plane     type of thermal plane
 * \param[in]  property         type of property to compute
 * \param[in]  n_vals           size of properties arrays
 * \param[in]  var1             array of pressure
 * \param[in]  var2             array of thermal properties
 * \param[out] val              array of property
 */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
void
cs_phys_prop_coolprop(char                              *CoolPropMaterial,
                      cs_phys_prop_thermo_plane_type_t   thermo_plane,
                      cs_phys_prop_type_t                property,
                      const cs_lnum_t                    n_vals,
                      double                             var1[],
                      double                             var2[],
                      cs_real_t                          val[])
{
  std::vector<std::string> fluids;
  std::vector<std::string> outputs;

  fluids.push_back(CoolPropMaterial);
  char *variable1;
  BFT_MALLOC(variable1, 2, char);
  char variable2[] = "P";

  if (thermo_plane == CS_PHYS_PROP_PLANE_PH)
    strcpy(variable1, "H");
  else if (thermo_plane == CS_PHYS_PROP_PLANE_PT)
    strcpy(variable1, "T");

  switch (property) {
    case CS_PHYS_PROP_PRESSURE:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: you can't choose pressure as an output variable"));
      break;
    case CS_PHYS_PROP_TEMPERATURE:
      /* temperature is in K */
      outputs.push_back("T");
      break;
    case CS_PHYS_PROP_ENTHALPY:
      outputs.push_back("H");
      break;
    case CS_PHYS_PROP_ENTROPY:
      outputs.push_back("S");
      break;
    case CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY:
      outputs.push_back("C");
      break;
    case CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: isochoric heat capacity not available yet"));
      break;
    case CS_PHYS_PROP_SPECIFIC_VOLUME:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: specific volume not available yet"));
      break;
    case CS_PHYS_PROP_DENSITY:
      outputs.push_back("D");
      break;
    case CS_PHYS_PROP_INTERNAL_ENERGY:
      outputs.push_back("U");
      break;
    case CS_PHYS_PROP_QUALITY:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: quality not available yet"));
      break;
    case CS_PHYS_PROP_THERMAL_CONDUCTIVITY:
      outputs.push_back("L");
      break;
    case CS_PHYS_PROP_DYNAMIC_VISCOSITY:
      outputs.push_back("V");
      break;
    case CS_PHYS_PROP_SPEED_OF_SOUND:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: speed of sound not available yet"));
      break;
  }

  //std::vector<double> T(1,298), p(1e5);
  //std::cout << PropsSImulti(outputs,"T", T, "P", p, "", fluids, z)[0][0] << std::endl;
  //std::cout << PropsSImulti(outputs,"T", T, "P", p, "HEOS", fluids, z)[0][0] << std::endl;
  //std::cout << PropsSImulti(outputs,"T", T, "P", p, "REFPROP", fluids, z)[0][0] << std::endl;

  std::vector<double> val1;
  std::vector<double> val2;
  std::vector<double> fractions;

  val1.clear();
  val2.clear();
  fractions.clear();

  for (int i = 0; i < n_vals; i++)
  {
    val1.push_back(var2[i]);
    val2.push_back(var1[i]);
  }
  fractions.push_back(1.0);

  std::vector<std::vector<double> > out = CoolProp::PropsSImulti(outputs,
                                                                 variable1,
                                                                 val1,
                                                                 variable2,
                                                                 val2,
                                                                 "HEOS",
                                                                 fluids,
                                                                 fractions);

  for (int i = 0; i < n_vals; i++)
    val[i] = out[i][0];

  BFT_FREE(variable1);
}

/*----------------------------------------------------------------------------*/



