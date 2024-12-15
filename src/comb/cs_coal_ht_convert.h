#ifndef CS_COAL_HT_CONVERT_H
#define CS_COAL_HT_CONVERT_H

/*============================================================================
 * Coal combustion model: enthaly to and from temperature conversion.
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
#include "cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Calculation of the gas temperature from gas enthalpy and
 *        concentrations at cells for coal combustion.
 *
 * \param[in]       location_id     mesh location id (cells or boundary faces)
 * \param[in]       eh              gas enthalpy (j/kg of gaseous mixture)
 * \param[in, out]  tp              gas temperature in kelvin
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_ht_convert_h_to_t_gas(int              location_id,
                              const cs_real_t  eh[],
                              cs_real_t        tp[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Calculation of the gas temperature from gas enthalpy and
 *        given mass fractions and average f1/f2 for coal combustion.
 *
 * \param[in]  eh            gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 * \param[in]  xesp          mass fraction (yi) of species
 * \param[in]  f1mc          average f1 per coal
 * \param[in]  f2mc          average f2 per coal
 *
 * \return  gas temperature (in kelvin)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_h_to_t_gas_by_yi_f1f2(cs_real_t        eh,
                                         const cs_real_t  xesp[],
                                         const cs_real_t  f1mc[],
                                         const cs_real_t  f2mc[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Calculation of the gas enthalpy from gas temperature and
 *        given mass fractions and average f1/f2 for coal combustion.
 *
 * \param[in]  tp            gas temperature (in kelvin)
 * \param[in]  xesp          mass fraction (yi) of species
 * \param[in]  f1mc          average f1 per coal
 * \param[in]  f2mc          average f2 per coal
 *
 * \return  gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_t_to_h_gas_by_yi_f1f2(cs_real_t        tp,
                                         const cs_real_t  xesp[],
                                         const cs_real_t  f1mc[],
                                         const cs_real_t  f2mc[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Calculation of the gas temperature from gas enthalpy and
 *        given mass fractions for coal combustion.
 *
 * \param[in]  eh            gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 * \param[in]  xesp          mass fraction (yi) of species
 * \param[in]  f1mc          average f1 per coal
 * \param[in]  f2mc          average f2 per coal
 *
 * \return  gas temperature (in kelvin)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_h_to_t_gas_by_yi(cs_real_t        eh,
                                    const cs_real_t  xesp[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the gas enthalpy from gas temperature and
 *        given mass fractions for coal combustion and 0 f1 and f2 values
 *
 * \param[in]  tp            gas temperature (in kelvin)
 * \param[in]  xesp          mass fraction (yi) of species
 *
 * \return  gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_t_to_h_gas_by_yi(cs_real_t        tp,
                                    const cs_real_t  xesp[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Calculation of the gas temperature from gas enthalpy and
 *        given mass fractions for coal combustion with drying.
 *
 * \param[in]  eh            gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 * \param[in]  xesp          mass fraction (yi) of species
 *
 * \return  gas temperature (in kelvin)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_h_to_t_gas_by_yi_with_drying(cs_real_t        eh,
                                                const cs_real_t  xesp[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Calculation of the gas enthalpy from gas temperature and
 *        given mass fractions for coal combustion with drying.
 *
 * \param[in]  tp            gas temperature (in kelvin)
 * \param[in]  xesp          mass fraction (yi) of species
 *
 * \return  gas enthalpy (\f$ j . kg^{-1} \f$ of mixed gas)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_t_to_h_gas_by_yi_with_drying(cs_real_t        tp,
                                                const cs_real_t  xesp[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Calculation of the particles temperature from particles enthalpy and
 *        concentrations at cells for coal combustion.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_ht_convert_h_to_t_particles(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Calculation of the particles temperature from particles enthalpy and
 *        given mass fractions for coal combustion.
 *
 * \remark  Function not called in code, so should probably be removed,
 *          unless useful for advanced postprocessing.
 *
 * \param[in]  enthal     mass enthalpy (\f$ j . kg^{-1} \f$)
 * \param[in]  class_id   class id (0 to n-1)
 * \param[in]  xesp       mass fraction of components
 *                        (size: cm->nsolid)
 * \param[in]  t1         coal inlet/boundary temperature
 *
 * \return   temperature (in kelvin)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_h_to_t_particles_by_yi(cs_real_t        enthal,
                                          int              class_id,
                                          const cs_real_t  xsolid[],
                                          cs_real_t        t1);

/*----------------------------------------------------------------------------*/
/*
 * \brief Calculation of the particles from particles temperature and
 *        given mass fractions for coal combustion.
 *
 * \param[in]  temper        temperature (in kelvin)
 * \param[in]  class_id      class id (0 to n-1)
 * \param[in]  xesp          mass fraction of components
 *                           (size: cm->nsolid)
 *
 * \return  mass enthalpy (\f$ j . kg^{-1} \f$)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_coal_ht_convert_t_to_h_particles_by_yi(cs_real_t        temper,
                                          int              class_id,
                                          const cs_real_t  xsolid[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Convert temperature to enthalpy at boundary for coal combustion.
 *
 * \param[in]   n_faces   number of faces in list
 * \param[in]   face_ids  list of boundary faces at which conversion
 *                        is requested (0-based numbering)
 * \param[in]   t_b       temperature at boundary
 * \param[out]  h_b       enthalpy at boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_ht_convert_t_to_h_faces(cs_lnum_t        n_faces,
                                const cs_lnum_t  face_ids[],
                                const cs_real_t  t_b[],
                                cs_real_t        h_b[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_COAL_HT_CONVERT_H */
