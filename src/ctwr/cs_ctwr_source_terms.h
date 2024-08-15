#ifndef __CS_CTWR_SOURCE_TERMS_H__
#define __CS_CTWR_SOURCE_TERMS_H__

/*============================================================================
 * Cooling towers related functions
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

#include "fvm_nodal.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS


/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change source terms - Exchange terms between the injected
 *        liquid and the water vapor phase in the bulk, humid air
 *
 * \param[in]     f_id          field id
 * \param[in,out] exp_st        Explicit source term
 * \param[in,out] imp_st        Implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_source_term(int              f_id,
                    cs_real_t        exp_st[],
                    cs_real_t        imp_st[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute volume mass injection for
 *   pressure (mass) equation.
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_volume_mass_injection_dof_func(cs_lnum_t         n_elts,
                                       const cs_lnum_t  *elt_ids,
                                       bool              dense_output,
                                       void             *input,
                                       cs_real_t        *retval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute volume mass injection for
 *   pressure (mass) equation for the rain.
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_volume_mass_injection_rain_dof_func(cs_lnum_t         n_elts,
                                            const cs_lnum_t  *elt_ids,
                                            bool              dense_output,
                                            void             *input,
                                            cs_real_t        *retval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute volume mass injection for
 *   yphp rain equation (enthalpy).
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_volume_mass_injection_yh_rain_dof_func(cs_lnum_t         n_elts,
                                               const cs_lnum_t  *elt_ids,
                                               bool              dense_output,
                                               void             *input,
                                               cs_real_t        *retval);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CTWR_SOURCE_TERMS_H__ */
