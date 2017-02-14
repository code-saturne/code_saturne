#ifndef __CS_BALANCE_BY_ZONE_H__
#define __CS_BALANCE_BY_ZONE_H__

/*============================================================================
 * Scalar balance on zones.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes the different terms of the balance of a scalar which name is
 * given as argument, on a volumic zone defined by the criterion also given as
 * argument. The different contributions to the balance are printed in the
 * listing.
 *
 * \param[in]     selection_crit      zone selection criterion
 * \param[in]     scalar_name         scalar name
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_by_zone(const char  *selection_crit,
                   const char  *scalar_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes one term of the head loss balance (pressure drop) on a
 * volumic zone defined by the criterion also given as argument.
 * The different contributions are printed in the listing.
 *
 * \param[in]     selection_crit      zone selection criterion
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_drop_by_zone(const char  *selection_crit);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes the surface balance of a scalar which name is
 * given as argument, on a oriented surface area defined by the criterion also
 * given as argument. The flux is counted negatively in the given
 * outwards-facing normal direction.
 *
 *
 * \param[in]     selection_crit      zone selection criterion
 * \param[in]     scalar_name         scalar name
 * \param[in]     normal              out normal surface
 */
/*----------------------------------------------------------------------------*/

void
cs_surface_balance(const char       *selection_crit,
                   const char       *scalar_name,
                   const cs_real_t   normal[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the face by face surface flux of a scalar which name is
 * given as argument, on a surface area defined by the criterion also given as
 * argument. Counted negatively through the normal.
 *
 * \param[in]     selection_crit      zone selection criterion
 * \param[in]     scalar_name         scalar name
 * \param[in]     normal              normal surface
 * \param[in,out] flux_b_faces        pointer surface flux boundary faces
 * \param[in,out] flux_i_faces        pointer surface flux internal faces
 * \param[out]    n_b_faces_sel       number of boundary faces
 * \param[out]    n_i_faces_sel       number of internal faces
 * \param[out]    b_face_sel_ids      ids of boundary faces
 * \param[out]    i_face_sel_ids      ids of internal faces
 */
/*----------------------------------------------------------------------------*/

void
cs_flux_through_surface(const char          *selection_crit,
                        const char          *scalar_name,
                        const cs_real_t      normal[3],
                        cs_real_t          *flux_b_faces,
                        cs_real_t          *flux_i_faces,
                        cs_lnum_t          *n_b_faces_sel,
                        cs_lnum_t          *n_i_faces_sel,
                        cs_lnum_t          *b_face_sel_ids,
                        cs_lnum_t          *i_face_sel_ids);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BALANCE_BY_ZONE_H__ */
