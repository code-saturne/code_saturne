#ifndef __CS_BALANCE_BY_ZONE_H__
#define __CS_BALANCE_BY_ZONE_H__

/*============================================================================
 * Scalar balance on zones.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 * given as argument, on a volumic zone defined by the criterium also given as
 * argument. The different contributions to the balance are printed in the
 * listing.
 *
 * \param[in]     selection_crit      zone selection criterium
 * \param[in]     scalar_name         scalar name
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_by_zone(const char *selection_crit,
                   const char *scalar_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes one term of the head loss balance (pressure drop) on a
 * volumic zone defined by the criterium also given as argument.
 * The different contributions are printed in the listing.
 *
 * \param[in]     selection_crit      zone selection criterium
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_drop_by_zone(const char *selection_crit);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BALANCE_BY_ZONE_H__ */
