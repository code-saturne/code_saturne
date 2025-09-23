#ifndef __CS_ATMO_IMBRICATION_H__
#define __CS_ATMO_IMBRICATION_H__

/*============================================================================
 * Main for atmospheric imbrication related functions
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------
 * Atmospheric model options descriptor
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Atmospheric imbrication option
 *----------------------------------------------------------------------------*/

typedef struct {

  /*! activation flag */
  bool imbrication_flag;
  bool imbrication_verbose;

  /*! Flags for activating the cressman interpolation for the boundary
      conditions */
  bool cressman_u;
  bool cressman_v;
  bool cressman_qw;
  bool cressman_nc;
  bool cressman_tke;
  bool cressman_eps;
  bool cressman_theta;

  /*! numerical parameters for the cressman interpolation formulas */
  cs_real_t vertical_influence_radius;
  cs_real_t horizontal_influence_radius;

  /*! additional variables */
  int id_u;
  int id_v;
  int id_qw;
  int id_nc;
  int id_tke;
  int id_eps;
  int id_theta;

} cs_atmo_imbrication_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to atmo imbrication structure */
extern cs_atmo_imbrication_t *cs_glob_atmo_imbrication;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*! ----------------------------------------------------------------
 * \brief  Final step for free arrays imbrication
 * ---------------------------------------------------------------- */

void
cs_finalize_imbrication(void);

/*----------------------------------------------------------------------------*/

/*!----------------------------------------------------------------------------
 * \brief Prepare data for imbrication by reading meteo files
 *
 * Warning : the list of files is supposed to be "imbrication_files_list.txt"
 * --------------------------------------------------------------------------- */

void
cs_activate_imbrication(void);

/*----------------------------------------------------------------------------*/

/*!----------------------------------------------------------------------------
 * \brief Prepare for the cressman interpolation of the variables
 *
 * \param[in]  the_time        current time
 * --------------------------------------------------------------------------- */

void
cs_summon_cressman(cs_real_t the_time);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ATMO_IMBRICATION_H__ */
