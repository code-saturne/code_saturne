#ifndef __CS_ZONE_H__
#define __CS_ZONE_H__

/*============================================================================
 * Base zones handling.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Zone description */

typedef struct {

  const char       *name;              /*!< zone name */

  int               id;                /*!< boundary zone id */
  int               type;              /*!< boundary zone type flag */

  int               location_id;       /*!< associated mesh location id */

  cs_lnum_t         n_elts;            /*!< local number of associated faces */
  const cs_lnum_t  *elt_ids;           /*!< associated face ids */

  bool              time_varying;      /*!< does the selected zone change
                                            with time ? */

  bool              allow_overlay;     /*!< allow overlaying of this zone ? */

  cs_real_t         measure;           /*!< Geometrical measure of the zone */
  cs_real_t         boundary_measure;  /*!< Boundary geometrical measure of the zone */

} cs_zone_t;

/*! Boundary zone description
  *
  * \deprecated  use cs_zone_t instead.
  */

typedef struct {

  const char       *name;           /*!< zone name */

  int               id;             /*!< boundary zone id */
  int               type;           /*!< boundary zone type flag */

  int               location_id;    /*!< associated mesh location id */

  cs_lnum_t         n_faces;        /*!< local number of associated faces */
  const cs_lnum_t  *face_ids;       /*!< associated face ids */


  bool              time_varying;   /*!< does the selected zone change
                                         with time ? */
  bool              allow_overlay;  /*!< allow overlaying of this zone ? */

  cs_real_t         surface;        /*!< Geometrical measure of the zone */
  cs_real_t         perimeter;      /*!< Boundary geometrical measure of the zone */

} cs_boundary_zone_t;

/*! Volume zone description
  *
  * \deprecated  use cs_zone_t instead.
  */

typedef struct {

  const char       *name;           /*!< zone name */

  int               id;             /*!< volume zone id */
  int               type;           /*!< volume zone type flag */

  int               location_id;    /*!< associated mesh location id */

  cs_lnum_t         n_cells;        /*!< local number of associated cells */
  const cs_lnum_t  *cell_ids;       /*!< associated cell ids */


  bool              time_varying;   /*!< does the selected zone change
                                      with time ? */
  bool              allow_overlay;  /*!< allow overlaying of this zone ? */

  cs_real_t         volume;         /*!< Geometrical measure of the zone */
  cs_real_t         surface;        /*!< Boundary geometrical measure of the zone */

} cs_volume_zone_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ZONE_H__ */
