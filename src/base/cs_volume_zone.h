#ifndef __CS_VOLUME_ZONE_H__
#define __CS_VOLUME_ZONE_H__

/*============================================================================
 * Volume zones handling.
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

/*!
 * @defgroup volume_zone_flags Flags specifying general volume zone type
 *
 * @{
 */

/*
 * Zone type
 */

/*! initialization zone */
#define CS_VOLUME_ZONE_INITIALIZATION       (1 << 0)

/*! porosity zone */
#define CS_VOLUME_ZONE_POROSITY             (1 << 1)

/*! head loss zone */
#define CS_VOLUME_ZONE_HEAD_LOSS            (1 << 2)

/*! source term (general) */
#define CS_VOLUME_ZONE_SOURCE_TERM          (1 << 3)

/*! source term (mass) */
#define CS_VOLUME_ZONE_MASS_SOURCE_TERM     (1 << 4)

/*! soil used in the groundwater flow module */
#define CS_VOLUME_ZONE_GWF_SOIL             (1 << 5)

/*! @} */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Volume zone description */

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

} cs_volume_zone_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize volume zone structures.
 *
 * This defines a default volume zone. This is the first function of
 * the volume zone handling functions which should be called, and it should
 * only be called after \ref cs_mesh_location_initialize.
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all volume zone structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of volume zones defined.
 */
/*----------------------------------------------------------------------------*/

int
cs_volume_zone_n_zones(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of volume zones which may vary in time.
 *
 * \return  number of zones which may vary in time
 */
/*----------------------------------------------------------------------------*/

int
cs_volume_zone_n_zones_time_varying(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update association of volume zones with a mesh.
 *
 * For time-varying zones, the associated mesh location is updated.
 *
 * \param[in]  mesh_modified  indicate if mesh has been modified
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_build_all(bool  mesh_modified);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new volume zone using a selection criteria string.
 *
 * \param[in]  name       name of location to define
 * \param[in]  criteria   selection criteria for associated elements
 * \param[in]  type_flag  mask of zone category values
 *
 * \return  id of newly defined volume zone
 */
/*----------------------------------------------------------------------------*/

int
cs_volume_zone_define(const char  *name,
                      const char  *criteria,
                      int          type_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new mesh location with an associated selection function.
 *
 * So as to define a subset of mesh entities of a given type, a pointer
 * to a selection function may be given.
 *
 * This requires more programming but allows finer control than selection
 * criteria, as the function has access to the complete mesh structure.
 *
 * \param[in]  name        name of location to define
 * \param[in]  func        pointer to selection function for associated elements
 * \param[in, out]  input  pointer to optional (untyped) value
 *                         or structure.
 * \param[in]  type_flag   mask of zone category values
 *
 * \return  id of newly defined created mesh location
 */
/*----------------------------------------------------------------------------*/

int
cs_volume_zone_define_by_func(const char                 *name,
                              cs_mesh_location_select_t  *func,
                              void                       *input,
                              int                         type_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a volume zone based on its id.
 *
 * This function requires that a volume zone of the given id is defined.
 *
 * \param[in]  id   zone id
 *
 * \return  pointer to the volume zone structure
 */
/*----------------------------------------------------------------------------*/

const cs_volume_zone_t  *
cs_volume_zone_by_id(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a volume zone based on its name if present.
 *
 * This function requires that a volume zone of the given name is defined.
 *
 * \param[in]  name  volume zone name
 *
 * \return  pointer to (read-only) zone structure
 */
/*----------------------------------------------------------------------------*/

const cs_volume_zone_t  *
cs_volume_zone_by_name(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a volume zone based on its name if present.
 *
 * If no volume zone of the given name is defined, NULL is returned.
 *
 * \param[in]  name  volume zone name
 *
 * \return  pointer to (read only) zone structure, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_volume_zone_t  *
cs_volume_zone_by_name_try(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set type flag for a given volume zone.
 *
 * \param[in]  id         volume zone id
 * \param[in]  type_flag  volume zone type flag
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_set_type(int   id,
                        int   type_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set time varying behavior for a given volume zone.
 *
 * \param[in]  id            volume zone id
 * \param[in]  time_varying  true if the zone's definition varies in time
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_set_time_varying(int   id,
                                bool  time_varying);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set overlay behavior for a given volume zone.
 *
 * \param[in]  id             volume zone id
 * \param[in]  allow_overlay  true if the zone may be overlayed by another
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_set_overlay(int   id,
                           bool  allow_overlay)
;
/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to zone id associated with each cell.
 *
 * In case of overlayed zones, the highest zone id associated with
 * a given cell is given.
 */
/*----------------------------------------------------------------------------*/

const int *
cs_volume_zone_cell_zone_id(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to a given volume zone to log file.
 *
 * \param[in]  z   pointer to volume zone structure
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_log_info(const cs_volume_zone_t  *z);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log setup information relative to defined volume zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of volume zones associated with a
 *        given zone flag.
 *
 * \param[in]  type_flag  flag to compare to zone type
 *
 * \return  number of zones matching the given type flag
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_volume_zone_n_type_zones(int  type_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of volume zone cells associated with a
 *        given zone flag.
 *
 * Note that in the case of overlapping zones, a cell may be accounted
 * for multiple times.
 *
 * \param[in]  type_flag  flag to compare to zone type
 *
 * \return  number of cells in zones matching the given type flag
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_volume_zone_n_type_cells(int  type_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select cells associated with volume zones of a given type.
 *
 * Note that in the case of overlapping zones, a cell may be accounted
 * for multiple times.
 *
 * \param[in]   type_flag  flag to compare to zone type
 * \param[out]  cell_ids   ids of selected cells (size: given by
 *                         \ref cs_volume_zone_n_type_cells)
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_select_type_cells(int        type_flag,
                                 cs_lnum_t  cell_ids[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_VOLUME_ZONE_H__ */
