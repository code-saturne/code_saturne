#ifndef __CS_BOUNDARY_ZONE_H__
#define __CS_BOUNDARY_ZONE_H__

/*============================================================================
 * Boundary zones handling.
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
 * @defgroup volume_zone_flags Flags specifying general boundary zone type
 *
 * @{
 */

/*
 * Zone type
 */

/*! initialization zone */
#define CS_BOUNDARY_ZONE_WALL       (1 << 0)

/*! @} */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Boundary zone description */

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


} cs_boundary_zone_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize boundary zone structures.
 *
 * This defines a default boundary zone. This is the first function of
 * the boundary zone handling functions which should be called, and it should
 * only be called after \ref cs_mesh_location_initialize.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_zone_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all boundary zone structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_zone_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of boundary zones defined.
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_zone_n_zones(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of boundary zones which may vary in time.
 *
 * \return  number of zones which may vary in time
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_zone_n_zones_time_varying(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update association of boundary zones with a mesh.
 *
 * For time-varying zones, the associated mesh location is updated.
 *
 * \param[in]  mesh_modified  indicate if mesh has been modified
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_zone_build_all(bool  mesh_modified);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new boundary zone using a selection criteria string.
 *
 * \param[in]  name       name of location to define
 * \param[in]  criteria   selection criteria for associated elements
 * \param[in]  type_flag  mask of zone category values
 *
 * \return  id of newly defined boundary zone
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_zone_define(const char  *name,
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
 * \param[in]  name  name of location to define
 * \param[in]  func  pointer to selection function for associated elements
 * \param[in, out]  input  pointer to optional (untyped) value
 *                         or structure.
 * \param[in]  type_flag  mask of zone category values
 *
 * \return  id of newly defined created mesh location
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_zone_define_by_func(const char                 *name,
                                cs_mesh_location_select_t  *func,
                                void                       *input,
                                int                         type_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a boundary zone based on its id.
 *
 * This function requires that a boundary zone of the given id is defined.
 *
 * \param[in]  id   zone id
 *
 * \return  pointer to the boundary zone structure
 */
/*----------------------------------------------------------------------------*/

const cs_boundary_zone_t  *
cs_boundary_zone_by_id(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a boundary zone based on its name if present.
 *
 * This function requires that a boundary zone of the given name is defined.
 *
 * \param[in]  name  boundary zone name
 *
 * \return  pointer to (read-only) zone structure
 */
/*----------------------------------------------------------------------------*/

const cs_boundary_zone_t  *
cs_boundary_zone_by_name(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a boundary zone based on its name if present.
 *
 * If no boundary zone of the given name is defined, NULL is returned.
 *
 * \param[in]  name  boundary zone name
 *
 * \return  pointer to (read only) zone structure, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_boundary_zone_t  *
cs_boundary_zone_by_name_try(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set type flag for a given boundary zone.
 *
 * \param[in]  id         boundary zone id
 * \param[in]  type_flag  volume zone type flag
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_zone_set_type(int   id,
                          int   type_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set time varying behavior for a given boundary zone.
 *
 * \param[in]  id            boundary zone id
 * \param[in]  time_varying  true if the zone's definition varies in time
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_zone_set_time_varying(int   id,
                                  bool  time_varying);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set overlay behavior for a given boundary zone.
 *
 * \param[in]  id             boundary zone id
 * \param[in]  allow_overlay  true if the zone may be overlayed by another
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_zone_set_overlay(int   id,
                             bool  allow_overlay)
;
/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to zone id associated with each boundary face.
 *
 * In case of overlayed zones, the highest zone id associated with
 * a given face is given.
 */
/*----------------------------------------------------------------------------*/

const int *
cs_boundary_zone_face_zone_id(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to a given boundary zone to log file.
 *
 * \param[in]  z   pointer to boundary zone structure
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_zone_log_info(const cs_boundary_zone_t  *z);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log setup information relative to defined boundary zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_zone_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of boundary zones associated with a
 *        given zone flag.
 *
 * \param[in]  type_flag  flag to compare to zone type
 *
 * \return  number of zones matching the given type flag
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_boundary_zone_n_type_zones(int  type_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to optional boundary face class ids.
 *
 * For each boundary face, a specific output (logging and postprocessing)
 * class id may be assigned. This allows realizing logging, postprocessing,
 * or otherwise extracting data based on this class.
 *
 * Using this function at a given point indicates that user-defined class
 * ids will be used. The face class ids are initially equal to the
 * face zone ids, but may be modified by the user.
 *
 * In the presence of a time-varying mesh or boundary zones, the face
 * class ids will be reset to the zone ids, so it may be necessary to
 * update the user definitions.
 *
 * The class id values are arbitrarily chosen by the user, but must be
 * positive integers; numbers do not need to be contiguous, but very high
 * numbers may also lead to higher memory consumption.
 *
 * \return  pointer to array of boundary face output zone ids;
 */
/*----------------------------------------------------------------------------*/

int *
cs_boundary_zone_face_class_id(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get read pointer to optional boundary face class or zone ids.
 *
 * If no face classes have been defined by \ref cs_boundary_zone_face_class_id
 * the boundary face zone id is returned instead.
 *
 * \return  pointer to array of boundary face output zone ids;
 */
/*----------------------------------------------------------------------------*/

const int *
cs_boundary_zone_face_class_or_zone_id(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update boundary face output class ids if present.
 *
 * Face class ids lower than 0 are replaced by the matching face zone id.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_zone_update_face_class_id(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the maximum defined face class or zone id.
 *
 * \return  maximum face class or zone id;
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_zone_max_class_or_zone_id(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_ZONE_H__ */
