#ifndef __CS_PARAMEDMEM_REMAPPER_HXX__
#define __CS_PARAMEDMEM_REMAPPER_HXX__

/*============================================================================
 * Parallel interpolation using ParaMEDMEM OverlapDEC
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

typedef struct _cs_paramedmem_remapper_t cs_paramedmem_remapper_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Creates a new cs_paramedmem_remapper_t instance
 *
 * \param[in] name          name of the remapper
 * \param[in] sel_criteria  cells selection criteria
 * \param[in] file_name     med file name
 * \param[in] mesh_name     name of the mesh in the med file
 * \param[in] center        center of bounding sphere
 * \param[in] radius        radius of bounding sphere
 *
 * \return  cs_paramedmem_remapper_t struct
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_remapper_t *
cs_paramedmem_remapper_create(char        *name,
                              const char  *sel_criteria,
                              char        *file_name,
                              char        *mesh_name,
                              cs_real_t    center[3],
                              cs_real_t    radius);

/*----------------------------------------------------------------------------*/
/*!
 * \brief get a remapper by its name
 *
 * \param[in] name  name of the remapper
 *
 * \return  pointer to cs_paramedmem_remapper_t struct
 */
/*----------------------------------------------------------------------------*/

cs_paramedmem_remapper_t *
cs_paramedmem_remapper_by_name_try(const char *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Interpolate a given field on the local mesh for a given time
 *
 * \param[in] r             pointer to cs_paramedmem_remapper_t struct
 * \param[in] field_name    name of the field to remap from the file
 * \param[in] default_val  default value for unmapped elements
 * \param[in] time_choice   Choice of the time interpolation.
 *                          0: Value of field interpolated at t=tval from the
 *                          med file.
 *                          1: Returns field values for the first time step in
 *                          the file. tval is then ignored.
 *                          2: Returns field values for the last time step in
 *                          the file. tval is then ignored.
 * \param[in] tval          requested time instant. If time choice is 0 and
 *                          tval outside of the file time bounds, return value
 *                          will be at the the first time step (if tval < tmin)
 *                          or last time step (if tval > tmax)
 *
 * \return  cs_real_t pointer containing the new values on target mesh
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_paramedmem_remap_field(cs_paramedmem_remapper_t *r,
                          char                     *field_name,
                          cs_real_t                 default_val,
                          int                       time_choice,
                          double                    tval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Remaps a field from the med file to the local mesh for a given time
 *
 * \param[in] r           pointer to cs_paramedmem_remapper_t struct
 * \param[in] field_name  name of the field to remap from the file
 * \param[in] default_val default value for unmapped elements
 * \param[in] dt          time value to use from the file
 * \param[in] it          time iteration to use from the file
 *
 * \return  cs_real_t pointer containing the new values on target mesh
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_paramedmem_remap_field_one_time(cs_paramedmem_remapper_t *r,
                                   char                     *field_name,
                                   cs_real_t                 default_val,
                                   int                       dt,
                                   int                       it);

/*----------------------------------------------------------------------------*/
/*!
 * \brief translate the mesh using a given vector
 *
 * \param[in] r            pointer to the cs_paramedmem_remapper_t struct
 * \param[in] translation  translation vector
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_remapper_translate(cs_paramedmem_remapper_t  *r,
                                 cs_real_t                  translation[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Rotate the mesh using a center point, axis and angle
 *
 * \param[in] r          pointer to the cs_paramedmem_remapper_t struct
 * \param[in] invariant  coordinates of the invariant point
 * \param[in] axis       rotation axis vector
 * \param[in] angle      rotation angle in radians
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_remapper_rotate(cs_paramedmem_remapper_t  *r,
                              cs_real_t                  invariant[3],
                              cs_real_t                  axis[3],
                              cs_real_t                  angle);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all remappers
 */
/*----------------------------------------------------------------------------*/

void
cs_paramedmem_remapper_destroy_all(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAMEDMEM_REMAPPER_HXX__ */
