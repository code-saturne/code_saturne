#ifndef __CS_MEDCOUPLING_INTERSECTOR_HXX__
#define __CS_MEDCOUPLING_INTERSECTOR_HXX__

/*============================================================================
 * Interpolation using MEDCoupling Intersector.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

BEGIN_C_DECLS

/*============================================================================
 * Structure definitions
 *============================================================================*/

typedef struct _cs_medcoupling_intersector_t cs_medcoupling_intersector_t;

/*============================================================================
 * Structure definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a MEDCoupling intersector.
 *
 * \param[in] name             name of the intersector
 * \param[in] medfile_path     path to the MED file
 * \param[in] interp_method    interpolation method (P0P0, P1P0, ..)
 * \param[in] select_criteria  selection criteria
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_add(const char  *name,
                               const char  *medfile_path,
                               const char  *interp_method,
                               const char  *select_criteria);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a given MEDCoupling intersector.
 *
 * \param[in]  mi  pointer to the cs_medcoupling_intersector_t struct
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_destroy(cs_medcoupling_intersector_t  *mi);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all allocated intersectors.
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a MEDCoupling intersector using its id.
 *
 * \param[in] id  id of the intersector
 *
 * \return pointer to the cs_medcoupling_intersector_t or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_intersector_t *
cs_medcoupling_intersector_by_id(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a MEDCoupling intersector by name.
 *
 * \param[in] name  name of the intersector
 *
 * \return pointer to the cs_medcoupling_intersector_t or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_intersector_t *
cs_medcoupling_intersector_by_name(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the intersection volumes between the source mesh and
 * code mesh
 *
 * \param[in] mi            pointer to the cs_medcoupling_intersector_t struct
 *
 * \return a pointer to the array containing the intersected volume of each cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_medcoupling_intersect_volumes(cs_medcoupling_intersector_t  *mi);

/*----------------------------------------------------------------------------*/
/*!
 * \brief translate the mesh using a given vector
 *
 * \param[in] mi           pointer to the cs_medcoupling_intersector_t struct
 * \param[in] translation  translation vector
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_translate(cs_medcoupling_intersector_t  *mi,
                                     cs_real_t  translation[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief rotate the mesh
 *
 * \param[in] mi         pointer to the cs_medcoupling_intersector_t struct
 * \param[in] invariant  Invariant point
 * \param[in] axis       Rotation axis
 * \param[in] angle      angle (in radians)
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_rotate(cs_medcoupling_intersector_t  *mi,
                                  cs_real_t                      invariant[3],
                                  cs_real_t                      axis[3],
                                  cs_real_t                      angle);

/*----------------------------------------------------------------------------*/
/*!
 * \brief translate the mesh using a given vector
 *
 * \param[in] mi      pointer to the cs_medcoupling_intersector_t struct
 * \param[in] prefix  subdir prefix
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_dump_mesh(cs_medcoupling_intersector_t  *mi,
                                     const char                    *prefix);

END_C_DECLS

#endif /* __CS_MEDCOUPLING_INTERSECTOR_HXX__ */
