#ifndef __CS_MEDCOUPLING_REMAPPER_HXX__
#define __CS_MEDCOUPLING_REMAPPER_HXX__

/*============================================================================
 * Interpolation using MEDCoupling Remapper.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * MED library headers
 *----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Structure definitions
 *============================================================================*/

typedef struct _cs_medcoupling_remapper_t cs_medcoupling_remapper_t;

/*============================================================================
 * Public C function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief get a remapper by its id
 *
 * \param[in] r_id  id of the remapper
 *
 * \return  pointer to cs_medcoupling_remapper_t struct
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_remapper_t *
cs_medcoupling_remapper_by_id(int  r_id);

/* -------------------------------------------------------------------------- */
/*!
 * \brief get a remapper by its name
 *
 * \param[in] name  name of the remapper
 *
 * \return  pointer to cs_medcoupling_remapper_t struct
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_remapper_t *
cs_medcoupling_remapper_by_name_try(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief initialize a remapper based on a set of given arguments
 *
 * \param[in] name             name of the new remapper
 * \param[in] elt_dim          element dimension
 * \param[in] select_criteria  selection criteria for the elements
 * \param[in] medfile_path     path to the med file
 * \param[in] n_fields         number of fields to load
 * \param[in] field_names      names of the fields to load
 * \param[in] iteration        time iteration to load
 * \param[in] order            iteration order to load
 *
 * \return  id of the new remapper
 */
/*----------------------------------------------------------------------------*/

int
cs_medcoupling_remapper_initialize(const char   *name,
                                   int           elt_dim,
                                   const char   *select_criteria,
                                   const char   *medfile_path,
                                   int           n_fields,
                                   const char  **field_names,
                                   int           iteration,
                                   int           order);

/*----------------------------------------------------------------------------*/
/*!
 * \brief set and load a given time iteration from the MED file
 *
 * \param[in] r            pointer to the cs_medcoupling_remapper_t struct
 * \param[in] iteration    time iteration to load
 * \param[in] order        iteration order to load
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_set_iteration(cs_medcoupling_remapper_t  *r,
                                      int                         iteration,
                                      int                         order);

/*----------------------------------------------------------------------------*/
/*!
 * \brief set non-default options for a remapper
 *
 * \param[in] r      pointer to the cs_medcoupling_remapper_t struct
 * \param[in] key    pointer to string representing key
 *                   currently handled: one of {Precision, IntersectionType}
 * \param[in] value  pointer to string representing value:
 *                   - for Precision: floating-point value (default: 1e-12)
 *                   - for IntersectionType: one of {Triangulation, Convex,
 *                     Geometric2D, PointLocator, Barycentric,
 *                     BarycentricGeo2D, MappedBarycentric}
 *                     (see MEDCoupling INTERP_KERNEL documentation)
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_set_options(cs_medcoupling_remapper_t  *r,
                                    const char                  key[],
                                    const char                  value[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief update the interpolation matrix of the remapper
 *
 * \param[in] r            pointer to the cs_medcoupling_remapper_t struct
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_setup(cs_medcoupling_remapper_t  *r);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Interpolate values for a given field
 *
 * \param[in] r            pointer to the cs_medcoupling_remapper_t struct
 * \param[in] field_id     id of the field to interpolate (in the list
 *                         given before)
 * \param[in] default_val  value to apply for elements not intersected by
 *                         source mesh
 *
 * \return  pointer to cs_real_t array containing the new values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_medcoupling_remapper_copy_values(cs_medcoupling_remapper_t  *r,
                                    int                         field_id,
                                    double                      default_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief translate the mesh using a given vector
 *
 * \param[in] r            pointer to the cs_medcoupling_remapper_t struct
 * \param[in] translation  translation vector
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_translate(cs_medcoupling_remapper_t  *r,
                                  cs_real_t                   translation[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Rotate the mesh using a center point, axis and angle
 *
 * \param[in] r          pointer to the cs_medcoupling_remapper_t struct
 * \param[in] invariant  coordinates of the invariant point
 * \param[in] axis       rotation axis vector
 * \param[in] angle      rotation angle in radians
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_rotate(cs_medcoupling_remapper_t  *r,
                               cs_real_t                   invariant[3],
                               cs_real_t                   axis[3],
                               cs_real_t                   angle);

/*----------------------------------------------------------------------------*/
/*! \brief Retrieve the two closest time steps indexes.
 *
 * The returned value is int[2].
 * If the requested time value if outside the time bounds stored in the file,
 * the both values are identical (first or last value), and a warning is printed
 * in the listing file.
 *
 * \param[in]      r    pointer to remapper object
 * \param[in]      t    requested time value
 * \param[in,out]  id1  first returned index
 * \param[in,out]  id2  second returned index
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_find_time_index(cs_medcoupling_remapper_t *r,
                                        cs_real_t                  t,
                                        int                       *id1,
                                        int                       *id2);

/*----------------------------------------------------------------------------*/
/*! \brief Retrieve the two closest time steps indexes.
 *
 * The returned value is int[2].
 * If the requested time value if outside the time bounds stored in the file,
 * the both values are identical (first or last value), and a warning is printed
 * in the listing file.
 *
 * \param[in]      r    pointer to remapper object
 * \param[in]      id   requested index
 * \param[in,out]  t    corresponding time value
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_get_time_from_index(cs_medcoupling_remapper_t *r,
                                            int                        id,
                                            cs_real_t                 *t);

/*----------------------------------------------------------------------------*/
/*! \brief Retrieve the two closest time steps indexes.
 *
 * The returned value is int[2].
 * If the requested time value if outside the time bounds stored in the file,
 * the both values are identical (first or last value), and a warning is output
 * in the lod file.
 *
 * \param[in]      r      pointer to remapper object
 * \param[in]      id     requested time index
 * \param[in,out]  it     index iteration
 * \param[in,out]  order  index iteration order
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_get_iter_order_from_index(cs_medcoupling_remapper_t *r,
                                                  int                        id,
                                                  int                       *it,
                                                  int                       *order);

/*----------------------------------------------------------------------------*/
/*! \brief Load the time value corresponding to id.
 *
 * \param[in]      r      pointer to remapper object
 * \param[in]      id     requested time index
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_update_time_value(cs_medcoupling_remapper_t *r,
                                          int                        id);
/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all remappers
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_destroy_all(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MEDCOUPLING_REMAPPER_HXX__ */
