#ifndef __CS_POST_UTIL_H__
#define __CS_POST_UTIL_H__

/*============================================================================
 * Postprocessing utility functions.
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_mesh_location.h"
#include "cs_field_operator.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Type Definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select cells cut by a given segment
 *
 * This selection function may be used as an elements selection function
 * for postprocessing.
 *
 * In this case, the input points to a real array containing the segment's
 * start and end coordinates.
 *
 * Note: the input pointer must point to valid data when this selection
 * function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * The caller is responsible for freeing the returned cell_ids array.
 * When passed to postprocessing mesh or probe set definition functions,
 * this is handled automatically.
 *
 * \deprecated Use cs_mesh_intersect_segment_cell_select (rename) instead.
 *
 * \param[in]   input     pointer to segment start and end:
 *                        [x0, y0, z0, x1, y1, z1]
 * \param[out]  n_cells   number of selected cells
 * \param[out]  cell_ids  array of selected cell ids (0 to n-1 numbering)
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_segment_intersect_select(void        *input,
                                 cs_lnum_t   *n_cells,
                                 cs_lnum_t  **cell_ids);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select cells cut by a line composed of segments
 *
 * This selection function may be used as an elements selection function
 * for postprocessing.
 *
 * In this case, the input points to a real array containing the segment's
 * start and end coordinates.
 *
 * Note: the input pointer must point to valid data when this selection
 * function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * The caller is responsible for freeing the returned cell_ids array.
 * When passed to postprocessing mesh or probe set definition functions,
 * this is handled automatically.
 *
 * \deprecated Use cs_mesh_intersect_polyline_cell_select (rename) instead.
 *
 * \param[in]   input     pointer to segments starts and ends:
 *                        [x0, y0, z0, x1, y1, z1]
 * \param[in]   n_points  number of vertices in the polyline
 * \param[out]  n_cells   number of selected cells
 * \param[out]  cell_ids  array of selected cell ids (0 to n-1 numbering)
 * \param[out]  seg_c_len array of length of the segment in the selected cells
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_polyline_intersect_select(void        *input,
                                  cs_lnum_t   n_points,
                                  cs_lnum_t   *n_cells,
                                  cs_lnum_t  **cell_ids,
                                  cs_real_t  **seg_c_len);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define probes based on the centers of cells intersected by
 *        a given segment.
 *
 * This selection function may be used as a probe set definition function
 * for postprocessing.
 *
 * In this case, the input points to a real array containing the segment's
 * start and end coordinates.
 *
 * Note: the input pointer must point to valid data when this selection
 * function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * The caller is responsible for freeing the returned cell_ids array.
 * When passed to postprocessing mesh or probe set definition functions,
 * this is handled automatically.
 *
 * \deprecated: higher-level cs_probe_set_create_from_segment function with
 *              n_probes argument set at 0 or lower includes equivalent
 *              function.
 *
 * \param[in]   input   pointer to segment start and end:
 *                      [x0, y0, z0, x1, y1, z1]
 * \param[out]  n_elts  number of selected coordinates
 * \param[out]  coords  coordinates of selected elements.
 * \param[out]  s       curvilinear coordinates of selected elements
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_segment_intersect_probes_define(void          *input,
                                        cs_lnum_t     *n_elts,
                                        cs_real_3_t  **coords,
                                        cs_real_t    **s);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a profile based on centers of faces defined by a given
 *        criterion
 *
 * Here, the input points to string describing a selection criterion.
 *
 * \param[in]   input   pointer to selection criterion
 * \param[out]  n_elts  number of selected coordinates
 * \param[out]  coords  coordinates of selected elements.
 * \param[out]  s       curvilinear coordinates of selected elements
 */
/*----------------------------------------------------------------------------*/

void
cs_b_face_criterion_probes_define(void          *input,
                                  cs_lnum_t     *n_elts,
                                  cs_real_3_t  **coords,
                                  cs_real_t    **s);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the head of a turbomachinery (total pressure increase)
 *
 * \param[in]   criteria_in   selection criteria of turbomachinery suction
 * \param[in]   location_in   mesh location of turbomachinery suction
 * \param[in]   criteria_out  selection criteria of turbomachinery discharge
 * \param[in]   location_out  mesh location of turbomachinery discharge
 *
 * \return turbomachinery head
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_turbomachinery_head(const char               *criteria_in,
                            cs_mesh_location_type_t   location_in,
                            const char               *criteria_out,
                            cs_mesh_location_type_t   location_out);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the magnitude of a moment of force torque) given an
 *         axis and the stress on a specific boundary.
 *
 * \param[in]   n_b_faces    number of faces
 * \param[in]   b_face_ids   list of faces (0 to n-1)
 * \param[in]   axis         axis
 *
 * \return couple about the axis
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_moment_of_force(cs_lnum_t        n_b_faces,
                        const cs_lnum_t  b_face_ids[],
                        cs_real_t        axis[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute tangential stress on a specific boundary.
 *
 * \param[in]   n_b_faces    number of faces
 * \param[in]   b_face_ids   list of faces (0 to n-1)
 * \param[out]  stress       tangential stress on the specific
 *                           boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_post_stress_tangential(cs_lnum_t        n_b_faces,
                          const cs_lnum_t  b_face_ids[],
                          cs_real_3_t      stress[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute pressure on a specific boundary region.
 *
 * \param[in]   n_b_faces    number of faces
 * \param[in]   b_face_ids   list of faces (0 to n-1)
 * \param[out]  pres         pressure on a specific boundary region
 */
/*----------------------------------------------------------------------------*/

void
cs_post_b_pressure(cs_lnum_t         n_b_faces,
                   const cs_lnum_t   b_face_ids[],
                   cs_real_t         pres[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute total pressure on a specific boundary region.
 *
 * \param[in]   n_b_faces    number of faces
 * \param[in]   b_face_ids   list of faces (0 to n-1)
 * \param[out]  pres         total pressure on a specific boundary region
 */
/*----------------------------------------------------------------------------*/

void
cs_post_b_total_pressure(cs_lnum_t         n_b_faces,
                         const cs_lnum_t   b_face_ids[],
                         cs_real_t         pres[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute Reynolds stresses in case of Eddy Viscosity Models
 *
 * \param[in]  interpolation_type interpolation type for turbulent kinetic
 *                                energy field
 * \param[in]  n_cells            number of points
 * \param[in]  cell_ids           cell location of points
 *                                 (indexed from 0 to n-1, or NULL if 1-to-1)
 * \param[in]  coords             point coordinates (or NULL for cell centers)
 * \param[out] rst                Reynolds stresses stored as vector
 *                                [r11, r22, r33, r12, r23, r13]
 */
/*----------------------------------------------------------------------------*/

void
cs_post_evm_reynolds_stresses(cs_field_interpolate_t  interpolation_type,
                              cs_lnum_t               n_cells,
                              const cs_lnum_t         cell_ids[],
                              const cs_real_3_t      *coords,
                              cs_real_6_t            *rst);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the invariant of the anisotropy tensor
 *
 * \param[in]  n_cells            number of points
 * \param[in]  cell_ids           cell location of points
 *                                (indexed from 0 to n-1)
 * \param[in]  coords             point coordinates
 * \param[out] inv                Anisotropy tensor invariant
 *                                [xsi, eta]
 */
/*----------------------------------------------------------------------------*/

void
cs_post_anisotropy_invariant(cs_lnum_t               n_cells,
                             const cs_lnum_t         cell_ids[],
                             const cs_real_t         coords[][3],
                             cs_real_t               inv[][2]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute scalar flux on a specific boundary region.
 *
 * The flux is counted negatively through the normal.
 *
 * \param[in]   scalar_name    scalar name
 * \param[in]   n_loc_b_faces  number of selected boundary faces
 * \param[in]   b_face_ids     ids of selected boundary faces
 * \param[out]  b_face_flux    surface flux through selected faces
 */
/*----------------------------------------------------------------------------*/

void
cs_post_boundary_flux(const char       *scalar_name,
                      cs_lnum_t         n_loc_b_faces,
                      const cs_lnum_t   b_face_ids[],
                      cs_real_t         b_face_flux[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute values at a selection of boundary faces of a given field
 *        located on cells.
 *
 * Field BCs are taken into account and boundary cell values are reconstructed
 * using the cell gradient.
 *
 * \param[in]   f              field pointer
 * \param[in]   n_loc_b_faces  number of selected boundary faces
 * \param[in]   b_face_ids     ids of selected boundary faces
 * \param[out]  b_val          values on boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_post_field_cell_to_b_face_values(const cs_field_t  *f,
                                    cs_lnum_t          n_loc_b_faces,
                                    const cs_lnum_t    b_face_ids[],
                                    cs_real_t         *b_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface integral of a scalar array of values located
 *        on boundary faces over a selection of boundary faces
 *
 * \param[in]   scalar_vals    array of scalar values of size n_b_faces
 * \param[in]   n_loc_b_faces  number of selected boundary faces
 * \param[in]   b_face_ids     ids of selected boundary faces
 *
 * \return  Value of computed surface integral
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_scalar_boundary_integral(const cs_real_t *scalar_vals,
                                 const cs_lnum_t  n_loc_b_faces,
                                 const cs_lnum_t  b_face_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface integral over a selection of boundary faces
 *        for a scalar array of values located over the selection.
 *
 * \param[in]   scalar_vals    array of scalar values of size n_b_faces
 * \param[in]   n_loc_b_faces  number of selected boundary faces
 * \param[in]   b_face_ids     ids of selected boundary faces
 *
 * \return  Value of computed surface integral
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_bnd_scalar_boundary_integral(const cs_real_t *scalar_vals,
                                     const cs_lnum_t  n_loc_b_faces,
                                     const cs_lnum_t  b_face_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface mean of a scalar array of values located
 *        on boundary faces over a selection of boundary faces. Weighting
 *        is done using total surface of boundary faces.
 *
 * \param[in]   scalar_vals    array of scalar values of size n_b_faces
 * \param[in]   n_loc_b_faces  number of selected boundary faces
 * \param[in]   b_face_ids     ids of selected boundary faces
 *
 * \return  Value of computed surface mean value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_scalar_boundary_mean(const cs_real_t *scalar_vals,
                             const cs_lnum_t  n_loc_b_faces,
                             const cs_lnum_t  b_face_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface mean over a selection of boundary faces
 *        for a scalar array of values located over the selection.
 *        Weighting is done using total surface of boundary faces.
 *
 * \param[in]   scalar_vals    array of scalar values of size n_b_faces
 * \param[in]   n_loc_b_faces  number of selected boundary faces
 * \param[in]   b_face_ids     ids of selected boundary faces
 *
 * \return  Value of computed surface mean value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_bnd_scalar_boundary_mean(const cs_real_t *scalar_vals,
                                 const cs_lnum_t  n_loc_b_faces,
                                 const cs_lnum_t  b_face_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface integral of a scalar array of values located
 *        on boundary faces over a boundary zone.
 *
 * \param[in]   z              pointer to boundary zone
 * \param[in]   scalar_vals    array of scalar values of size n_b_faces
 *
 * \return  Value of computed surface integral
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_scalar_b_zone_integral(const cs_zone_t *z,
                               const cs_real_t *scalar_vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface integral of a scalar over a boundary zone faces,
 *        for an array of values located on the zone's faces.
 *
 * \param[in]   z              pointer to boundary zone
 * \param[in]   scalar_vals    array of scalar values of size n_b_faces
 *
 * \return  Value of computed surface integral
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_bnd_scalar_b_zone_integral(const cs_zone_t *z,
                                   const cs_real_t *scalar_vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface mean of a scalar array of values located
 *        on boundary faces over a boundary zone. Weighting
 *        is done using total surface of boundary faces.
 *
 * \param[in]   z              pointer to boundary zone
 * \param[in]   scalar_vals    array of scalar values of size n_b_faces
 *
 * \return  Value of computed surface mean value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_scalar_b_zone_mean(const cs_zone_t *z,
                           const cs_real_t *scalar_vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface mean of a scalar over a boundary zone faces,
 *        for an array of values located on the zone's faces.
 *        Weighting is done using total surface of boundary faces.
 *
 * \param[in]   z              pointer to boundary zone
 * \param[in]   scalar_vals    array of scalar values of size n_b_faces
 *
 * \return  Value of computed surface mean value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_bnd_scalar_b_zone_mean(const cs_zone_t *z,
                               const cs_real_t *scalar_vals);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POST_UTIL_H__ */
