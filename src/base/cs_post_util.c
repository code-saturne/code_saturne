/*============================================================================
 * Postprocessing utility functions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_array_reduce.h"
#include "cs_balance_by_zone.h"
#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_geom.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_intersect.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_physical_model.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_selector.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"
#include "cs_velocity_pressure.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_post_util.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
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
                                 cs_lnum_t  **cell_ids)
{
  cs_mesh_intersect_segment_cell_select(input,
                                        n_cells,
                                        cell_ids);
}

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
                                  cs_real_t  **seg_c_len)
{
  cs_mesh_intersect_polyline_cell_select(input,
                                         n_points,
                                         n_cells,
                                         cell_ids,
                                         seg_c_len);
}

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
                                        cs_real_t    **s)
{
  cs_real_t *sx = (cs_real_t *)input;

  const cs_real_t dx1[3] = {sx[3]-sx[0], sx[4]-sx[1], sx[5]-sx[2]};
  const cs_real_t s_norm2 = cs_math_3_square_norm(dx1);

  const cs_real_3_t  *cell_cen
    = (const cs_real_3_t *)(cs_glob_mesh_quantities->cell_cen);

  cs_lnum_t n_cells = 0;
  cs_lnum_t *cell_ids = NULL;
  cs_real_t *seg_c_len = NULL;

  /* This version is better than cs_cell_segment_intersect_select
     because it gives the cell
     if the segment is included in this cell */
  cs_cell_polyline_intersect_select(input, 2, &n_cells, &cell_ids, &seg_c_len);

  cs_real_3_t *_coords;
  cs_real_t *_s;
  BFT_MALLOC(_coords, n_cells, cs_real_3_t);
  BFT_MALLOC(_s, n_cells, cs_real_t);

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    cs_real_t dx[3], coo[3];
    for (cs_lnum_t j = 0; j < 3; j++) {
      coo[j] = cell_cen[cell_ids[i]][j];
      dx[j] = coo[j] - sx[j];
      _coords[i][j] = coo[j];
    }
    _s[i] = cs_math_3_dot_product(dx, dx1) / s_norm2;
  }

  BFT_FREE(cell_ids);
  BFT_FREE(seg_c_len);

  /* Set return values */

  *n_elts = n_cells;
  *coords = _coords;
  *s = _s;
}

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
                                  cs_real_t    **s)
{
  const char *criterion = (const char *)input;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_lnum_t   n_faces;
  cs_lnum_t  *face_ids;

  BFT_MALLOC(face_ids, m->n_b_faces, cs_lnum_t);
  cs_selector_get_b_face_list(criterion, &n_faces, face_ids);

  cs_real_3_t *_coords;
  cs_real_t *_s;
  BFT_MALLOC(_coords, n_faces, cs_real_3_t);
  BFT_MALLOC(_s, n_faces, cs_real_t);

  for (cs_lnum_t i = 0; i < n_faces; i++) {
    for (cs_lnum_t j = 0; j < 3; j++)
      _coords[i][j] = mq->b_face_cog[face_ids[i]*3 + j];
    _s[i] = _coords[i][0];
  }

  BFT_FREE(face_ids);

  /* Set return values */

  *n_elts = n_faces;
  *coords = _coords;
  *s = _s;
}

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
                            cs_mesh_location_type_t   location_out)
{
  cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;

  cs_real_t *total_pressure = cs_field_by_name("total_pressure")->val;
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_t *density = CS_F_(rho)->val;

  cs_real_t pabs_in = 0.;
  cs_real_t sum_in = 0.;
  cs_real_t pabs_out = 0.;
  cs_real_t sum_out = 0.;

  for (int _n = 0; _n < 2; _n++) {

    cs_lnum_t n_elts = 0;
    cs_lnum_t *elt_list = NULL;
    cs_real_t pabs = 0.;
    cs_real_t sum = 0.;

    cs_mesh_location_type_t location;
    const char *criteria = NULL;

    if (_n == 0) {
      location = location_in;
      criteria = criteria_in;
    } else {
      location = location_out;
      criteria = criteria_out;
    }

    switch(location) {
    case CS_MESH_LOCATION_CELLS:

      BFT_MALLOC(elt_list, mesh->n_cells, cs_lnum_t);
      cs_selector_get_cell_list(criteria, &n_elts, elt_list);

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t cell_id = elt_list[i];
        cs_real_t weight = mesh_quantities->cell_vol[cell_id];
        pabs += weight*(total_pressure[cell_id] + 0.5*density[cell_id]*
                        cs_math_3_square_norm(vel[cell_id]));
        sum += weight;
      }
      BFT_FREE(elt_list);
      break;

    case CS_MESH_LOCATION_BOUNDARY_FACES:

      BFT_MALLOC(elt_list, mesh->n_b_faces, cs_lnum_t);
      cs_selector_get_b_face_list(criteria, &n_elts, elt_list);

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t face_id = elt_list[i];
        cs_lnum_t cell_id = mesh->b_face_cells[face_id];
        cs_real_t surf = mesh_quantities->b_face_surf[face_id];
        pabs += surf*(total_pressure[cell_id] + 0.5*density[cell_id]
                      *cs_math_3_square_norm(vel[cell_id]));
        sum += surf;
      }
      BFT_FREE(elt_list);
      break;

    case CS_MESH_LOCATION_INTERIOR_FACES:

      BFT_MALLOC(elt_list, mesh->n_i_faces, cs_lnum_t);
      cs_selector_get_i_face_list(criteria, &n_elts, elt_list);

      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t face_id = elt_list[i];
        cs_lnum_t c_i = mesh->i_face_cells[face_id][0];
        cs_lnum_t c_j = mesh->i_face_cells[face_id][1];
        cs_real_t w = mesh_quantities->i_face_surf[face_id];

        cs_real_t pt = w*total_pressure[c_i] + (1.-w)*total_pressure[c_j];
        cs_real_t r = w*density[c_i] + (1.-w)*density[c_j];
        cs_real_3_t v = {w*vel[c_i][0] + (1.-w)*vel[c_j][0],
                         w*vel[c_i][1] + (1.-w)*vel[c_j][1],
                         w*vel[c_i][2] + (1.-w)*vel[c_j][2]};
        pabs += w*(pt + 0.5*r*cs_math_3_square_norm(v));
        sum += w;
      }
      BFT_FREE(elt_list);
      break;

    default:
      pabs = 0.;
      sum = 1.;
      bft_printf
        (_("Warning: while post-processing the turbomachinery head.\n"
           "         Mesh location %d is not supported, so the computed head\n"
           "         is erroneous.\n"
           "         The %s parameters should be checked.\n"),
           location, __func__);
      break;
    }

    if (_n == 0) {
      pabs_in = pabs;
      sum_in = sum;
    } else {
      pabs_out = pabs;
      sum_out = sum;
    }

  }

  double _s[4] = {pabs_in, pabs_out, sum_in, sum_out};
  cs_parall_sum(4, CS_DOUBLE, _s);

  pabs_in  = _s[0] / _s[2];
  pabs_out = _s[1] / _s[3];

  return pabs_out - pabs_in;
}

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
                        cs_real_t        axis[3])
{
  const cs_real_3_t *b_face_cog
    = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;
  const cs_real_3_t *b_forces
    = (const cs_real_3_t *)cs_field_by_name("boundary_forces")->val;

  cs_real_3_t moment = {0., 0., 0.};

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    cs_real_3_t m;
    cs_lnum_t face_id = b_face_ids[i];
    cs_math_3_cross_product(b_face_cog[face_id], b_forces[face_id], m);

    /* b_forces is the stress on the solid boundary,
       thus it comes with a '-' sign here */
    for (int j = 0; j < 3; j++)
      moment[j] -= m[j];
  }
  cs_parall_sum(3, CS_DOUBLE, moment);

  return cs_math_3_dot_product(moment, axis);
}

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
                          cs_real_3_t      stress[])
{
  const cs_real_3_t *b_face_normal =
    (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_normal;
  const cs_real_t *b_face_surf = cs_glob_mesh_quantities->b_face_surf;
  const cs_real_3_t *forbr =
    (const cs_real_3_t *)cs_field_by_name("boundary_forces")->val;
  cs_lnum_t ifac;
  cs_real_t srfbn, srfnor[3], fornor;

  for (cs_lnum_t iloc = 0 ; iloc < n_b_faces; iloc++) {
    ifac = b_face_ids[iloc];
    srfbn = b_face_surf[ifac];
    srfnor[0] = b_face_normal[ifac][0] / srfbn;
    srfnor[1] = b_face_normal[ifac][1] / srfbn;
    srfnor[2] = b_face_normal[ifac][2] / srfbn;
    fornor =   forbr[ifac][0]*srfnor[0]
             + forbr[ifac][1]*srfnor[1]
             + forbr[ifac][2]*srfnor[2];
    stress[iloc][0] = (forbr[ifac][0] - fornor*srfnor[0]) / srfbn;
    stress[iloc][1] = (forbr[ifac][1] - fornor*srfnor[1]) / srfbn;
    stress[iloc][2] = (forbr[ifac][2] - fornor*srfnor[2]) / srfbn;
  }
}

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
                   cs_real_t         pres[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *diipb = (const cs_real_3_t *)mq->diipb;
  cs_real_3_t *gradp;

  BFT_MALLOC(gradp, m->n_cells_with_ghosts, cs_real_3_t);

  int hyd_p_flag = cs_glob_velocity_pressure_param->iphydr;
  cs_real_3_t *f_ext = (hyd_p_flag == 1) ?
    (cs_real_3_t *)cs_field_by_name_try("volume_forces")->val:NULL;

  cs_field_gradient_potential(CS_F_(p),
                              false, /* use_previous_t */
                              1,     /* inc */
                              hyd_p_flag,
                              f_ext,
                              gradp);

  for (cs_lnum_t iloc = 0 ; iloc < n_b_faces; iloc++) {
    cs_lnum_t face_id = b_face_ids[iloc];
    cs_lnum_t cell_id = m->b_face_cells[face_id];

    cs_real_t pip =   CS_F_(p)->val[cell_id]
                    + cs_math_3_dot_product(gradp[cell_id],
                                            diipb[face_id]);
    pres[iloc] =   CS_F_(p)->bc_coeffs->a[face_id]
                 + CS_F_(p)->bc_coeffs->b[face_id]*pip;


  }
  BFT_FREE(gradp);

  const cs_turb_model_t  *turb_model = cs_glob_turb_model;

  if (   turb_model->itytur == 2
      || turb_model->itytur == 5
      || turb_model->itytur == 6) {
    cs_real_3_t *gradk;
    BFT_MALLOC(gradk, m->n_cells_with_ghosts, cs_real_3_t);

    cs_field_gradient_scalar(CS_F_(k),
                             false,  /* use_previous_t */
                             1,      /* inc */
                             gradk);

    for (cs_lnum_t iloc = 0 ; iloc < n_b_faces; iloc++) {
      cs_lnum_t face_id = b_face_ids[iloc];
      cs_lnum_t cell_id = m->b_face_cells[face_id];

      cs_real_t kip =   CS_F_(k)->val[cell_id]
        + cs_math_3_dot_product(gradk[cell_id],
                                diipb[face_id]);
      pres[iloc] -= 2./3.*CS_F_(rho_b)->val[face_id]
                         *(  CS_F_(k)->bc_coeffs->a[face_id]
                           + CS_F_(k)->bc_coeffs->b[face_id]*kip);
    }
    BFT_FREE(gradk);
  }
}

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
                         cs_real_t         pres[])
{
  /* Compute boundary values for resolved pressure */
  cs_post_b_pressure(n_b_faces, b_face_ids, pres);

  /* Add components related to referenece pressure and gravity */
  const cs_real_t Pref = cs_glob_fluid_properties->p0
                       - cs_glob_fluid_properties->pred0;

  const cs_real_t rho0 = cs_glob_fluid_properties->ro0;

  const cs_real_t *g = cs_glob_physical_constants->gravity;

  const cs_real_t *xyzp0 = cs_glob_fluid_properties->xyzp0;

  const cs_real_3_t *xyz = (cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;

  for (cs_lnum_t iloc = 0; iloc < n_b_faces; iloc++) {
    cs_lnum_t face_id = b_face_ids[iloc];
    pres[face_id] += Pref + rho0 * cs_math_3_distance_dot_product(xyzp0,
                                                                  xyz[face_id],
                                                                  g);
  }
}

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
                              cs_real_6_t            *rst)
{
  const cs_turb_model_t  *turb_model = cs_glob_turb_model;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  if (   turb_model->itytur != 2
      && turb_model->itytur != 6
      && turb_model->itytur != 5)
    bft_error(__FILE__, __LINE__, 0,
              _("This post-processing utility function is only available for "
                "Eddy Viscosity Models."));

  /* velocity gradient */

  cs_real_33_t *gradv;
  BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

  cs_field_gradient_vector(CS_F_(vel),
                           false,  /* use_previous_t */
                           1,      /* inc */
                           gradv);

  const cs_real_t *xk = CS_F_(k)->val;
  cs_real_t *_xk = NULL;

  if (cell_ids != NULL) {
    BFT_MALLOC(_xk, n_cells, cs_real_t);

    if (coords != NULL)
      cs_field_interpolate(CS_F_(k),
                           interpolation_type,
                           n_cells,
                           cell_ids,
                           coords,
                           _xk);
    else {
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        _xk[i] = xk[cell_ids[i]];
      }
    }

    xk = _xk;
  }

  /* Compute Reynolds stresses */

  const cs_real_t d2s3 = 2./3.;
  const cs_real_t *cpro_mu_t = CS_F_(mu_t)->val;
  const cs_real_t *cpro_rho = CS_F_(rho)->val;
  for (cs_lnum_t iloc = 0; iloc < n_cells; iloc++) {
    cs_lnum_t iel = cell_ids[iloc];

    cs_real_t divu = gradv[iel][0][0] + gradv[iel][1][1] + gradv[iel][2][2];
    cs_real_t nut = cpro_mu_t[iel] / cpro_rho[iel];

    cs_real_t xdiag = d2s3*(xk[iloc]+ nut*divu);
    rst[iloc][0] =  xdiag - 2.*nut*gradv[iel][0][0];
    rst[iloc][1] =  xdiag - 2.*nut*gradv[iel][1][1];
    rst[iloc][2] =  xdiag - 2.*nut*gradv[iel][2][2];
    rst[iloc][3] = -nut*(gradv[iel][1][0]+gradv[iel][0][1]);
    rst[iloc][4] = -nut*(gradv[iel][2][1]+gradv[iel][1][2]);
    rst[iloc][5] = -nut*(gradv[iel][2][0]+gradv[iel][0][2]);
  }

  BFT_FREE(gradv);
  BFT_FREE(_xk);
}

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
                             cs_real_2_t             inv[])
{
  const cs_turb_model_t  *turb_model = cs_glob_turb_model;

  if (   turb_model->itytur != 2
      && turb_model->itytur != 3
      && turb_model->itytur != 6
      && turb_model->itytur != 5)
    bft_error(__FILE__, __LINE__, 0,
              _("This post-processing utility function is only available for "
                "RANS Models."));

  cs_real_6_t *rij = NULL;
  BFT_MALLOC(rij, n_cells, cs_real_6_t);
  cs_field_interpolate_t interpolation_type = CS_FIELD_INTERPOLATE_MEAN;

  /* Compute the Reynolds Stresses if we are using EVM */
  if (   turb_model->order == CS_TURB_FIRST_ORDER
      && turb_model->type  == CS_TURB_RANS) {
    cs_post_evm_reynolds_stresses(interpolation_type,
                                  n_cells,
                                  cell_ids,
                                  coords, /* coords */
                                  rij);
  } else {
    cs_real_6_t *cvar_rij = (cs_real_6_t *)CS_F_(rij)->val;
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      cs_lnum_t c_id = cell_ids[i];
      for (cs_lnum_t j = 0; j < 6; j++)
        rij[i][j] = cvar_rij[c_id][j];
    }
  }

  /* Compute Invariants */

  const cs_real_t d1s3 = 1./3.;
  for (cs_lnum_t iloc = 0; iloc < n_cells; iloc++) {
    cs_lnum_t iel = cell_ids[iloc];

    cs_real_t xk = 0.5*(rij[iel][0]+rij[iel][1]+rij[iel][2]);
    cs_real_t bij[3][3];
    cs_real_t xeta, xksi;

    bij[0][0] = rij[iel][0]/(2.0*xk) - d1s3;
    bij[1][1] = rij[iel][1]/(2.0*xk) - d1s3;
    bij[2][2] = rij[iel][2]/(2.0*xk) - d1s3;
    bij[0][1] = rij[iel][3]/(2.0*xk);
    bij[1][2] = rij[iel][4]/(2.0*xk);
    bij[0][2] = rij[iel][5]/(2.0*xk);
    bij[1][0] = bij[0][1];
    bij[2][1] = bij[1][2];
    bij[2][0] = bij[0][2];

    xeta = 0.;
    xksi = 0.;
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        xeta += cs_math_pow2(bij[i][j]);
        for (cs_lnum_t k = 0; k < 3; k++)
          xksi += bij[i][j]*bij[j][k]*bij[k][i];
      }
    }

    inv[iloc][0] = sqrt(xeta/6.0);
    inv[iloc][1] = cbrt(xksi/6.0);
  }

  BFT_FREE(rij);
}

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
                      cs_real_t         b_face_flux[])
{
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;

  cs_real_t normal[] = {0, 0, 0};

  cs_flux_through_surface(scalar_name,
                          normal,
                          n_loc_b_faces,
                          0,
                          b_face_ids,
                          NULL,
                          NULL,
                          b_face_flux,
                          NULL);

  if (b_face_ids != NULL) {
    for (cs_lnum_t i = 0; i < n_loc_b_faces; i++) {
      cs_lnum_t f_id = b_face_ids[i];
      b_face_flux[i] /= b_face_surf[f_id];
    }
  }
  else {
    for (cs_lnum_t f_id = 0; f_id < n_loc_b_faces; f_id++) {
      b_face_flux[f_id] /= b_face_surf[f_id];
    }
  }
}

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
                                    cs_real_t         *b_val)
{
  if (f->location_id != CS_MESH_LOCATION_CELLS)
    bft_error(__FILE__, __LINE__, 0,
              _("Postprocessing face boundary values for field %s :\n"
                " not implemented for fields on location %s."),
              f->name, cs_mesh_location_type_name[f->location_id]);

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  const cs_lnum_t dim = f->dim;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  int coupled = 0;
  if (f->type & CS_FIELD_VARIABLE && f->dim > 1) {
    int coupled_key_id = cs_field_key_id_try("coupled");
    if (coupled_key_id > -1)
      coupled = cs_field_get_key_int(f, coupled_key_id);
  }

  if (dim == 1) {

    cs_real_3_t *grad;
    BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

    cs_field_gradient_scalar(f,
                             true, /* use_previous_t */
                             1,    /* not an increment */
                             grad);

    for (cs_lnum_t ii = 0; ii < n_loc_b_faces; ii++) {
      cs_lnum_t face_id = b_face_ids[ii];
      cs_lnum_t cell_id = b_face_cells[face_id];

      cs_real_t c =   f->val[cell_id]
                    + cs_math_3_dot_product(diipb[face_id], grad[cell_id]);

      b_val[ii] =   f->bc_coeffs->a[face_id]
                  + f->bc_coeffs->b[face_id] * c;
    }

    BFT_FREE(grad);

  }

  else if (dim == 3 || dim == 6) {

    cs_real_t *grad;
    BFT_MALLOC(grad, 3*dim*n_cells_ext, cs_real_t);

    if (dim == 3)
      cs_field_gradient_vector(f,
                               true, /* use_previous_t */
                               1,    /* not an increment */
                               (cs_real_33_t *)grad);
    else if (dim == 6)
      cs_field_gradient_tensor(f,
                               true, /* use_previous_t */
                               1,    /* not an increment */
                               (cs_real_63_t *)grad);

    for (cs_lnum_t ii = 0; ii < n_loc_b_faces; ii++) {
      cs_lnum_t face_id = b_face_ids[ii];
      cs_lnum_t cell_id = b_face_cells[face_id];

      cs_real_t val_ip[3];
      for (cs_lnum_t j = 0; j < dim; j++) {
        cs_lnum_t k = (cell_id*dim + j)*3;
        val_ip[j] =   f->val[cell_id*dim + j]
                    + diipb[face_id][0] * grad[k]
                    + diipb[face_id][1] * grad[k+1]
                    + diipb[face_id][2] * grad[k+2];
      }

      for (int j = 0; j < dim; j++) {
        b_val[ii*dim + j] = f->bc_coeffs->a[dim*face_id + j];

        if (coupled)
          for (int k = 0; k < dim; k++)
            b_val[ii*dim + j]
              += f->bc_coeffs->b[dim*dim*face_id+j*dim+k]*val_ip[k];
        else
          b_val[ii*dim + j] += f->bc_coeffs->b[dim*face_id + j]*val_ip[j];
      }
    }

    BFT_FREE(grad);
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("Postprocessing face boundary values for field %s"
                " of dimension %d:\n not implemented."),
              f->name, f->dim);
}

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
                                 const cs_lnum_t  b_face_ids[])
{
  cs_real_t retval = 0.;

  double _l_sum = 0.;

  if (n_loc_b_faces > 0) {
    const cs_real_t *restrict b_face_surf
      = cs_glob_mesh_quantities->b_face_surf;

    /* Compute sum on rank */
    cs_array_reduce_wsum_l(n_loc_b_faces, /* Number of values */
                           1,             /* Dimension - 1 for scalar */
                           b_face_ids,    /* list of ids for values + weights */
                           NULL,          /* No ids for weights only */
                           scalar_vals,   /* Values to sum */
                           b_face_surf,   /* Array of weights */
                           &_l_sum);      /* local sum on rank */
  }

  /* Parallel sum */
  cs_parall_sum(1, CS_DOUBLE, &_l_sum);
  retval = _l_sum;

  /* Return value */
  return retval;
}

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
                                     const cs_lnum_t  b_face_ids[])
{
  cs_real_t retval = 0.;

  double _l_sum = 0.;

  if (n_loc_b_faces > 0) {
    const cs_real_t *restrict b_face_surf
      = cs_glob_mesh_quantities->b_face_surf;

    /* Compute sum on rank */
    cs_array_reduce_wsum_l(n_loc_b_faces, /* Number of values */
                           1,             /* Dimension - 1 for scalar */
                           NULL,          /* no ids list for values */
                           b_face_ids,    /* list of ids for weights */
                           scalar_vals,   /* Values to sum */
                           b_face_surf,   /* Array of weights */
                           &_l_sum);      /* local sum on rank */
  }

  /* Parallel sum */
  cs_parall_sum(1, CS_DOUBLE, &_l_sum);
  retval = _l_sum;

  /* Return value */
  return retval;
}

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
                             const cs_lnum_t  b_face_ids[])
{
  cs_real_t retval = 0.;

  double _l_sum  = 0.;
  double _l_wtot = 0.;

  if (n_loc_b_faces > 0) {
    const cs_real_t *restrict b_face_surf
      = cs_glob_mesh_quantities->b_face_surf;

    /* Compute sum on rank */
    cs_array_reduce_wsum_components_l(n_loc_b_faces, /* number of values */
                                      1,             /* dim. - 1 for scalar */
                                      b_face_ids,    /* list of ids for values
                                                        and weights */
                                      NULL,          /* no ids for weight only */
                                      scalar_vals,   /* values to sum */
                                      b_face_surf,   /* weights */
                                      &_l_sum,       /* local sum on rank */
                                      &_l_wtot);     /* local total weight */
  }

  /* Parallel sum */
  double _w[2] = {_l_sum, _l_wtot};
  cs_parall_sum(2, CS_DOUBLE, _w);

  retval = _w[0] / _w[1];

  /* Return value */
  return retval;
}

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
                                 const cs_lnum_t  b_face_ids[])
{
  cs_real_t retval = 0.;

  double _l_sum  = 0.;
  double _l_wtot = 0.;

  if (n_loc_b_faces > 0) {
    const cs_real_t *restrict b_face_surf
      = cs_glob_mesh_quantities->b_face_surf;

    /* Compute sum on rank */
    cs_array_reduce_wsum_components_l(n_loc_b_faces, /* number of values */
                                      1,             /* dim - 1 for scalar */
                                      NULL,          /* no ids list for values */
                                      b_face_ids,    /* ids for weights */
                                      scalar_vals,   /* values to sum */
                                      b_face_surf,   /* weights */
                                      &_l_sum,       /* local sum on rank */
                                      &_l_wtot);     /* local total weight */
  }

  /* Parallel sum */
  double _w[2] = {_l_sum, _l_wtot};
  cs_parall_sum(2, CS_DOUBLE, _w);

  retval = _w[0] / _w[1];

  /* Return value */
  return retval;
}

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
                               const cs_real_t *scalar_vals)
{
  return cs_post_scalar_boundary_integral(scalar_vals,
                                          z->n_elts,
                                          z->elt_ids);
}

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
                                   const cs_real_t *scalar_vals)
{
  return cs_post_bnd_scalar_boundary_integral(scalar_vals,
                                              z->n_elts,
                                              z->elt_ids);
}

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
                           const cs_real_t *scalar_vals)
{
  /* Since a boundary zone surface is known, no need to recompute
   * the total surface using the cs_array_reduce_wsum_components_l
   * function. We can directly divide the integral by the zone
   * surface (measure).
   */
  cs_real_t retval = cs_post_scalar_boundary_integral(scalar_vals,
                                                      z->n_elts,
                                                      z->elt_ids);
  retval /= z->measure;

  return retval;
}

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
                               const cs_real_t *scalar_vals)
{
  /* Since a boundary zone surface is known, no need to recompute
   * the total surface using the cs_array_reduce_wsum_components_l
   * function. We can directly divide the integral by the zone
   * surface (measure).
   */
  cs_real_t retval = cs_post_bnd_scalar_boundary_integral(scalar_vals,
                                                          z->n_elts,
                                                          z->elt_ids);
  retval /= z->measure;

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
