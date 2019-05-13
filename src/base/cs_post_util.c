/*============================================================================
 * Postprocessing utility functions.
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

#include "fvm_selector.h"

#include "cs_interface.h"

#include "cs_base.h"
#include "cs_balance_by_zone.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_geom.h"
#include "cs_gradient.h"
#include "cs_gradient.h"
#include "cs_gradient_perio.h"
#include "cs_join.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_math.h"
#include "cs_matrix_default.h"
#include "cs_mesh.h"
#include "cs_mesh_coherency.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_physical_constants.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_renumber.h"
#include "cs_rotation.h"
#include "cs_stokes_model.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"

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

/*! Status of post utilities */

int cs_glob_post_util_flag[CS_POST_UTIL_N_TYPES]
  = {-1, -1};

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
  cs_real_t *sx = (cs_real_t *)input;

  const cs_real_t sx0[3] = {sx[0], sx[1], sx[2]};
  const cs_real_t sx1[3] = {sx[3], sx[4], sx[5]};

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_lnum_t _n_cells = m->n_cells;
  cs_lnum_t *_cell_ids = NULL;

  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  BFT_MALLOC(_cell_ids, _n_cells, cs_lnum_t); /* Allocate selection list */

  /* Mark for each cell */
  /*--------------------*/

  for (cs_lnum_t cell_id = 0; cell_id < _n_cells; cell_id++) {
    _cell_ids[cell_id] = -1;
  }

  const cs_real_3_t *vtx_coord= (const cs_real_3_t *)m->vtx_coord;

  /* Contribution from interior faces;
     note the to mark cells, we could use a simple loop,
     as thread races would not lead to a incorrect result, but
     even if is slightly slower, we prefer to have a clean
     behavior under thread debuggers. */

  for (int g_id = 0; g_id < n_i_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {

      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        int n_crossings[2] = {0, 0};

        cs_lnum_t vtx_start = m->i_face_vtx_idx[face_id];
        cs_lnum_t vtx_end = m->i_face_vtx_idx[face_id+1];
        cs_lnum_t n_vertices = vtx_end - vtx_start;
        const cs_lnum_t *vertex_ids = m->i_face_vtx_lst + vtx_start;

        const cs_real_t *face_center = fvq->i_face_cog + (3*face_id);

        double t = cs_geom_segment_intersect_face(0,
                                                  n_vertices,
                                                  vertex_ids,
                                                  vtx_coord,
                                                  face_center,
                                                  sx0,
                                                  sx1,
                                                  n_crossings,
                                                  NULL);

        if (t >= 0 && t <= 1) {
          cs_lnum_t  c_id0 = m->i_face_cells[face_id][0];
          cs_lnum_t  c_id1 = m->i_face_cells[face_id][1];
          if (c_id0 < _n_cells)
            _cell_ids[c_id0] = 1;
          if (c_id1 < _n_cells)
            _cell_ids[c_id1] = 1;
        }

      }

    }

  }

  /* Contribution from boundary faces*/

  for (int g_id = 0; g_id < n_b_groups; g_id++) {

#   pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {

      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        int n_crossings[2] = {0, 0};

        cs_lnum_t vtx_start = m->b_face_vtx_idx[face_id];
        cs_lnum_t vtx_end = m->b_face_vtx_idx[face_id+1];
        cs_lnum_t n_vertices = vtx_end - vtx_start;
        const cs_lnum_t *vertex_ids = m->b_face_vtx_lst + vtx_start;

        const cs_real_t *face_center = fvq->b_face_cog + (3*face_id);

        double t = cs_geom_segment_intersect_face(0,
                                                  n_vertices,
                                                  vertex_ids,
                                                  vtx_coord,
                                                  face_center,
                                                  sx0,
                                                  sx1,
                                                  n_crossings,
                                                  NULL);

        if (t >= 0 && t <= 1) {
          cs_lnum_t  c_id = m->b_face_cells[face_id];
          _cell_ids[c_id] = 1;
        }

      }

    }

  }

  /* Now check marked cells */

  _n_cells = 0;
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    if (_cell_ids[cell_id] >= 0)
      _cell_ids[_n_cells++] = cell_id;
  }

  BFT_REALLOC(_cell_ids, _n_cells, cs_lnum_t); /* Adjust size (good practice,
                                                  but not required) */

  /* Set return values */

  *n_cells = _n_cells;
  *cell_ids = _cell_ids;
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

  cs_cell_segment_intersect_select(input, &n_cells, &cell_ids);

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
 *----------------------------------------------------------------------------*/

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
  const cs_real_3_t *surfbo =
    (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_normal;
  const cs_real_t *surfbn = cs_glob_mesh_quantities->b_face_surf;
  const cs_real_3_t *forbr =
    (const cs_real_3_t *)cs_field_by_name("boundary_forces")->val;
  cs_lnum_t ifac;
  cs_real_t srfbn, srfnor[3], fornor;

  for (cs_lnum_t iloc = 0 ; iloc < n_b_faces; iloc++) {
    ifac = b_face_ids[iloc];
    srfbn = surfbn[ifac];
    srfnor[0] = surfbo[ifac][0] / srfbn;
    srfnor[1] = surfbo[ifac][1] / srfbn;
    srfnor[2] = surfbo[ifac][2] / srfbn;
    fornor = forbr[ifac][0]*srfnor[0]
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
 * \param[in]   hyd_p_flag   flag for hydrostatic pressure
 * \param[in]   f_ext        exterior force generating
 *                           the hydrostatic pressure
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

  int hyd_p_flag = cs_glob_stokes_model->iphydr;
  cs_real_3_t *f_ext = (hyd_p_flag == 1) ?
    (cs_real_3_t *)cs_field_by_name_try("volume_forces"):NULL;

  bool use_previous_t = false;
  int inc = 1;
  int _recompute_cocg = 1;
  cs_field_gradient_potential(CS_F_(p),
                              use_previous_t,
                              inc,
                              _recompute_cocg,
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
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute Reynolds stresses in case of Eddy Viscosity Models
 *
 * \param[in]  n_cells     number of cells
 * \param[in]  cell_ids    list of cells (0 to n-1) containing given coordinates
 * \param[in]  coords      coordinates
 * \param[out] rst         Reynolds stresses stored as vector
 *                         [r11,r22,r33,r12,r23,r13]
 */
/*----------------------------------------------------------------------------*/

void
cs_post_evm_reynolds_stresses(cs_lnum_t          n_cells,
                              const cs_lnum_t    cell_ids[],
                              const cs_real_3_t *coords,
                              cs_real_6_t       *rst)
{
  const cs_turb_model_t *turb_model = cs_glob_turb_model;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_real_3_t *cell_cen =
    (cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;

  if (   turb_model->itytur != 2
      && turb_model->itytur != 6
      && turb_model->itytur != 5)
    bft_error(__FILE__, __LINE__, 0,
              _("This post-processing utility function is only available for "
                "Eddy Viscosity Models."));

  /* velocity gradient */

  cs_real_33_t *gradv;
  BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

  bool use_previous_t = false;
  int inc = 1;
  cs_field_gradient_vector(CS_F_(vel),
                           use_previous_t,
                           inc,
                           gradv);

  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);

  /* turbulent kinetic energy gradient for reconstruction */

  cs_real_3_t *gradk;
  if (var_cal_opt.ircflu > 0 && coords != NULL) {
    BFT_MALLOC(gradk, n_cells_ext, cs_real_3_t);
    bool recompute_cocg = true;
    cs_field_gradient_scalar(CS_F_(k),
                             use_previous_t,
                             inc,
                             recompute_cocg,
                             gradk);
  }

  /* Compute Reynolds stresses */

  const cs_real_t d2s3 = 2./3.;
  for (cs_lnum_t iloc = 0; iloc < n_cells; iloc++) {
    cs_lnum_t iel = cell_ids[iloc];

    cs_real_t divu = gradv[iel][0][0] + gradv[iel][1][1] + gradv[iel][2][2];
    cs_real_t nut = CS_F_(mu_t)->val[iel]/CS_F_(rho)->val[iel];
    cs_real_t xk = CS_F_(k)->val[iel];

    if (var_cal_opt.ircflu > 0 && coords != NULL) {
      for (int ii = 0; ii < 3; ii++) {
        xk += (coords[iloc][ii] - cell_cen[iel][ii])*gradk[iel][ii];
      }
    }

    cs_real_t xdiag = d2s3*(xk+ nut*divu);
    rst[iloc][0] =  xdiag - 2.*nut*gradv[iel][0][0];
    rst[iloc][1] =  xdiag - 2.*nut*gradv[iel][1][1];
    rst[iloc][2] =  xdiag - 2.*nut*gradv[iel][2][2];
    rst[iloc][3] = -nut*(gradv[iel][1][0]+gradv[iel][0][1]);
    rst[iloc][4] = -nut*(gradv[iel][2][1]+gradv[iel][1][2]);
    rst[iloc][5] = -nut*(gradv[iel][2][0]+gradv[iel][0][2]);
  }

  BFT_FREE(gradv);
  if (var_cal_opt.ircflu > 0 && coords != NULL) BFT_FREE(gradk);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the Q-criterion from Hunt et. al over each cell of a specified
 *        volume region.
 *
 * \f[
 *    Q = \tens{\Omega}:\tens{\Omega} -
 *    \deviator{ \left(\tens{S} \right)}:\deviator{ \left(\tens{S} \right)}
 * \f]
 * where \f$\tens{\Omega}\f$ is the vorticity tensor and
 * \f$\deviator{ \left(\tens{S} \right)}\f$ the deviatoric of the rate of strain
 * tensor.
 *
 * \param[in]  n_loc_cells  number of cells
 * \param[in]  cell_ids     list of cells (0 to n-1)
 * \param[out] q_crit       Q-criterion over the specified volume region.
 */
/*----------------------------------------------------------------------------*/

void
cs_post_q_criterion(const cs_lnum_t  n_loc_cells,
                    const cs_lnum_t  cell_ids[],
                    cs_real_t        q_crit[])
{
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  cs_real_33_t *gradv;

  BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

  bool use_previous_t = false;
  int inc = 1;
  cs_field_gradient_vector(CS_F_(vel),
                           use_previous_t,
                           inc,
                           gradv);

  for (cs_lnum_t i = 0; i < n_loc_cells; i++) {
    cs_lnum_t c_id = cell_ids[i];
    q_crit[i] = -1./6. * (   cs_math_sq(gradv[c_id][0][0])
                          +  cs_math_sq(gradv[c_id][1][1])
                          +  cs_math_sq(gradv[c_id][2][2]))
                - gradv[c_id][0][1]*gradv[c_id][1][0]
                - gradv[c_id][0][2]*gradv[c_id][2][0]
                - gradv[c_id][1][2]*gradv[c_id][2][1];
  }

  BFT_FREE(gradv);
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

END_C_DECLS
