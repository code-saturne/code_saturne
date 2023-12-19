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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal.h"
#include "fvm_nodal_extract.h"
#include "fvm_point_location.h"

#include "cs_base.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_boundary_zone.h"
#include "cs_coupling.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_file_csv_parser.h"
#include "cs_geom.h"
#include "cs_halo.h"
#include "cs_io.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_porous_model.h"
#include "cs_equation_iterative_solve.h"
#include "cs_physical_constants.h"
#include "cs_post.h"
#include "cs_timer.h"

#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_porosity_from_scan.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

static cs_porosity_from_scan_opt_t _porosity_from_scan_opt = {
  .compute_porosity_from_scan = false,
  .file_names = NULL,
  .output_name = NULL,
  .postprocess_points = true,
  .transformation_matrix = {{1., 0., 0., 0.},
                            {0., 1., 0., 0.},
                            {0., 0., 1., 0.}},
  .nb_sources = 0,
  .sources = NULL,
  .source_c_ids = NULL,
  .threshold = 4
};

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_porosity_from_scan_opt_t *cs_glob_porosity_from_scan_opt
= &_porosity_from_scan_opt;

static  ple_locator_t  *_locator = NULL;  /* PLE locator */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_porosity_from_scan_get_pointer(bool  **compute_porosity_from_scan);

/*============================================================================
 * Private function definitions
 *============================================================================*/
/*----------------------------------------------------------------------------
 * Incremental local solid plane computation at cells from scan points.
 * So that the summed squared distance to all points is minimized.
 *
 * parameters:
 *   m               <-- pointer to mesh
 *   n_points        <-- total number of points
 *   elt_ids         <-- point to cell id
 *   n_points_cell   <-- number of points in cell
 *   cen_cell        <-- center of gravity of cell
 *   cov_mat         --> incremental covariance matrix
 *----------------------------------------------------------------------------*/

static void
_incremental_solid_plane_from_points(const cs_mesh_t   *m,
                                     cs_lnum_t          n_points,
                                     const cs_lnum_t    elt_ids[],
                                     const cs_real_t    n_points_cell[],
                                     const cs_real_3_t  cen_cell[],
                                     const cs_real_3_t  point_coords[],
                                     cs_real_33_t       cov_mat[])
{
  cs_real_33_t  *c, *d, *z, *t;
  BFT_MALLOC(c, m->n_cells, cs_real_33_t);
  BFT_MALLOC(d, m->n_cells, cs_real_33_t);
  BFT_MALLOC(z, m->n_cells, cs_real_33_t);
  BFT_MALLOC(t, m->n_cells, cs_real_33_t);
  const cs_real_t threshold = _porosity_from_scan_opt.threshold;

  // Initialization

  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        c[c_id][i][j] = 0.;
        d[c_id][i][j] = 0.;
        t[c_id][i][j] = 0.;
        z[c_id][i][j] = 0.;
      }
    }
  }

  // Loop over points

  for (cs_lnum_t p_id = 0; p_id < n_points; p_id++) {

    cs_lnum_t cell_id = elt_ids[p_id];

    if (n_points_cell[cell_id] > threshold) { // At least three points required

      cs_real_t point_local[3];
      for (cs_lnum_t i = 0; i < 3; i++)
        point_local[i] = point_coords[p_id][i] - cen_cell[cell_id][i];

      // Kahan summation
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          z[cell_id][i][j] = point_local[i] * point_local[j]
                            - c[cell_id][i][j];
          t[cell_id][i][j] = d[cell_id][i][j] + z[cell_id][i][j];
          c[cell_id][i][j] = (t[cell_id][i][j] - d[cell_id][i][j])
                            - z[cell_id][i][j];
          d[cell_id][i][j] = t[cell_id][i][j];
        }
      }
    }

  } // Loop over points

  BFT_FREE(z);
  BFT_FREE(t);
  BFT_FREE(c);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    for (cs_lnum_t i = 0; i < 3; i++)
      for (cs_lnum_t j = 0; j < 3; j++)
          cov_mat[cell_id][i][j] += d[cell_id][i][j];

  } // Loop over cells

  BFT_FREE(d);
}

/*----------------------------------------------------------------------------
 * Compute local solid planes at cells from scan points file.
 * So that the summed squared distance to all points is minimized.
 *
 * parameters:
 *   m               <-- pointer to mesh
 *   n_points_cell   <-- number of points in cell
 *   cen_points      <-- center of gravity of points relative to cell center
 *   cov_mat         <-- Covariance matrix of points in cell
 *   c_w_face_normal --> normal vector to the solid plane
 *----------------------------------------------------------------------------*/

static void
_solid_plane_from_points(const cs_mesh_t   *m,
                         const cs_real_t    n_points_cell[],
                         const cs_real_3_t  cen_points[],
                         const cs_real_33_t cov_mat[],
                         cs_real_3_t        c_w_face_normal[])
{
  const cs_real_t threshold = _porosity_from_scan_opt.threshold;

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {

    /* At least three points required */
    if (n_points_cell[cell_id] > threshold) {

      cs_real_33_t cv;
      cs_real_3_t sx;

      cs_real_33_t a_x, a_y, a_z;
      cs_real_33_t b;

      // Initialization
      for (cs_lnum_t i = 0; i < 3; i++) {
        sx[i] = cen_points[cell_id][i]/n_points_cell[cell_id];
        for (cs_lnum_t j = 0; j < 3; j++) {
          cv[i][j] = cov_mat[cell_id][i][j]/n_points_cell[cell_id];
        }
      }

      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          a_x[i][j] = cv[i][j];
          a_y[i][j] = cv[i][j];
          a_z[i][j] = cv[i][j];
        }
      }
      a_x[0][0] = 1.;
      a_x[0][1] = sx[1];
      a_x[0][2] = sx[2];
      a_x[1][0] = sx[1];
      a_x[2][0] = sx[2];

      a_y[0][1] = sx[0];
      a_y[1][0] = sx[0];
      a_y[1][1] = 1.;
      a_y[1][2] = sx[2];
      a_y[2][1] = sx[2];

      a_z[0][2] = sx[0];
      a_z[1][2] = sx[1];
      a_z[2][0] = sx[0];
      a_z[2][1] = sx[1];
      a_z[2][2] = 1.;


      cs_real_t det_x =  cs_math_33_determinant(a_x);
      cs_real_t det_y =  cs_math_33_determinant(a_y);
      cs_real_t det_z =  cs_math_33_determinant(a_z);

      /* RHS of linear system */
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          b[i][j] = -cv[i][j];
          if (i == j)
            b[i][j] = -sx[i];
        }
      }

      cs_real_t det_max = fmax(fmax(det_x, det_y), det_z);
      if (det_max > 0.0) {
        // Pick path with best conditioning:
        if (CS_ABS(det_max - det_x) < cs_math_epzero*det_max) {
          cs_math_33_inv_cramer_in_place(a_x);
          cs_math_33_3_product(a_x, b[0], c_w_face_normal[cell_id]);
          c_w_face_normal[cell_id][0] = 1.;
        }
        else if (CS_ABS(det_max - det_y) < cs_math_epzero*det_max) {
          cs_math_33_inv_cramer_in_place(a_y);
          cs_math_33_3_product(a_y, b[1], c_w_face_normal[cell_id]);
          c_w_face_normal[cell_id][1] = 1.;
        }
        else {
          cs_math_33_inv_cramer_in_place(a_z);
          cs_math_33_3_product(a_z, b[2], c_w_face_normal[cell_id]);
          c_w_face_normal[cell_id][2] = 1.;
        }
        cs_math_3_normalize(c_w_face_normal[cell_id],
                            c_w_face_normal[cell_id]);
      }

      // c_w_face_normal is forced to be 0 vector
      else {
        c_w_face_normal[cell_id][0] = 0.;
        c_w_face_normal[cell_id][1] = 0.;
        c_w_face_normal[cell_id][2] = 0.;
      }
    }

    // If not enough points, c_w_face_normal is forced to be 0 vector
    else {
      c_w_face_normal[cell_id][0] = 0.;
      c_w_face_normal[cell_id][1] = 0.;
      c_w_face_normal[cell_id][2] = 0.;
    }

  } // Loop over cells

}
/*----------------------------------------------------------------------------
 * Prepare computation of porosity from scan points file.
 *
 * It read the points file, count the point per cells and penalize cell
 * with points
 *
 * parameters:
 *   m  <-- pointer to mesh
 *   mq <-- pointer to mesh quantities
 *----------------------------------------------------------------------------*/

static void
_prepare_porosity_from_scan(const cs_mesh_t             *m,
                            const cs_mesh_quantities_t  *mq) {
  char line[512];

  cs_real_t *restrict cell_f_vol = mq->cell_f_vol;

  /* Open file */
  bft_printf(_("\n\n  Compute the porosity from a scan points file:\n"
               "    %s\n\n"),
             _porosity_from_scan_opt.file_names);

  bft_printf(_("  Transformation       %12.5g %12.5g %12.5g %12.5g\n"
               "  matrix:              %12.5g %12.5g %12.5g %12.5g\n"
               "                       %12.5g %12.5g %12.5g %12.5g\n"
               "    (last column is translation vector)\n\n"),
             _porosity_from_scan_opt.transformation_matrix[0][0],
             _porosity_from_scan_opt.transformation_matrix[0][1],
             _porosity_from_scan_opt.transformation_matrix[0][2],
             _porosity_from_scan_opt.transformation_matrix[0][3],
             _porosity_from_scan_opt.transformation_matrix[1][0],
             _porosity_from_scan_opt.transformation_matrix[1][1],
             _porosity_from_scan_opt.transformation_matrix[1][2],
             _porosity_from_scan_opt.transformation_matrix[1][3],
             _porosity_from_scan_opt.transformation_matrix[2][0],
             _porosity_from_scan_opt.transformation_matrix[2][1],
             _porosity_from_scan_opt.transformation_matrix[2][2],
             _porosity_from_scan_opt.transformation_matrix[2][3]);


  /* Bounding box */
  cs_real_t min_vec_tot[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
  cs_real_t max_vec_tot[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};

  /* Points local center of gravity */
  cs_real_t *cen_points =
    (cs_real_t *)cs_field_by_name("cell_scan_points_cog")->val;

  /* Points local averaged color */
  cs_real_t *cell_color =
    (cs_real_t *)cs_field_by_name("cell_scan_points_color")->val;

  /* Immersed solid roughness */
  cs_real_t *c_w_face_rough =
    (cs_real_t *)cs_field_by_name("solid_roughness")->val;

  /* Covariance matrix for solid plane computation */
  cs_real_33_t *cov_mat;
  BFT_MALLOC(cov_mat, m->n_cells, cs_real_33_t);
  memset(cov_mat, 0., m->n_cells * sizeof(cs_real_33_t));

  /* Loop on file_names */
  char *tok;
  const char sep[4] = ";";
  // parse names

  char *file_names;
  BFT_MALLOC(file_names, strlen(_porosity_from_scan_opt.file_names)+1, char);
  strcpy(file_names, _porosity_from_scan_opt.file_names);

  tok = strtok(file_names, sep);

  while (tok != NULL) {
    char *f_name;
    BFT_MALLOC(f_name,
        strlen(tok) + 1,
        char);
    strcpy(f_name, tok);

    bft_printf(_("\n\n  Open file:\n"
          "    %s\n\n"),
        f_name);
    FILE *file = fopen(f_name, "rt");
    if (file == NULL)
      bft_error(__FILE__,__LINE__, 0,
          _("Porosity from scan: Could not open file."));

    /* next file to be read */
    tok = strtok(NULL, sep);
    long int n_read_points = 0;
    long int n_points = 0;
    if (fscanf(file, "%ld\n", &n_read_points) != 1)
      bft_error(__FILE__,__LINE__, 0,
                _("Porosity from scan: Could not read the number of lines."));

    bft_printf(_("  Porosity from scan: %ld points to be read.\n\n"),
               n_read_points);

    /* Pointer to field */
    cs_field_t *f_nb_scan = cs_field_by_name_try("nb_scan_points");

    /* Location mesh where points will be localized */
    fvm_nodal_t *location_mesh =
      cs_mesh_connect_cells_to_nodal(m,
                                     "pts_location_mesh",
                                     false, // no family info
                                     m->n_cells,
                                     NULL);

    fvm_nodal_make_vertices_private(location_mesh);

    /* Read multiple scan files
     * ------------------------ */

    for (int n_scan = 0; n_read_points > 0; n_scan++) {
      n_points = n_read_points;
      bft_printf(_("  Immersed boundary from scan: scan %d.\n"
                   "                               %ld points to be read.\n\n"),
                 n_scan, n_points);


      cs_real_3_t *point_coords;
      float *colors;
      BFT_MALLOC(point_coords, n_points, cs_real_3_t);
      BFT_MALLOC(colors, 3*n_points, float);

      cs_real_3_t min_vec = { HUGE_VAL,  HUGE_VAL,  HUGE_VAL};
      cs_real_3_t max_vec = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};

      /* Read points */
      for (int i = 0; i < n_points; i++ ) {
        int num, green, red, blue;
        cs_real_4_t xyz;
        for (int j = 0; j < 3; j++)
          point_coords[i][j] = 0.;

        if (fscanf(file, "%lf", &(xyz[0])) != 1)
          bft_error
            (__FILE__,__LINE__, 0,
             _("Porosity from scan: Error while reading dataset. Line %d\n"), i);
        if (fscanf(file, "%lf", &(xyz[1])) != 1)
          bft_error
            (__FILE__,__LINE__, 0,
             _("Porosity from scan: Error while reading dataset."));
        if (fscanf(file, "%lf", &(xyz[2])) != 1)
          bft_error(__FILE__,__LINE__, 0,
                    _("Porosity from scan: Error while reading dataset."));

        /* Translation and rotation */
        xyz[3] = 1.;
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 4; k++)
            point_coords[i][j]
              += _porosity_from_scan_opt.transformation_matrix[j][k] * xyz[k];

          /* Compute bounding box*/
          min_vec[j] = CS_MIN(min_vec[j], point_coords[i][j]);
          max_vec[j] = CS_MAX(max_vec[j], point_coords[i][j]);
        }

        /* Intensities */
        if (fscanf(file, "%d", &num) != 1)
          bft_error(__FILE__,__LINE__, 0,
                    _("Porosity from scan: Error while reading dataset."));

        /* Red */
        if (fscanf(file, "%d", &red) != 1)
          bft_error(__FILE__,__LINE__, 0,
                    _("Porosity from scan: Error while reading dataset (red). "
                      "npoints read %d\n"), i);
        /* Green */
        if (fscanf(file, "%d", &green) != 1)
          bft_error(__FILE__,__LINE__, 0,
                    _("Porosity from scan: Error while reading dataset (green). "
                      "npoints read %d\n"), i);
        /* Blue */
        if (fscanf(file, "%d\n", &blue) != 1)
          bft_error(__FILE__,__LINE__, 0,
                    _("Porosity from scan: Error while reading dataset (blue). "
                      "npoints read %d\n"), i);

        /* When colors are written as int, Paraview interprets them in [0, 255]
         * when they are written as float, Paraview interprets them in [0., 1.]
         * */
        colors[3*i + 0] = red/255.;
        colors[3*i + 1] = green/255.;
        colors[3*i + 2] = blue/255.;
      }

      /* Check EOF was correctly reached */
      if (fgets(line, sizeof(line), file) != NULL)
        n_read_points = strtol(line, NULL, 10);
      else
        n_read_points = 0;

      /* Bounding box*/
      bft_printf(_("  Bounding box [%f, %f, %f], [%f, %f, %f].\n\n"),
                 min_vec[0], min_vec[1], min_vec[2],
                 max_vec[0], max_vec[1], max_vec[2]);

      /* Update global bounding box */
      for (int j = 0; j < 3; j++) {
        min_vec_tot[j] = CS_MIN(min_vec[j], min_vec_tot[j]);
        max_vec_tot[j] = CS_MAX(max_vec[j], max_vec_tot[j]);
      }

      if (n_read_points > 0)
        bft_printf
          (_("  Porosity from scan: %ld additional points to be read.\n\n"),
           n_read_points);

      /* FVM meshes for writers */
      if (_porosity_from_scan_opt.postprocess_points) {
        char *fvm_name;
        if (_porosity_from_scan_opt.output_name == NULL) {
          BFT_MALLOC(fvm_name,
                     strlen(f_name) + 3 + 1,
                     char);
          strcpy(fvm_name, f_name);
        }
        else {
          BFT_MALLOC(fvm_name,
                     strlen(_porosity_from_scan_opt.output_name) + 3 + 1,
                     char);
          strcpy(fvm_name, _porosity_from_scan_opt.output_name);
        }
        char suffix[13];
        sprintf(suffix, "_%02d", n_scan);
        strcat(fvm_name, suffix);

        /* Build FVM mesh from scanned points */
        fvm_nodal_t *pts_mesh = fvm_nodal_create(fvm_name, 3);

        /* Only the first rank writes points for now */
        cs_gnum_t *vtx_gnum = NULL;

        if (cs_glob_rank_id < 1) {
          /* Update the points set structure */
          fvm_nodal_define_vertex_list(pts_mesh, n_points, NULL);
          fvm_nodal_set_shared_vertices(pts_mesh, (cs_coord_t *)point_coords);

          BFT_MALLOC(vtx_gnum, n_points, cs_gnum_t);
          for (cs_lnum_t i = 0; i < n_points; i++)
            vtx_gnum[i] = i + 1;

        }
        fvm_nodal_init_io_num(pts_mesh, vtx_gnum, 0);

        /* Free if allocated */
        BFT_FREE(vtx_gnum);

        /* Create default writer */
        fvm_writer_t *writer
          = fvm_writer_init(fvm_name,
                            "postprocessing",
                            cs_post_get_default_format(),
                            cs_post_get_default_format_options(),
                            FVM_WRITER_FIXED_MESH);

        fvm_writer_export_nodal(writer, pts_mesh);

        const void *var_ptr[1] = {NULL};

        var_ptr[0] = colors;

        fvm_writer_export_field(writer,
                                pts_mesh,
                                "color",
                                FVM_WRITER_PER_NODE,
                                3,
                                CS_INTERLACE,
                                0,
                                0,
                                CS_FLOAT,
                                -1,
                                0.0,
                                (const void * *)var_ptr);

        /* Free and destroy */
        fvm_writer_finalize(writer);
        pts_mesh = fvm_nodal_destroy(pts_mesh);
        BFT_FREE(fvm_name);
        BFT_FREE(f_name);
      }

      /* Now build locator
       * Locate points on this location mesh */
      /*-------------------------------------*/

      int options[PLE_LOCATOR_N_OPTIONS];
      for (int i = 0; i < PLE_LOCATOR_N_OPTIONS; i++)
        options[i] = 0;
      options[PLE_LOCATOR_NUMBERING] = 0; /* base 0 numbering */

#if defined(PLE_HAVE_MPI)
      _locator = ple_locator_create(cs_glob_mpi_comm,
                                    cs_glob_n_ranks,
                                    0);
#else
      _locator = ple_locator_create();
#endif

      cs_lnum_t _n_points = (cs_glob_rank_id < 1) ? n_points : 0;
      cs_real_t *_point_coords = (cs_glob_rank_id < 1) ? (cs_real_t *)point_coords : NULL;

      ple_locator_set_mesh(_locator,
                           location_mesh,
                           options,
                           0., /* tolerance_base */
                           0.1, /* tolerance */
                           3, /* dim */
                           _n_points,
                           NULL,
                           NULL, /* point_tag */
                           _point_coords,
                           NULL, /* distance */
                           cs_coupling_mesh_extents,
                           cs_coupling_point_in_mesh_p);

      /* Shift from 1-base to 0-based locations */
      ple_locator_shift_locations(_locator, -1);

      /* dump locator */
#if 0
      ple_locator_dump(_locator);
#endif

      /* Number of distant points located on local mesh. */
      cs_lnum_t n_points_dist = ple_locator_get_n_dist_points(_locator);

#if 0
      bft_printf("ple_locator_get_n_dist_points = %d, n_points = %ld\n",
                 n_points_dist, n_points);
#endif

      const cs_lnum_t *dist_loc = ple_locator_get_dist_locations(_locator);
      const ple_coord_t *dist_coords = ple_locator_get_dist_coords(_locator);

      float *dist_colors = NULL;
      BFT_MALLOC(dist_colors, 3*n_points_dist, float);

      ple_locator_exchange_point_var(_locator,
                                     dist_colors,
                                     colors,
                                     NULL,
                                     sizeof(float),
                                     3,
                                     1);

      for (cs_lnum_t i = 0; i < n_points_dist; i++) {
        cs_lnum_t c_id = dist_loc[i];
        f_nb_scan->val[c_id] += 1.;
        for (cs_lnum_t idim = 0; idim < 3; idim++) {
          cen_points[c_id*3+idim] += (dist_coords[i*3 + idim]
                                    - mq->cell_cen[c_id*3+idim]);
          cell_color[c_id*3+idim] += dist_colors[i*3 + idim];
        }
      }

      BFT_FREE(dist_colors);

      _incremental_solid_plane_from_points(m,
                                           n_points_dist,
                                           dist_loc,
                                           (const cs_real_t   *)f_nb_scan->val,
                                           (const cs_real_3_t *)mq->cell_cen,
                                           (const cs_real_3_t *)dist_coords,
                                           cov_mat);

      /* Solid face roughness from point cloud is computed as the RMS of points
       *  distance to the reconstructed plane */
      cs_real_t vec_w_point[3]  = {0., 0., 0.};
      cs_real_t w_point_dist =  0.;

      for (cs_lnum_t i = 0; i < n_points_dist; i++) {
        cs_lnum_t c_id = dist_loc[i];

        if (f_nb_scan->val[c_id] > 1.)  { // at least 2 points to compute distance

          for (cs_lnum_t idim = 0; idim < 3; idim++)
            vec_w_point[idim] = dist_coords[i*3 + idim]
                              - cen_points[c_id*3+idim];
          //TODO compute roughness incrementally
          //w_point_dist = cs_math_3_dot_product(vec_w_point, c_w_face_normal[c_id]);

          c_w_face_rough[c_id] += 2 * sqrt(w_point_dist * w_point_dist)
                                / f_nb_scan->val[c_id];
        }
      }
      //TODO compute the minimum distance between point to suggest a minimum resolution

      /* Free memory */
      _locator = ple_locator_destroy(_locator);
      BFT_FREE(point_coords);
      BFT_FREE(colors);

    } /* End loop on multiple scans */

    if (fclose(file) != 0)
      bft_error(__FILE__,__LINE__, 0,
                _("Porosity from scan: Could not close the file."));

    /* Nodal mesh is not needed anymore */
    location_mesh = fvm_nodal_destroy(location_mesh);

  } /* End of multiple files */

  /* Finalization */

  const cs_field_t *f_nb_scan = cs_field_by_name_try("nb_scan_points");

  // Normal vector to the solid plane
  cs_real_3_t *restrict c_w_face_normal
    = (cs_real_3_t *restrict)mq->c_w_face_normal;

  _solid_plane_from_points(m,
                           (const cs_real_t   *)f_nb_scan->val,
                           (const cs_real_3_t *)cen_points,
                           cov_mat,
                           c_w_face_normal);

  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {
    if (f_nb_scan->val[c_id] > 0) {
      for (cs_lnum_t idim = 0; idim < 3; idim++) {
        cen_points[c_id*3+idim] /= f_nb_scan->val[c_id];
        cen_points[c_id*3+idim] += mq->cell_cen[c_id*3+idim];
        cell_color[c_id*3+idim] /= f_nb_scan->val[c_id];
      }
    }
  }

  /* Free memory */
  BFT_FREE(cov_mat);
  BFT_FREE(file_names);

  /* Parallel synchronisation */
  cs_mesh_sync_var_scal(mq->cell_vol);
  if (m->halo != NULL) {
    cs_halo_sync_var_strided(m->halo, CS_HALO_EXTENDED,
                             (cs_real_t *)cen_points, 3);
    cs_halo_sync_var_strided(m->halo, CS_HALO_EXTENDED,
                             (cs_real_t *)cell_color, 3);
    cs_halo_sync_var_strided(m->halo, CS_HALO_EXTENDED,
                             (cs_real_t *)c_w_face_normal, 3);
    if (m->n_init_perio > 0)
      cs_halo_perio_sync_coords(m->halo, CS_HALO_EXTENDED,
                                (cs_real_t *)cen_points);
    cs_parall_min(3, CS_REAL_TYPE, min_vec_tot);
    cs_parall_max(3, CS_REAL_TYPE, max_vec_tot);
  }

  /* Bounding box */
  bft_printf(_("  Global bounding box [%f, %f, %f], [%f, %f, %f].\n\n"),
             min_vec_tot[0], min_vec_tot[1], min_vec_tot[2],
             max_vec_tot[0], max_vec_tot[1], max_vec_tot[2]);


  /* Solid cells should have enough points */
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++) {
    cell_f_vol[cell_id] = mq->cell_vol[cell_id];
    if (cs_math_3_norm(c_w_face_normal[cell_id]) > 0.) {
      cell_f_vol[cell_id] = 0.;
      mq->c_disable_flag[cell_id] = 1;
    }
  }

  /* Parallel synchronisation */
  cs_mesh_sync_var_scal(cell_f_vol);

  cs_real_3_t *restrict i_face_normal
    =  (cs_real_3_t *restrict)mq->i_face_normal;
  cs_real_3_t *restrict b_face_normal
    =  (cs_real_3_t *restrict)mq->b_face_normal;
  cs_real_t *restrict i_face_surf
    = (cs_real_t *restrict)mq->i_face_surf;
  cs_real_t *restrict b_face_surf
    = (cs_real_t *restrict)mq->b_face_surf;

  cs_real_3_t *restrict i_f_face_normal
    = (cs_real_3_t *restrict)mq->i_f_face_normal;
  cs_real_3_t *restrict b_f_face_normal
    = (cs_real_3_t *restrict)mq->b_f_face_normal;
  cs_real_t *restrict i_f_face_surf
    = (cs_real_t *restrict)mq->i_f_face_surf;
  cs_real_t *restrict b_f_face_surf
    = (cs_real_t *restrict)mq->b_f_face_surf;

  /* Penalization of cells with points */

  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
    cs_lnum_t cell_id1 = m->i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = m->i_face_cells[face_id][1];
    if (cell_f_vol[cell_id1] <= 0. || cell_f_vol[cell_id2] <= 0.) {
      i_f_face_normal[face_id][0] = 0.;
      i_f_face_normal[face_id][1] = 0.;
      i_f_face_normal[face_id][2] = 0.;
      i_f_face_surf[face_id] = 0.;
    }
    else {
      i_f_face_normal[face_id][0] = i_face_normal[face_id][0];
      i_f_face_normal[face_id][1] = i_face_normal[face_id][1];
      i_f_face_normal[face_id][2] = i_face_normal[face_id][2];
      i_f_face_surf[face_id]      = i_face_surf[face_id]     ;
    }
  }

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
    cs_lnum_t cell_id = m->b_face_cells[face_id];
    if (cell_f_vol[cell_id] <= 0.) {
      b_f_face_normal[face_id][0] = 0.;
      b_f_face_normal[face_id][1] = 0.;
      b_f_face_normal[face_id][2] = 0.;
      b_f_face_surf[face_id] = 0.;
    }
    else {
      b_f_face_normal[face_id][0] = b_face_normal[face_id][0];
      b_f_face_normal[face_id][1] = b_face_normal[face_id][1];
      b_f_face_normal[face_id][2] = b_face_normal[face_id][2];
      b_f_face_surf[face_id]      = b_face_surf[face_id];
    }
  }
}

/*----------------------------------------------------------------------------
 * Get pointer
 *----------------------------------------------------------------------------*/

void
cs_f_porosity_from_scan_get_pointer(bool **compute_porosity_from_scan)
{
  *compute_porosity_from_scan
    = &(_porosity_from_scan_opt.compute_porosity_from_scan);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the file name of points for the computation of the
 *        porosity from scan.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_set_file_name(const char  *file_name)
{
  if (file_name == NULL) {
    _porosity_from_scan_opt.compute_porosity_from_scan = false;
    return;
  }

  _porosity_from_scan_opt.compute_porosity_from_scan = true;
  /* Force porous model */
  cs_glob_porous_model = 3;

  if (_porosity_from_scan_opt.file_names == NULL) {
    BFT_MALLOC(_porosity_from_scan_opt.file_names,
               strlen(file_name) + 1 + 1,
               char);
    sprintf(_porosity_from_scan_opt.file_names, "%s;", file_name);
  }
  else {
    int length = strlen(_porosity_from_scan_opt.file_names);
    BFT_REALLOC(_porosity_from_scan_opt.file_names,
               length + strlen(file_name) + 1 + 1,
               char);
    sprintf(_porosity_from_scan_opt.file_names, "%s%s;",
                    _porosity_from_scan_opt.file_names,
                    file_name);
  }

  bft_printf("Add file %s to the list %s\n",
                  file_name, _porosity_from_scan_opt.file_names);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the output name for the FVM writer of scan points.
 *
 * \param[in] output_name  name of the output (a suffix will be added)
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_set_output_name(const char  *output_name)
{
  if (output_name == NULL) {
    _porosity_from_scan_opt.postprocess_points = false;
    return;
  }

  _porosity_from_scan_opt.postprocess_points = true;

  BFT_MALLOC(_porosity_from_scan_opt.output_name,
             strlen(output_name) + 1,
             char);

  sprintf(_porosity_from_scan_opt.output_name, "%s", output_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a scanner source point.
 *
 * \param[in] source     source vector
 * \param[in] transform  flag to apply the transformation matrix to the source
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_add_source(const cs_real_t  source[3],
                                 bool             transform)
{
  /* Add a source */
  const int s_id = _porosity_from_scan_opt.nb_sources;
  _porosity_from_scan_opt.nb_sources++;

  BFT_REALLOC(_porosity_from_scan_opt.source_c_ids,
              _porosity_from_scan_opt.nb_sources,
              cs_lnum_t);

  BFT_REALLOC(_porosity_from_scan_opt.sources,
              _porosity_from_scan_opt.nb_sources,
              cs_real_3_t);

  if (transform) {
    /* Apply translation and rotation */
    for (int i = 0; i < 3; i++) {
      _porosity_from_scan_opt.sources[s_id][i] = 0;
      for (int j = 0; j < 3; j++)
        _porosity_from_scan_opt.sources[s_id][i]
          += _porosity_from_scan_opt.transformation_matrix[i][j] * source[j];
      _porosity_from_scan_opt.sources[s_id][i]
        += _porosity_from_scan_opt.transformation_matrix[i][3];
    }
  }
  else {
    for (int i = 0; i < 3; i++)
      _porosity_from_scan_opt.sources[s_id][i] = source[i];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the scanner sources from csv file to fill fluid space.
 *
 * \param[in] csv file containing the (x,y,z) coordinates of each scanner
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_add_sources_by_file_name(const char *file_name)
{
  if (file_name == NULL)
    bft_error(__FILE__,__LINE__, 0, _("Could not read scanner sources file"));

  /* read the csv file */
  const int s_col_idx[3] = {0, 1, 2}; /* columns to read */
  int nb_scan = 0, nb_cols = 0;
  char ***csv_data = cs_file_csv_parse(file_name,
                                       ",", /* separator */
                                       0, /* n_headers */
                                       3, /* 3 columns to read */
                                       s_col_idx,
                                       true, /* ignore_missing_tokens */
                                       &nb_scan,
                                       &nb_cols);

  /* loop on the scanner sources */
  for (int i = 0; i < nb_scan; i++) {
    const char *origin_x = csv_data[i][0];
    const char *origin_y = csv_data[i][1];
    const char *origin_z = csv_data[i][2];

    /* char to double conversion */
    cs_real_3_t source = {atof(origin_x), atof(origin_y), atof(origin_z)};

    /* Add source */
    bool transform = true;
    cs_porosity_from_scan_add_source(source, transform);
  }
  // Free data which is no longer needed.
  for (int i = 0; i < nb_scan; i++) {
    for (int j = 0; j < nb_cols; j++)
      BFT_FREE(csv_data[i][j]);
    BFT_FREE(csv_data[i]);
  }
  BFT_FREE(csv_data);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes the porosity which is equal to one from
 *        a source, radiating sphericaly, and is 0 when touching points
 *        of the scan.
 *
 *  This function solves the following transport equation on \f$ \varia \f$:
 *  \f[
 *  \dfrac{\partial \varia}{\partial t} + \divs \left( \varia \vect{e}_r \right)
 *      - \divs \left( \vect{e}_r \right) \varia = 0
 *  \f]
 *  where \f$ \vect{e}_r
 *          = \dfrac{\vect{x} - \vect{x}_0}{\norm{\vect{x} - \vect{x}_0}} \f$
 *  is the radial direction from the source \f$\vect{x}_0 \f$.
 *
 *  The boundary conditions on \f$ \varia \f$ is an homogeneous Neumann, and
 *  a penalisation term is impose in the cell of center \f$ \vect{x}_0\f$.
 *  \f[
 *   \dfrac{\partial \varia}{\partial n} = 0 \textrm{everywhere}
 *  \f]
 *
 *  Remarks:
 *  - a steady state is looked for.
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_porosity_from_scan(void)
{
  /* Initialization */

  const cs_domain_t *domain = cs_glob_domain;
  const cs_mesh_t *m = domain->mesh;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;
  cs_real_t *restrict cell_f_vol = mq->cell_f_vol;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)mq->cell_cen;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)mq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)mq->b_face_cog;
  cs_real_3_t *restrict i_f_face_normal =
     (cs_real_3_t *restrict)mq->i_f_face_normal;
  cs_real_3_t *restrict b_f_face_normal =
     (cs_real_3_t *restrict)mq->b_f_face_normal;
  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *b_face_cells
    = (const cs_lnum_t *)m->b_face_cells;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)mq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)mq->b_face_normal;
  cs_real_t *restrict i_f_face_surf
    = (cs_real_t *restrict)mq->i_f_face_surf;
  cs_real_t *restrict b_f_face_surf
    = (cs_real_t *restrict)mq->b_f_face_surf;

  /* Pointer to porosity field */
  cs_field_t *f = cs_field_by_name_try("porosity");

  cs_equation_param_t *eqp = cs_field_get_equation_param(f);
  eqp->verbosity = 2;
  cs_lnum_t *source_c_ids = _porosity_from_scan_opt.source_c_ids;
  int nb_sources = _porosity_from_scan_opt.nb_sources;

  /* First pass to put fluid surfaces to 0 for all faces' cells with
   * at least some points */
  _prepare_porosity_from_scan(m, mq);

  bft_printf(_("  Positions of the %d given sources:\n"), nb_sources);
  const cs_real_3_t *s_pos = _porosity_from_scan_opt.sources;
  for (int s_id = 0; s_id < nb_sources; s_id++) {
    bft_printf(_("   - Source %3d: [ %12.5e, %12.5e, %12.5e ]\n"),
              s_id, s_pos[s_id][0], s_pos[s_id][1], s_pos[s_id][2]);
  } /* Loop on s_id */
  bft_printf(_("\n"));

  cs_real_t *restrict i_massflux
    = cs_field_by_id
        (cs_field_get_key_int(f, cs_field_key_id("inner_mass_flux_id")))->val;
  cs_real_t *restrict b_massflux
    = cs_field_by_id
        (cs_field_get_key_int(f, cs_field_key_id("boundary_mass_flux_id")))->val;

  cs_var_cal_opt_t vcopt;
  cs_field_get_key_struct(f, cs_field_key_id("var_cal_opt"), &vcopt);

  /* Local variables */
  cs_real_t *rovsdt, *pvar, *dpvar, *rhs;
  BFT_MALLOC(rovsdt, m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(pvar, m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(dpvar, m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(rhs, m->n_cells_with_ghosts, cs_real_t);

  /* Loop over sources */
  for (int s_id = 0; s_id < nb_sources; s_id++) {

    int rank_source;
    cs_geom_closest_point(m->n_cells,
                          (const cs_real_3_t *)mq->cell_cen,
                          _porosity_from_scan_opt.sources[s_id],
                          &(source_c_ids[s_id]),
                          &rank_source);

    cs_real_t source_cen[3];

    if (source_c_ids[s_id] > -1) {
      for (int i = 0; i < 3; i++)
        source_cen[i] = cell_cen[source_c_ids[s_id]][i];
    }

    cs_parall_bcast(rank_source, 3, CS_REAL_TYPE, source_cen);

    /* Compute the mass flux due to V = e_r
     *=====================================*/

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      cs_real_t x0xf[3] = {
        i_face_cog[face_id][0] - source_cen[0],
        i_face_cog[face_id][1] - source_cen[1],
        i_face_cog[face_id][2] - source_cen[2]
      };
      i_massflux[face_id] = cs_math_3_dot_product(x0xf,
                                                  i_f_face_normal[face_id]);
    }

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
      cs_real_t x0xf[3] = {
        b_face_cog[face_id][0] - source_cen[0],
        b_face_cog[face_id][1] - source_cen[1],
        b_face_cog[face_id][2] - source_cen[2]
      };

      b_massflux[face_id] = cs_math_3_dot_product(x0xf,
                                                  b_f_face_normal[face_id]);
    }

    /* Boundary conditions
     *====================*/

    /* homogeneous Neumann except if ingoing flux */
    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
      /* Dirichlet BCs */
      if (b_massflux[face_id] <= 0.) {

       vcopt.ndircl = 1;
       cs_real_t hint = 1. / mq->b_dist[face_id];
       cs_real_t pimp = 0.;

       cs_boundary_conditions_set_dirichlet_scalar(&(f->bc_coeffs->a[face_id]),
                                                   &(f->bc_coeffs->af[face_id]),
                                                   &(f->bc_coeffs->b[face_id]),
                                                   &(f->bc_coeffs->bf[face_id]),
                                                   pimp,
                                                   hint,
                                                   cs_math_infinite_r);
      }
      else {

        cs_real_t hint = 1. / mq->b_dist[face_id];
        cs_real_t qimp = 0.;

        cs_boundary_conditions_set_neumann_scalar(&(f->bc_coeffs->a[face_id]),
                                                  &(f->bc_coeffs->af[face_id]),
                                                  &(f->bc_coeffs->b[face_id]),
                                                  &(f->bc_coeffs->bf[face_id]),
                                                  qimp,
                                                  hint);

      }
    }

    /* Matrix
     *=======*/

    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++)
      rovsdt[cell_id] = 0.;

    /* Penalisation term for the source
     * in parallel, only one rank takes that
     * */
    if (source_c_ids[s_id] > -1)
      rovsdt[source_c_ids[s_id]] = mq->cell_vol[source_c_ids[s_id]];

    /* Even if there is no Dirichlet, the system is well-posed
     * so no need of diagonal reinforcement */
    vcopt.ndircl = 1;

    /* Right hand side and initial guess
     *==================================*/

    for (cs_lnum_t cell_id = 0; cell_id< m->n_cells_with_ghosts; cell_id++) {
      rhs[cell_id] = 0.;
      pvar[cell_id] = 0.;
    }

    /* in parallel, only one rank takes that */
    if (source_c_ids[s_id] > -1)
      rhs[source_c_ids[s_id]] = mq->cell_vol[source_c_ids[s_id]];

    /* Norm
     *======*/

    cs_real_t norm = mq->tot_vol;

    /* Solving
     *=========*/

    /* In case of a theta-scheme, set theta = 1;
       no relaxation in steady case either */

    cs_equation_iterative_solve_scalar(0,   /* idtvar: no steady state algo */
                                       -1,  /* no over loops */
                                       f->id,
                                       NULL,
                                       0,   /* iescap */
                                       0,   /* imucpp */
                                       norm,
                                       &vcopt,
                                       f->val_pre,
                                       f->val,
                                       f->bc_coeffs->a,
                                       f->bc_coeffs->b,
                                       f->bc_coeffs->af,
                                       f->bc_coeffs->bf,
                                       i_massflux,
                                       b_massflux,
                                       i_massflux, /* viscosity, not used */
                                       b_massflux, /* viscosity, not used */
                                       i_massflux, /* viscosity, not used */
                                       b_massflux, /* viscosity, not used */
                                       NULL,
                                       NULL,
                                       NULL,
                                       0, /* icvflb (upwind) */
                                       NULL,
                                       rovsdt,
                                       rhs,
                                       pvar,
                                       dpvar,
                                       NULL,
                                       NULL);

    for (cs_lnum_t cell_id = 0; cell_id< m->n_cells; cell_id++)
      f->val[cell_id] = CS_MAX(f->val[cell_id], pvar[cell_id]);

    /* Parallel synchronisation */
    cs_mesh_sync_var_scal(f->val);

  } /* End loop over sources */

  /* Finalization */

  /* Compute the correct c_w_face_normal orientation
   * using the porosity gradient with the cell data to point data
   * gradient */

  /* Enable solid cells in fluid_solid mode */
  cs_porous_model_set_has_disable_flag(0);
  // Correction of orientation c_w_face_normal to point towards the solid region
  cs_real_3_t *grdporo;
  BFT_MALLOC(grdporo, m->n_cells_with_ghosts, cs_real_3_t);

  int imrgra = eqp->imrgra;
  eqp->imrgra = 7;

  cs_field_gradient_scalar(f,
                           false, /* use_previous_t */
                           1,
                           grdporo);

  eqp->imrgra = imrgra;

  /* Disable solid cells in fluid_solid mode */
  cs_porous_model_set_has_disable_flag(1);

  cs_real_3_t *restrict c_w_face_normal
    = (cs_real_3_t *restrict)mq->c_w_face_normal;

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++) {

    if (cs_math_3_dot_product(grdporo[cell_id], c_w_face_normal[cell_id]) >
        0.0) {
      for (cs_lnum_t i = 0; i < 3; i++)
        c_w_face_normal[cell_id][i] = -c_w_face_normal[cell_id][i];
    }
  }

  // Porosity
  for (cs_lnum_t c_id = 0; c_id < m->n_cells_with_ghosts; c_id++) {
    cs_real_t porosity = f->val[c_id];
    if (porosity < cs_math_epzero) {
      f->val[c_id] = 0.;
      mq->c_disable_flag[c_id] = 1;
    }
    cell_f_vol[c_id] =  f->val[c_id] * mq->cell_vol[c_id];
  }

  /* Set interior face values */
  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

    cs_lnum_t c_id0 = i_face_cells[face_id][0];
    cs_lnum_t c_id1 = i_face_cells[face_id][1];

    cs_real_t face_porosity = CS_MIN(f->val[c_id0], f->val[c_id1]);

    for (cs_lnum_t i = 0; i < 3; i++)
      i_f_face_normal[face_id][i] = face_porosity * i_face_normal[face_id][i];

    mq->i_f_face_surf[face_id] = cs_math_3_norm(i_f_face_normal[face_id]);

    if (mq->i_f_face_factor != NULL) {
      //FIXME
      //if (face_porosity > cs_math_epzero) {
      //  mq->i_f_face_factor[face_id][0] = cpro_porosi[c_id0] / face_porosity;
      //  mq->i_f_face_factor[face_id][1] = cpro_porosi[c_id1] / face_porosity;
      //}
      //else {
        mq->i_f_face_factor[face_id][0] = 1.;
        mq->i_f_face_factor[face_id][1] = 1.;
      //}
    }

  }

  /* Set boundary face values */
  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

    cs_lnum_t c_id = b_face_cells[face_id];

    cs_real_t face_porosity = f->val[c_id];

    for (cs_lnum_t i = 0; i < 3; i++)
      b_f_face_normal[face_id][i] = face_porosity * b_face_normal[face_id][i];

    mq->b_f_face_surf[face_id] = cs_math_3_norm(b_f_face_normal[face_id]);

    if (mq->b_f_face_factor != NULL) {
      //FIXME
      //if (face_porosity > cs_math_epzero) {
      //  mq->b_f_face_factor[face_id] = cpro_porosi[c_id] / face_porosity;
      //}
      //else {
        mq->b_f_face_factor[face_id] = 1.;
      //}
    }
  }

  /* Set interior face values */
  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

    cs_lnum_t c_id0 = i_face_cells[face_id][0];
    cs_lnum_t c_id1 = i_face_cells[face_id][1];

    if (mq->c_disable_flag[c_id0] == 1 || mq->c_disable_flag[c_id1] == 1) {
      i_f_face_normal[face_id][0] = 0.;
      i_f_face_normal[face_id][1] = 0.;
      i_f_face_normal[face_id][2] = 0.;
      i_f_face_surf[face_id] = 0.;
    }
  }

  /* Set boundary face values */
  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

    cs_lnum_t c_id0 = b_face_cells[face_id];

    if (mq->c_disable_flag[c_id0] == 1) {
      b_f_face_normal[face_id][0] = 0.;
      b_f_face_normal[face_id][1] = 0.;
      b_f_face_normal[face_id][2] = 0.;
      b_f_face_surf[face_id] = 0.;
    }
  }

  /* Free memory */

  BFT_FREE(grdporo);
  BFT_FREE(_porosity_from_scan_opt.output_name);
  BFT_FREE(_porosity_from_scan_opt.file_names);
  BFT_FREE(_porosity_from_scan_opt.sources);
  BFT_FREE(_porosity_from_scan_opt.source_c_ids);

  BFT_FREE(pvar);
  BFT_FREE(dpvar);
  BFT_FREE(rhs);
  BFT_FREE(rovsdt);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
