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
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_coupling.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_geom.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_io.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_equation_iterative_solve.h"
#include "cs_physical_constants.h"
#include "cs_post.h"
#include "cs_preprocessor_data.h"
#include "cs_restart.h"
#include "cs_selector.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"

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
  .file_name = NULL,
  .output_name = NULL,
  .postprocess_points = true,
  .transformation_matrix =
  {{1., 0., 0., 0.},
    {0., 1., 0., 0.},
    {0., 0., 1., 0.}},
  .nb_sources = 0,
  .sources = NULL,
  .source_c_ids = NULL
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
 * Function
 *
 * parameters:
 *----------------------------------------------------------------------------*/

static void
_count_from_file(const cs_mesh_t *m,
                 const cs_mesh_quantities_t *mq) {

  char line[512];

  cs_real_t *restrict cell_f_vol = mq->cell_f_vol;

  /* Open file */
  bft_printf(_(" Compute the porosity from a scan points file:\n    %s\n"),
             _porosity_from_scan_opt.file_name);

  bft_printf(_("  Transformation       %12.5g %12.5g %12.5g %12.5g\n"
               "  matrix:              %12.5g %12.5g %12.5g %12.5g\n"
               "                       %12.5g %12.5g %12.5g %12.5g\n"
               "    (last column is translation vector)\n"),
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

  FILE* file = fopen(_porosity_from_scan_opt.file_name, "rt");
  if (file == NULL)
    bft_error(__FILE__,__LINE__, 0, _("Porosity from scan: Could not open file."));

  int n_points = 0;
  int n_read_points = 0;

  if (fscanf(file, "%d", &n_read_points) != 1)
    bft_error(__FILE__,__LINE__, 0, _("Porosity from scan: Could not read the number of lines."));

  bft_printf(_("Porosity from scan: %d points to be read.\n"), n_points);

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

  /* Read multiple scan file
   * ----------------------- */
  for (int n_scan = 0; n_read_points != 0; n_scan++) {
    n_points = n_read_points;
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
        bft_error(__FILE__,__LINE__, 0, _("Porosity from scan: Error while reading dataset. Line %d\n"), i);
      if (fscanf(file, "%lf", &(xyz[1])) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Porosity from scan: Error while reading dataset."));
      if (fscanf(file, "%lf", &(xyz[2])) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Porosity from scan: Error while reading dataset."));

      /* Translation and rotation */
      xyz[3] = 1.;
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 4; k++)
          point_coords[i][j] += _porosity_from_scan_opt.transformation_matrix[j][k] * xyz[k];

        /* Compute bounding box*/
        min_vec[j] = CS_MIN(min_vec[j], point_coords[i][j]);
        max_vec[j] = CS_MAX(max_vec[j], point_coords[i][j]);
      }

      /* Intensities */
      if (fscanf(file, "%d", &num) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Porosity from scan: Error while reading dataset."));

      /* Red */
      if (fscanf(file, "%d", &red) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Porosity from scan: Error while reading dataset."));
      /* Green */
      if (fscanf(file, "%d", &green) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Porosity from scan: Error while reading dataset."));
      /* Blue */
      if (fscanf(file, "%d\n", &blue) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Porosity from scan: Error while reading dataset."));

      /* When colors are written as int, Paraview intreprates them in [0, 255]
       * when they are written as float, Paraview interprates them in [0., 1.]
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

    if (n_read_points > 0)
      bft_printf(_("  Porosity from scan: %d additional points to be read.\n\n"),
                 n_read_points);

    /* FVM meshes for writers */
    if (_porosity_from_scan_opt.postprocess_points) {
      char *fvm_name;
      if (_porosity_from_scan_opt.output_name == NULL) {
        BFT_MALLOC(fvm_name,
                   strlen(_porosity_from_scan_opt.file_name) + 3 + 1,
                   char);
        strcpy(fvm_name, _porosity_from_scan_opt.file_name);
      } else {
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
      fvm_writer_t *writer = fvm_writer_init(fvm_name,
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

    ple_locator_set_mesh(_locator,
                         location_mesh,
                         options,
                         0., /* tolerance_base */
                         0.1, /* tolerance */
                         3, /* dim */
                         n_points,
                         NULL,
                         NULL, /* point_tag */
                         (cs_real_t *)point_coords,
                         NULL, /* distance */
                         cs_coupling_mesh_extents,
                         cs_coupling_point_in_mesh_p);

    /* Shift from 1-base to 0-based locations */
    ple_locator_shift_locations(_locator, -1);

    /* dump locator */
#if 0
    ple_locator_dump(_locator);
#endif

    /* Get the element ids (list of points on the local rank) */
    cs_lnum_t n_points_loc = ple_locator_get_n_dist_points(_locator);

#if 0
    bft_printf("ple_locator_get_n_dist_points = %d, n_points = %d\n",
               n_points_loc, n_points);
#endif

    const cs_lnum_t *elt_ids = ple_locator_get_dist_locations(_locator);

    for (int i = 0; i < n_points_loc; i++) {
      if (elt_ids[i] >= 0) { /* Found */
        /* Could be improved with a prallel reading */
        f_nb_scan->val[elt_ids[i]] += 1./cs_glob_n_ranks;
      }
    }

    /* Free memory */
    _locator = ple_locator_destroy(_locator);
    BFT_FREE(point_coords);
    BFT_FREE(colors);

  } /* End loop on multiple scans */

  if (fclose(file) != 0)
    bft_error(__FILE__,__LINE__, 0, _("Porosity from scan: Could not close the file."));

  /* Nodal mesh is not needed anymore */
  location_mesh = fvm_nodal_destroy(location_mesh);

  int nb_points_thresholds = 10;
  /* Synchronize nb_scan_points */
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    cell_f_vol[cell_id] = mq->cell_vol[cell_id];
    if (f_nb_scan->val[cell_id] > nb_points_thresholds) {
      cell_f_vol[cell_id] = 0.;
      mq->c_disable_flag[cell_id] = 1;
    }
  }

  cs_real_3_t *restrict i_face_normal =
     (cs_real_3_t *restrict)mq->i_face_normal;
  cs_real_3_t *restrict b_face_normal =
     (cs_real_3_t *restrict)mq->b_face_normal;
  cs_real_t *restrict i_face_surf =
     (cs_real_t *restrict)mq->i_face_surf;
  cs_real_t *restrict b_face_surf =
     (cs_real_t *restrict)mq->b_face_surf;

  cs_real_3_t *restrict i_f_face_normal =
     (cs_real_3_t *restrict)mq->i_f_face_normal;
  cs_real_3_t *restrict b_f_face_normal =
     (cs_real_3_t *restrict)mq->b_f_face_normal;
  cs_real_t *restrict i_f_face_surf =
     (cs_real_t *restrict)mq->i_f_face_surf;
  cs_real_t *restrict b_f_face_surf =
     (cs_real_t *restrict)mq->b_f_face_surf;

  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
    cs_lnum_t cell_id1 = m->i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = m->i_face_cells[face_id][1];
    if (cell_f_vol[cell_id1] <= 0. || cell_f_vol[cell_id2] <= 0.) {
      i_f_face_normal[face_id][0] = 0.;
      i_f_face_normal[face_id][1] = 0.;
      i_f_face_normal[face_id][2] = 0.;
      i_f_face_surf[face_id] = 0.;
    } else {
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
    } else {
      b_f_face_normal[face_id][0] = b_face_normal[face_id][0];
      b_f_face_normal[face_id][1] = b_face_normal[face_id][1];
      b_f_face_normal[face_id][2] = b_face_normal[face_id][2];
      b_f_face_surf[face_id]      = b_face_surf[face_id]     ;
    }
  }

  return;
}

/*----------------------------------------------------------------------------
 * Get pointer
 *----------------------------------------------------------------------------*/

void
cs_f_porosity_from_scan_get_pointer(bool **compute_porosity_from_scan)
{
  *compute_porosity_from_scan =
    &(_porosity_from_scan_opt.compute_porosity_from_scan);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of points for the computation of the
 * porosity from scan.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_set_file_name(const char *file_name)
{
  if (file_name == NULL) {
    _porosity_from_scan_opt.compute_porosity_from_scan = false;
    return;
  }

  _porosity_from_scan_opt.compute_porosity_from_scan = true;

  BFT_MALLOC(_porosity_from_scan_opt.file_name,
             strlen(file_name) + 1,
             char);

  sprintf(_porosity_from_scan_opt.file_name, "%s", file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function sets the output name for the FVM writer of scan points.
 *
 * \param[in] output_name  name of the output (a suffix will be added)
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_set_output_name(const char *output_name)
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
 * \brief This function add a scanner source point
 *
 * \param[in] source     source vector
 * \param[in] transform  flag to apply the transformation matrix to the source
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_add_source(const cs_real_3_t source,
                                 const bool transform)
{
  /* Add a source */
  const int s_id = _porosity_from_scan_opt.nb_sources;
  _porosity_from_scan_opt.nb_sources++;

  BFT_REALLOC(
      _porosity_from_scan_opt.source_c_ids,
      _porosity_from_scan_opt.nb_sources,
      cs_lnum_t);

  BFT_REALLOC(
      _porosity_from_scan_opt.sources,
      _porosity_from_scan_opt.nb_sources,
      cs_real_3_t);

  if (transform) {
    /* Apply translation and rotation */
    for (int i = 0; i < 3; i++) {
      _porosity_from_scan_opt.sources[s_id][i] = 0;
      for (int j = 0; j < 3; j++)
        _porosity_from_scan_opt.sources[s_id][i] +=
          _porosity_from_scan_opt.transformation_matrix[i][j] * source[j];
      _porosity_from_scan_opt.sources[s_id][i] +=
        _porosity_from_scan_opt.transformation_matrix[i][3];
    }
  } else {
    for (int i = 0; i < 3; i++)
      _porosity_from_scan_opt.sources[s_id][i] = source[i];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes the porosity which is equal to one from
 *        a source, radiating sphericaly, and is 0 when touching points
 *        of the scan.
 *
 *  This function solves the following transport equation on \f$ \varia \f$:
 *  \f[
 *  \dfrac{\partial \varia}{\partial t} + \divs \left( \varia \vect{e}_r \right)
 *      - \divs \left( \vect{e}_r \right) \varia = 0
 *  \f]
 *  where \f$ \vect{e}_r = \dfrac{\vect{x} - \vect{x}_0}{\norm{\vect{x} - \vect{x}_0}} \f$
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
  /* Initialization
   *===============*/

  const cs_domain_t *domain = cs_glob_domain;
  const cs_mesh_t *m = domain->mesh;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;

  const cs_real_3_t *restrict cell_cen =
     (const cs_real_3_t *restrict)mq->cell_cen;
  const cs_real_3_t *restrict i_face_normal =
     (const cs_real_3_t *restrict)mq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal =
     (const cs_real_3_t *restrict)mq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog =
     (const cs_real_3_t *restrict)mq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog =
     (const cs_real_3_t *restrict)mq->b_face_cog;

  /* Pointer to porosity field */
  cs_field_t *f = cs_field_by_name_try("porosity_w_field");

  cs_lnum_t *source_c_ids = _porosity_from_scan_opt.source_c_ids;
  int nb_sources = _porosity_from_scan_opt.nb_sources;

  /* First pass to put fluid surfaces to 0 for all faces' cells with
   * at least 3 (???) points */
  _count_from_file(m, mq);

  cs_real_t *restrict i_massflux =
    cs_field_by_id(
        cs_field_get_key_int(f, cs_field_key_id("inner_mass_flux_id")))->val;
  cs_real_t *restrict b_massflux =
    cs_field_by_id(
        cs_field_get_key_int(f, cs_field_key_id("boundary_mass_flux_id")))->val;

  cs_var_cal_opt_t vcopt;
  cs_field_get_key_struct(f, cs_field_key_id("var_cal_opt"), &vcopt);

  /* Local variables */
  cs_real_t *rovsdt, *pvar, *dpvar;
  cs_real_t *rhs;
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
      cs_real_3_t normal;
      /* Normal direction is given by x-x0 */
      cs_math_3_normalise(x0xf, normal);

      i_massflux[face_id] = cs_math_3_dot_product(normal, i_face_normal[face_id]);
    }

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
      cs_real_t x0xf[3] = {
        b_face_cog[face_id][0] - source_cen[0],
        b_face_cog[face_id][1] - source_cen[1],
        b_face_cog[face_id][2] - source_cen[2]
      };
      cs_real_3_t normal;
      /* Normal direction is given by x-x0 */
      cs_math_3_normalise(x0xf, normal);

      b_massflux[face_id] = cs_math_3_dot_product(normal, b_face_normal[face_id]);
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
      } else {

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

    /* Penalisation terme for the source
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
                                       f->name,
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

  } /* End loop over sources */

  /* Free memory */
  BFT_FREE(pvar);
  BFT_FREE(dpvar);
  BFT_FREE(rhs);
  BFT_FREE(rovsdt);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
