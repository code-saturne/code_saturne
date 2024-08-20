/*============================================================================
 * This file define functions needed for Calculate mean profile
 * (scalars or vectors) over the mesh.
 *
 * Example in cs_user_extra_operations_mean_profiles.c
 *============================================================================*/

/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/* User header */

#include "cs_user_profile.h"

/*----------------------------------------------------------------------------
 * MEDCOUPLING library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MEDCOUPLING)

#include "cs_medcoupling_intersector.h"
#include "cs_medcoupling_mesh.hxx"

#include <MEDCoupling_version.h>

#include <MEDCouplingCMesh.hxx>
#include <MEDCouplingUMesh.hxx>

#include <MEDCouplingField.hxx>
#include <MEDCouplingFieldDouble.hxx>
#include <MEDCouplingFieldFloat.hxx>

#include <MEDCouplingRemapper.hxx>

#if defined(HAVE_MEDCOUPLING_LOADER)

#include <MEDFileMesh.hxx>
#include <MEDFileField1TS.hxx>
#include <MEDFileFieldMultiTS.hxx>

#include <MEDLoader.hxx>

#endif

#include "Interpolation3D.hxx"
#include <MEDCouplingNormalizedUnstructuredMesh.txx>

using namespace MEDCoupling;
#endif

/*----------------------------------------------------------------------------
 * OpenTurns library headers
 *----------------------------------------------------------------------------*/

#if HAVE_OT == 1
#include "openturns/OT.hxx"
#endif

/*
 * Private functions and global variables
 */

static user_profile_info_t _profile_list = { NULL, 0 };

/*----------------------------------------------------------------------------
 *  Intersector structure for MedCoupling
 *----------------------------------------------------------------------------*/

struct _user_profile_med_t {

#if defined(HAVE_MEDCOUPLING)
  MEDCouplingUMesh     **layer_mesh;
  cs_medcoupling_mesh_t *local_mesh;
#else
  void **layer_mesh;
  void  *local_mesh;
#endif
};

/*----------------------------------------------------------------------------
 * Compute the minimum and the maximum of a vector (locally).
 *
 * parameters:
 *   n_vals    <-- local number of elements
 *   var       <-- pointer to vector
 *   min       --> minimum
 *   max       --> maximum
 *----------------------------------------------------------------------------*/

static void
_compute_min_max(cs_lnum_t       n_vals,
                 const cs_real_t var[],
                 cs_real_t      *min,
                 cs_real_t      *max)
{
  cs_real_t _min = DBL_MAX, _max = -DBL_MAX;

  for (cs_lnum_t i = 0; i < n_vals; i++) {
    _min = CS_MIN(_min, var[i]);
    _max = CS_MAX(_max, var[i]);
  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    MPI_Allreduce(&_min, min, 1, CS_MPI_REAL, MPI_MIN, cs_glob_mpi_comm);

    MPI_Allreduce(&_max, max, 1, CS_MPI_REAL, MPI_MAX, cs_glob_mpi_comm);
  }

#endif

  if (cs_glob_n_ranks == 1) {
    *min = _min;
    *max = _max;
  }
}

/*----------------------------------------------------------------------------
 * Allocate memory or mesh med structure.
 * Create a med mesh structure, allocate memory and return pointer to struct
 *
 * parameters:
 *   n_layers   <-- number of layer
 *----------------------------------------------------------------------------*/

static user_profile_med_t *
_allocate_med_mesh_struct([[maybe_unused]] cs_lnum_t n_layers)
{
  user_profile_med_t *med_t = NULL;
  BFT_MALLOC(med_t, 1, user_profile_med_t);
#if defined(HAVE_MEDCOUPLING)
  BFT_MALLOC(med_t->layer_mesh, n_layers, MEDCouplingUMesh *);
  for (int layer_id = 0; layer_id < n_layers; layer_id++) {
    char name[10];
    snprintf(name, 10, "%s_%d", "layer", layer_id);
    MEDCoupling::MEDCouplingUMesh *mesh
      = MEDCoupling::MEDCouplingUMesh::New(name, 3);
    med_t->layer_mesh[layer_id] = mesh;
  }

  med_t->local_mesh = cs_medcoupling_mesh_from_base(cs_glob_mesh,
                                                    "Case_global_mesh",
                                                    "all[]",
                                                    3,
                                                    1);

#endif

  return med_t;
}

/*----------------------------------------------------------------------------
 * Create a 1d sample (same formalism than OpenTurns) with weigths for a layer
 *
 * parameters:
 *   profile        <-- pointer to the current profile structure
 *   n_elts_sample  --> number of elements in the sample
 *   sample         --> preallocated array, populated with selected cells
 *   weights        --> preallocated array, populated with cells weigth
 *   layer_id       <-- id of layer to be selected
 *----------------------------------------------------------------------------*/

static void
_create_1d_sample_(user_profile_t  *profile,
                   cs_lnum_t       *n_elts_sample,
                   cs_real_t       *sample,
                   cs_real_t       *weights,
                   int              layer_id)
{
  /* Profile shorter variables:

     Also available but not used here:
       cs_lnum_t   n_layers = profile->n_layers;
  */
  const char *field    = profile->field;
  const char *weighted = profile->weighted;

  /* Get mesh quantities */
  const cs_lnum_t  n_cells  = cs_glob_mesh->n_cells;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  /* Get physical fields */
  const cs_field_t *f       = cs_field_by_name_try(field);
  const cs_field_t *density = cs_field_by_name_try("density");
  double            ro0     = cs_glob_fluid_properties->ro0;

  /* Reset n_elts_sample */
  *n_elts_sample = 0;

  if (f == NULL) {
    bft_printf("Warning field %s set for profile calc does not exist\n",
               field);
    return;
  }
  cs_lnum_t f_dim = f->dim; // dimension of fields

  // Define pointer and variable for cs_selector
  cs_lnum_t  n_selected_cells = 0;
  cs_lnum_t *selected_cells   = NULL;

  // Allocate memory for the cells list which will be populated by cs_selector
  BFT_MALLOC(selected_cells, n_cells, cs_lnum_t);

  cs_selector_get_cell_list(profile->criteria,
                            &n_selected_cells,
                            selected_cells);

  cs_real_t sel_cells_weight = 0.0;

  for (cs_lnum_t ii = 0; ii < n_selected_cells; ii++) {
    cs_lnum_t c_id = selected_cells[ii];

    cs_real_t cell_percent_layer = profile->cells_layer_vol[layer_id][c_id];
    cs_real_t weight_cell        = 0.0;

    if (strcmp(weighted, "MASS") == 0) {
      if (density != NULL)
        weight_cell = cell_vol[c_id] * density->val[c_id];
      else
        weight_cell = cell_vol[c_id] * ro0;
    }
    else if (strcmp(weighted, "VOLUME") == 0)
      weight_cell = cell_vol[c_id];
    else if (strcmp(weighted, "NO") == 0) {
      weight_cell = 1.0; /*Unweighted profile*/
    }

    sel_cells_weight += weight_cell;
    weight_cell = weight_cell * cell_percent_layer;

    if (weight_cell > 0.0) {
      cs_real_t f_val = 0.0;
      if (f_dim == 1)
        f_val = f->val[c_id];
      else {
        cs_real_t f_val_norm = 0.0;
        for (int j = 0; j < f_dim; j++)
          f_val_norm += pow(f->val[f_dim * c_id + j], 2.0);

        f_val_norm = pow(f_val_norm, 0.5);
        f_val      = f_val_norm;
      }

      /*Fill Sample and weight if the cell belong to the layer*/
      sample[*n_elts_sample]  = f_val;
      weights[*n_elts_sample] = weight_cell;
      *n_elts_sample += 1;
    }
  }

  /*Total selected cells weigth could be calculated once and not for each layer
   * but keep as this because of current code structure (avoid de create a
   * dedicated function and re compute each cell weight*/
  cs_parall_sum(1, CS_REAL_TYPE, &sel_cells_weight);
  profile->sel_cells_weigth = sel_cells_weight;

  BFT_FREE(selected_cells);
}

/*----------------------------------------------------------------------------
 * Compute mean and standard of a given weighted 1d sample
 *
 * parameters:
 *   sample           <-- sample on which compute moment
 *   weights          <-> weights associated to sample, return normalized
 *   n_elts_sample    <-- number of elements in the sample
 *   moment           --> mu, sigma of the weighted sample
 *   min_max          --> min and max of the sample
 *----------------------------------------------------------------------------*/

static void
_compute_sample_moment(cs_real_t *sample,
                       cs_real_t *weights,
                       cs_real_t  n_elts_sample,
                       cs_real_t *moment,
                       cs_real_t *min_max)
{
  /* Reset moment
     It is assumed moment[0] = mu moment[1] = sigma */

  moment[0]          = 0.0;
  moment[1]          = 0.0;
  cs_real_t mu       = 0.0;
  cs_real_t sigma    = 0.0;
  cs_real_t variance = 0.0;

  /* Normalize weight */
  cs_real_t *w_n;
  cs_real_t  w_tot = 0.0;

  BFT_MALLOC(w_n, n_elts_sample, cs_real_t);

  for (cs_lnum_t iel = 0; iel < n_elts_sample; iel++)
    w_tot += weights[iel];

  /* Sum over all MPI ranks */
  cs_parall_sum(1, CS_REAL_TYPE, &w_tot);

  /* Normalized weights and update provided weights */
  for (cs_lnum_t iel = 0; iel < n_elts_sample; iel++) {
    w_n[iel]     = weights[iel] / w_tot;
    weights[iel] = w_n[iel];
  }

  /* Compute sample mean */
  for (cs_lnum_t iel = 0; iel < n_elts_sample; iel++)
    mu += w_n[iel] * sample[iel];

  /* Sum over all MPI ranks */
  cs_parall_sum(1, CS_REAL_TYPE, &mu);

  /* Compute sample variance */
  for (cs_lnum_t iel = 0; iel < n_elts_sample; iel++)
    variance += w_n[iel] * pow(sample[iel] - mu, 2.0);

  /* Sum over all MPI ranks */
  cs_parall_sum(1, CS_REAL_TYPE, &variance);

  sigma = pow(variance, 1.0 / 2.0);

  /* Update sample moment */
  moment[0] = mu;
  moment[1] = sigma;

  /* Get sample min and max */
  cs_real_t min_sample = 0.0;
  cs_real_t max_sample = 0.0;

  _compute_min_max(n_elts_sample, sample, &min_sample, &max_sample);

  min_max[0] = min_sample;
  min_max[1] = max_sample;

  BFT_FREE(w_n);
}

/*----------------------------------------------------------------------------
 * Fill histogram classes (width and height) for a provided bandwidth.
 * sample weights are assumed to be normalized.
 * bandwidth is updated given the final int number of bins.
 *
 * parameters:
 *   histogram        <-> pointer to histogram structure
 *   sample           <-- 1d sample on which the histogram is built
 *   weights          <-- weights associated to sample
 *   n_elts_sample    <-- number of elements in the sample
 *   bandwidth        <-> targeted bandwidth
 *----------------------------------------------------------------------------*/

static void
_fill_histogram_classes_u_bandwidth(user_histogram_t  *histogram,
                                    cs_real_t         *sample,
                                    cs_real_t         *weights,
                                    cs_lnum_t          n_elts_sample,
                                    cs_real_t         *bandwidth)
{
  /* Histogram shorter variables

     Also available but not used here:
       cs_real_t mu        = histogram->mean;
       cs_real_t sigma     = histogram->sd;
  */
  cs_lnum_t n_bins_max = histogram->n_bins_max;
  cs_real_t min        = histogram->min;
  cs_real_t max        = histogram->max;

  max              = (max - min) * 0.001 + max;
  histogram->max = max; /* ensure all values are caught */

  cs_lnum_t n_bins = n_bins_max;

  if (*bandwidth > DBL_EPSILON) {
    n_bins = (cs_lnum_t)((max - min) / (*bandwidth));
  }
  else
    n_bins = n_bins_max;

  if (n_bins > n_bins_max)
    n_bins = n_bins_max;

  *bandwidth = (max - min) / ((cs_real_t)n_bins);

  histogram->l_i[0] = *bandwidth;
  histogram->c_i[0] = min + histogram->l_i[0] / 2.0;

  /* Populate histogram classes */
  for (int b_id = 1; b_id < n_bins; b_id++) {
    histogram->l_i[b_id] = *bandwidth;
    histogram->c_i[b_id] =   histogram->c_i[b_id - 1]
                           + histogram->l_i[b_id - 1] / 2.0
                           + histogram->l_i[b_id] / 2.0;
  }

  /* Populate histogram class with weights */
  for (int b_id = 0; b_id < n_bins; b_id++)
    histogram->h_i[b_id] = 0.0;

  /* Compute height of each class on each rank */

  for (cs_lnum_t iel = 0; iel < n_elts_sample; iel++) {
    for (int b_id = 0; b_id < n_bins; b_id++) {
      cs_real_t b_min = histogram->c_i[b_id] - histogram->l_i[b_id] / 2.0;
      cs_real_t b_max = histogram->c_i[b_id] + histogram->l_i[b_id] / 2.0;
      if (sample[iel] >= b_min && sample[iel] < b_max) {
        histogram->h_i[b_id] += weights[iel];
        break;
      }
    }
  }

  /*Sum height classes over all MPI ranks*/

  cs_parall_sum(n_bins, CS_REAL_TYPE, histogram->h_i);

  /*update histogram*/
  histogram->n_bins = n_bins;
}

/*----------------------------------------------------------------------------
 * Compute histogram quantiles given a provided quantile array
 *
 * parameters:
 *   histogram <-> pointer to histogram structure
 *   quantile    <-> pointer to quantile array {X,proba}
 *   n_quantile  <-- number of quantiles to be computed
 *----------------------------------------------------------------------------*/

static void
_compute_histogram_quantile(user_histogram_t  *histogram,
                            cs_real_t          quantile[][2],
                            cs_lnum_t          n_quantile)
{
  /* Profile shorter variables */
  cs_real_t *h_i = histogram->h_i;
  cs_real_t *c_i = histogram->c_i;
  cs_real_t *l_i = histogram->l_i;

  cs_lnum_t n_bins = histogram->n_bins;
  cs_real_t prob_accum = 0.0;
  cs_lnum_t q_id_start = 0;

  for (int b_id = 0; b_id < n_bins; b_id++) {
    prob_accum += histogram->h_i[b_id];

    for (int q_id = q_id_start; q_id < n_quantile; q_id++) {
      if (prob_accum > quantile[q_id][1]) {
        q_id_start++;
        cs_real_t dp_add      = quantile[q_id][1] - (prob_accum - h_i[b_id]);
        cs_real_t percent_l_i = dp_add / h_i[b_id];
        quantile[q_id][0]
          = c_i[b_id] - l_i[b_id] / 2.0 + l_i[b_id] * percent_l_i;

      } /* End of quantile detection */
    }  /* End of loop on quantiles */

  } /* End of loop on classes */

  /* Update histogram proba tot:
     it is possible to check the correct sum of h_i */
  histogram->ptot = prob_accum;
}

/*----------------------------------------------------------------------------
 * Fill histogram classes (width and height). Optimal bandwith is estimated
 * based on gathered histogram from the sample dispatched on the n MPI ranks.
 * A first bandwidth assuming gaussian distribution is used
 * (Scott rule: see OpenTurns Histogram factory class).
 * Then the bandwitdh is iteratively precised using: AMISE-optimal bandwith,
 * known as Freedman and Diaconis rule [freedman1981]: see OpenTurns
 *
 * parameters:
 *   histogram     <-> pointer to histogram structure
 *   sample        <-- 1d sample on which the histogram has to be build
 *   weights       <-- weights associated to sample
 *   n_elts_sample <-- number of elements in the sample
 *----------------------------------------------------------------------------*/

static void
_compute_histogram(user_histogram_t  *histogram,
                   cs_real_t         *sample,
                   cs_real_t         *weights,
                   cs_lnum_t          n_elts_sample)
{
  /* Histogram shorter variable */
  cs_real_t sigma      = histogram->sd;
  cs_lnum_t n_bins_max = histogram->n_bins_max;

  cs_lnum_t n_gelts_sample = n_elts_sample;

  /* Get the total number of elements in sample over mpi ranks */
  cs_parall_sum(1, CS_INT_TYPE, &n_gelts_sample);

  /* Compute optimal bandwidth assuming gaussian distribution of sample
     - Scott rule */
  cs_real_t bandwidth = 0.0;

  bandwidth = sigma * pow(24.0 * pow(3.14, 0.5) / n_gelts_sample, 1.0 / 3.0);

  _fill_histogram_classes_u_bandwidth(histogram,
                                      sample,
                                      weights,
                                      n_elts_sample,
                                      &bandwidth);

  /*compute quartile based on gathered histogram and update bandwidth*/

  cs_real_t quantile[][2] = {{0.0, 0.25}, {0.0, 0.5}, {0.0, 0.75}};

  cs_lnum_t n_quantile = sizeof(quantile) / (2 * sizeof(cs_real_t *));

  _compute_histogram_quantile(histogram, quantile, n_quantile);

  /*Update histogram quantile*/

  histogram->Q1 = quantile[0][0];
  histogram->Q2 = quantile[1][0];
  histogram->Q3 = quantile[2][0];

  cs_real_t bandwidth_update = 0.0;
  cs_lnum_t n_bins           = n_bins_max;
  cs_real_t min              = histogram->min;
  cs_real_t max              = histogram->max;

  cs_real_t IQR              = histogram->Q3 - histogram->Q1;
  cs_real_t IQR_pre          = IQR;
  cs_real_t IQR_var          = 1.0;

  /*AMISE-optimal bandwith, known as Freedman and Diaconis rule [freedman1981]:
   * see OpenTurns*/
  bandwidth_update
    = IQR / (2 * 0.75) * pow(24.0 * pow(3.14, 0.5) / n_gelts_sample, 1.0 / 3.0);

  if (bandwidth_update > DBL_EPSILON) {
    n_bins = (cs_lnum_t)((max - min) / (bandwidth_update));
  }
  else
    n_bins = n_bins_max;

  if (n_bins > n_bins_max)
    n_bins = n_bins_max;

  bandwidth_update = (max - min) / ((cs_real_t)n_bins);

  cs_lnum_t n_loop = 1;
  cs_real_t bandwidth_variation
    = CS_ABS(bandwidth - bandwidth_update) / bandwidth;

  /*perform a loop to try top optimize bandwidth*/
  while (n_loop < 10 && bandwidth_variation > 0.05 && IQR > DBL_EPSILON
         && IQR_var > 0.02 && histogram->n_bins < histogram->n_bins_max) {

    bandwidth = bandwidth_update;

    _fill_histogram_classes_u_bandwidth(
      histogram, sample, weights, n_elts_sample, &bandwidth);

    _compute_histogram_quantile(histogram, quantile, n_quantile);

    /*Update histogram quantile*/

    histogram->Q1 = quantile[0][0];
    histogram->Q2 = quantile[1][0];
    histogram->Q3 = quantile[2][0];

    IQR_pre = IQR;
    IQR     = histogram->Q3 - histogram->Q1;

    bandwidth_update = IQR / (2 * 0.75)
                       * pow(24.0 * pow(3.14, 0.5) / n_gelts_sample, 1.0 / 3.0);

    n_bins = n_bins_max;
    min    = histogram->min;
    max    = histogram->max;

    if (bandwidth_update > DBL_EPSILON) {
      n_bins = (cs_lnum_t)((max - min) / (bandwidth_update));
    }
    else
      n_bins = n_bins_max;

    if (n_bins > n_bins_max)
      n_bins = n_bins_max;

    bandwidth_update = (max - min) / ((cs_real_t)n_bins);

    bandwidth_variation = CS_ABS(bandwidth - bandwidth_update) / bandwidth;
    IQR_var             = CS_ABS(IQR - IQR_pre) / IQR;

    n_loop++;
  }

  histogram->n_iter = n_loop;
}

/*----------------------------------------------------------------------------
 * Output histogram graphs using OpenTunrs library (need to be link to CS)
 *
 * parameters:
 *   histogram <-- pointer to the current histogram structure
 *   dirname   <-- pointer to directory output path
 *----------------------------------------------------------------------------*/

static void
_output_histogram_ot([[maybe_unused]]user_histogram_t  *histogram,
                     [[maybe_unused]]const char        *dirname)
{
#if HAVE_OT == 1

  /* Histogram shorter variable */
  cs_lnum_t n_bins = histogram->n_bins;
  cs_real_t min    = histogram->min;
  cs_real_t max    = histogram->max;

  cs_real_t min_li = DBL_MAX;

  OT::Point h_i(n_bins), l_i(n_bins);

  for (int b_id = 0; b_id < n_bins; b_id++) {
    h_i[b_id] = histogram->h_i[b_id];
    l_i[b_id] = histogram->l_i[b_id];
    min_li    = CS_MIN(l_i[b_id], min_li);
  }

  if (min_li < DBL_EPSILON * 1e4)
    return;

  /* Ensure surface histogram is higher than zero */
  OT::Histogram histogram_ot(min, l_i, h_i);

  /* Set description and name */
  char title_h[200];

  sprintf(title_h,
          "%s \n mu: %.2f - sigma: %.2f",
          histogram->name,
          histogram->mean,
          histogram->sd);

  histogram_ot.setDescription(OT::Description(1, histogram->field));

  /* Get a sample from the histogram */
  OT::UnsignedInteger s_size = 500;

  OT::Sample sample = histogram_ot.getSample(s_size);

  OT::KernelSmoothing ks;

  /* Fit the sample with KernelSmoothing distribution */
  OT::Distribution fitteddist = ks.build(sample);

  OT::Graph graph, graph_ks;

  graph = histogram_ot.drawPDF();
  graph.setLegendPosition("");
  graph_ks = fitteddist.drawPDF();

  graph_ks.setColors(OT::Description(1, "blue"));
  graph_ks.setLegendPosition("");

  graph.add(graph_ks);

  graph.setTitle(title_h);

  char filename[300];
  sprintf(filename, "%s/%s.png", dirname, histogram->name);

  graph.draw(filename);

#endif
}

/*----------------------------------------------------------------------------
 * Calculate min max sel mesh quantities in a orthonormal base vith dir_v
 * normalized as k vector
 *
 * The folowing members of the histogram_t structure are updated
 *   i_v: set of an abritray i_v vector otthorgonal to dir
 *   min_i: minimum mesh sel in i direction
 *   max_i: maximum mesh sel in i direction
 *   j_v: thirf vector to form orthormal base
 *   min_j: minimum mesh sel in j direction
 *   max_j: maximum mesh sel in j direction
 *   min_dir: minimum mesh sel in profile dir vector
 *   max_dir: maximum mesh sel in profile dir vector
 *
 * parameters:
 *   profile <-> pointer to the current profile structure
 *----------------------------------------------------------------------------*/

static void
_calculate_min_max_dir(user_profile_t  *profile)
{
  // Get mesh quantities
  const cs_mesh_t            *m        = cs_glob_mesh;
  const cs_lnum_t  n_vertices          = cs_glob_mesh->n_vertices;
  const cs_real_t *vtx_coord           = cs_glob_mesh->vtx_coord;
  // const cs_mesh_quantities_t *mq    = cs_glob_mesh_quantities;

  // Define pointer and variable for cs_selector
  cs_lnum_t  n_selected_cells = 0;
  cs_lnum_t *selected_cells   = NULL;

  // Allocate memory for the cells list which will be populated by cs_selector
  BFT_MALLOC(selected_cells, m->n_cells, cs_lnum_t);

  cs_selector_get_cell_list(profile->criteria,
                            &n_selected_cells,
                            selected_cells);

  cs_lnum_t  n_selected_vertices = 0;
  cs_lnum_t *selected_vertices   = NULL;

  BFT_MALLOC(selected_vertices, n_vertices, cs_lnum_t);
  cs_selector_get_cell_vertices_list_by_ids(n_selected_cells,
                                            selected_cells,
                                            &n_selected_vertices,
                                            selected_vertices);

  // Vector profile quantities
  cs_real_t dir_norm = cs_math_3_norm(profile->dir_v);

  cs_real_t n_v[3];
  n_v[0] = profile->dir_v[0] / dir_norm;
  n_v[1] = profile->dir_v[1] / dir_norm;
  n_v[2] = profile->dir_v[2] / dir_norm;

  // Set to arbitrary orthoganal vector within plane
  cs_real_t i_v[3];
  cs_real_t j_v[3];
  if (fabs(n_v[0]) > 1.e-20) {
    i_v[1] = 1.0;
    i_v[2] = 0.0;
    i_v[0] = (-i_v[1] * n_v[1] - i_v[2] * n_v[2]) / n_v[0];
  }
  else if (fabs(n_v[1]) > 1.e-20) {
    i_v[0] = 1.0;
    i_v[2] = 0.0;
    i_v[1] = (-i_v[0] * n_v[0] - i_v[2] * n_v[2]) / n_v[1];
  }
  else if (fabs(n_v[2]) > 1.e-20) {
    i_v[0] = 1.0;
    i_v[1] = 0.0;
    i_v[2] = (-i_v[0] * n_v[0] - i_v[1] * n_v[1]) / n_v[2];
  }

  cs_real_t i_v_norm;
  i_v_norm = pow(pow(i_v[0], 2.0) + pow(i_v[1], 2.0) + pow(i_v[2], 2.0), 0.5);

  // Normalized i_v_n
  i_v[0]   = i_v[0] / i_v_norm;
  i_v[1]   = i_v[1] / i_v_norm;
  i_v[2]   = i_v[2] / i_v_norm;
  i_v_norm = 1.0;

  // j_v is calculated by j_v = n_v x i_v (both normalized vector)
  j_v[0] = n_v[1] * i_v[2] - n_v[2] * i_v[1];
  j_v[1] = n_v[2] * i_v[0] - n_v[0] * i_v[2];
  j_v[2] = n_v[0] * i_v[1] - n_v[1] * i_v[0];

  cs_real_t *base_orth[] = { i_v, j_v, n_v };

  cs_real_t min_mesh_ijn[3] = { 0.0, 0.0, 0.0 };
  cs_real_t max_mesh_ijn[3] = { 0.0, 0.0, 0.0 };

  cs_real_t **vtx_dist = NULL;

  BFT_MALLOC(vtx_dist, 3, cs_real_t *);
  for (int k = 0; k < 3; k++)
    BFT_MALLOC(vtx_dist[k], n_selected_vertices, cs_real_t);

  for (cs_lnum_t ii = 0; ii < n_selected_vertices; ii++) {
    cs_lnum_t vtx_id = selected_vertices[ii];
    for (int v_id = 0; v_id < 3;
         v_id++) { /* Loop on the 3 vector of the new base */
      vtx_dist[v_id][ii] = 0.0;
      for (int jj = 0; jj < 3; jj++) /* Loop on x y z initial direction */
        vtx_dist[v_id][ii]
          += vtx_coord[3 * vtx_id + jj]
             * base_orth[v_id][jj]; /* Project each point in the new base */
    }
  }

  for (int k = 0; k < 3; k++) { /* Find min max in the new base*/
    _compute_min_max(
      n_selected_vertices, vtx_dist[k], &min_mesh_ijn[k], &max_mesh_ijn[k]);
  }

  /* Update i_v and j_v vector in structure */
  for (int k = 0; k < 3; k++)
    profile->i_v[k] = i_v[k];

  for (int k = 0; k < 3; k++)
    profile->j_v[k] = j_v[k];

  /* Update min max dimension of the selected mesh part in the structure */
  profile->min_i   = min_mesh_ijn[0];
  profile->max_i   = max_mesh_ijn[0];
  profile->min_j   = min_mesh_ijn[1];
  profile->max_j   = max_mesh_ijn[1];
  profile->min_dir = min_mesh_ijn[2];
  profile->max_dir = max_mesh_ijn[2];

  BFT_FREE(selected_cells);
  BFT_FREE(selected_vertices);
  for (int k = 0; k < 3; k++)
    BFT_FREE(vtx_dist[k]);
  BFT_FREE(vtx_dist);
}

/*----------------------------------------------------------------------------
 * Comput thickness of each layer profile given its law
 *
 * parameters:
 *   profile  <-> pointer to the current profile structure
 *----------------------------------------------------------------------------*/

static void
_compute_layer_thickness(user_profile_t  *profile)
{
  /* Shorter variable */

  const char *law         = profile->progression_law;
  cs_real_t   progression = profile->progression;
  cs_lnum_t   n_layers    = profile->n_layers;
  cs_real_t   dist_min    = profile->min_dir;
  cs_real_t   dist_max    = profile->max_dir;
  cs_real_t   len_p       = dist_max - dist_min;

  if (strcmp(law, "CONSTANT") == 0) {
    profile->progression = 1.0;
    for (int l_id = 0; l_id < n_layers; l_id++)
      profile->l_thick[l_id] = len_p / ((cs_real_t)n_layers);
  }
  else if (strcmp(law, "GEOMETRIC") == 0) {
    cs_real_t progression_n = pow(progression, n_layers);
    cs_real_t dx0           = len_p * (progression - 1.) / (progression_n - 1.);
    profile->l_thick[0]   = dx0;
    for (int l_id = 0; l_id < n_layers - 1; l_id++)
      profile->l_thick[l_id + 1] = profile->l_thick[l_id] * progression;
  }
  else if (strcmp(law, "PARABOLIC") == 0) {

    /* We need to disguish the case of even or odd number of layers */
    int       is_even = (n_layers % 2 == 0);
    int       np      = 0;
    cs_real_t dx0     = 0.;

    if (is_even) {
      np                       = n_layers / 2;
      cs_real_t progression_np = pow(progression, np);
      dx0 = 0.5 * len_p * (progression - 1.) / (progression_np - 1.);
    }
    else {

      np                        = (n_layers - 1) / 2;
      cs_real_t progression_np  = pow(progression, np);
      cs_real_t progression_np1 = progression_np * progression;
      dx0
        = len_p * (progression - 1.) / (progression_np1 + progression_np - 2.);
    }

    profile->l_thick[0]            = dx0;
    profile->l_thick[n_layers - 1] = dx0;
    for (int l_id = 0; l_id < np - 1; l_id++) {
      profile->l_thick[l_id + 1] = profile->l_thick[l_id] * progression;
      profile->l_thick[n_layers - 1 - 1 - l_id]
        = profile->l_thick[n_layers - 1 - l_id] * progression;
    }

    if (not(is_even))
      profile->l_thick[np] = profile->l_thick[np - 1] * progression;
  }
  else {
    bft_error(
      __FILE__,
      __LINE__,
      0,
      _("Error: Method must be CONSTANT, GEOMETRIC or PARABOLIC for '%s'\n"),
      __func__);
  }
}

/*----------------------------------------------------------------------------
 * Calculate position of the center of each layer in dir direction
 * (absolute and normalized)
 *
 * parameters:
 *   profile    <-> pointer to the current profile structure
 *----------------------------------------------------------------------------*/

static void
_calculate_pos_dir(user_profile_t  *profile)
{
  cs_real_t min_dir = profile->min_dir;
  cs_real_t max_dir = profile->max_dir;

  cs_real_t *l_thick = profile->l_thick;

  profile->pos[0]   = min_dir + 1.0 / 2.0 * l_thick[0];
  profile->pos_n[0] = profile->pos[0] / (max_dir - min_dir);

  for (int layer_id = 1; layer_id < profile->n_layers; layer_id++) {
    profile->pos[layer_id]
      = profile->pos[layer_id - 1]
        + (l_thick[layer_id - 1] + l_thick[layer_id]) * 1.0 / 2.0;
    profile->pos_n[layer_id] = profile->pos[layer_id] / (max_dir - min_dir);
  }
}

/*----------------------------------------------------------------------------
 * Calculate the coordinates of the four edge of a square plane for stl mesh
 * structure STL mesh is triangles base, the face is split into two triangles
 *
 * parameters:
 *   i_v                <-- pointer to i_v
 *   n_v                <-- pointer to normalized dir profile vector
 *   O_coord            <-- point to the center of the plane coordinates
 *   i_size             <-- size of the plane edge (i vector)
 *   j_size             <-- size of the plane edge (j vector)
 *   triangles_coord    --> coordinates of the triangles of the planes
 *
 *----------------------------------------------------------------------------*/

static void
_square_triangle_coords(cs_real_t   i_v[3],
                        cs_real_t   n_v[3],
                        cs_real_t   O_coord[3],
                        cs_real_t   i_size,
                        cs_real_t   j_size,
                        cs_real_3_t triangles_coord[6])
{
  cs_real_t T_i = i_size / 2.0; // Transalaltion compare to origin along i_v
  cs_real_t T_j = j_size / 2.0; // Transalaltion compare to origin along j_v

  cs_real_t j_v[3];
  // j_v is calculated by j_v = n_v x i_v (both normalized vector)
  j_v[0] = n_v[1] * i_v[2] - n_v[2] * i_v[1];
  j_v[1] = n_v[2] * i_v[0] - n_v[0] * i_v[2];
  j_v[2] = n_v[0] * i_v[1] - n_v[1] * i_v[0];

  /* order of the point are important for future calculation of normal of the
     face giben the base, rotation is the opposite of trigonometric
     orientation*/
  for (int k = 0; k < 3; k++)
    triangles_coord[0][k] = O_coord[k] - T_i * i_v[k] - T_j * j_v[k];

  for (int k = 0; k < 3; k++)
    triangles_coord[1][k] = O_coord[k] + T_i * i_v[k] + T_j * j_v[k];

  for (int k = 0; k < 3; k++)
    triangles_coord[2][k] = O_coord[k] - T_i * i_v[k] + T_j * j_v[k];

  for (int k = 0; k < 3; k++)
    triangles_coord[3][k] = O_coord[k] + T_i * i_v[k] - T_j * j_v[k];

  for (int k = 0; k < 3; k++)
    triangles_coord[4][k] = O_coord[k] + T_i * i_v[k] + T_j * j_v[k];

  for (int k = 0; k < 3; k++)
    triangles_coord[5][k] = O_coord[k] - T_i * i_v[k] - T_j * j_v[k];
}

/*----------------------------------------------------------------------------
 * Set seeds for each layer, one on each side of the box
 *
 * Alternative method using cs_geom_closest_point based on coordinates
 * of selected_cells center over all mpi ranks
 *
 * The following members of the structure are updated
 *   mesh_list->seed_coords: coordinates
 *
 * parameters:
 *   profile   <-> pointer to the current profile structure
 *   layer_id  <-- layer number on which the seeds needs to be set
 *----------------------------------------------------------------------------*/

static void
_set_stl_layers_seeds(user_profile_t  *profile,
                      cs_lnum_t        layer_id)
{
  // Shorter profile variable
  cs_stl_mesh_t *stl_mesh = profile->mesh_list[layer_id];

  cs_real_t n_v[3];
  cs_real_t dir_norm = cs_math_3_norm(profile->dir_v);
  for (int k = 0; k < 3; k++)
    n_v[k] = profile->dir_v[k] / dir_norm;

  cs_real_t *i_v = profile->i_v;
  cs_real_t *j_v = profile->j_v;

  cs_lnum_t n_layers = profile->n_layers;
  cs_real_t dist_min = profile->min_dir;
  cs_real_t dist_max = profile->max_dir;

  cs_real_t layer_thickness = profile->l_thick[layer_id];
  cs_real_t l_center_nCoord = profile->pos[layer_id];

  cs_real_t i_translate = (profile->max_i + profile->min_i) / 2.0;
  cs_real_t j_translate = (profile->max_j + profile->min_j) / 2.0;

  // Get mesh quantities
  const cs_mesh_quantities_t *mq         = cs_glob_mesh_quantities;
  const cs_real_t            *cell_cen   = mq->cell_cen;

  // define pointer and variable for cs_selector
  cs_lnum_t  n_selected_cells = 0;
  cs_lnum_t *selected_cells   = NULL;
  // Allocate memory for the cells list which will be populated by cs_selector
  BFT_MALLOC(selected_cells, cs_glob_mesh->n_cells_with_ghosts, cs_lnum_t);

  cs_selector_get_cell_list(
    profile->criteria, &n_selected_cells, selected_cells);

  // Calculate center lower plane (plane 0 in _set_layers_stl_mesh)
  cs_real_3_t O_planes_coord[1];
  cs_real_3_t target_point[2];

  /* For each layer Calculate seed target point as the middle of the segment
   * between layer top and bottom and the mesh sel min/max in the profile
   * direction*/

  for (int k = 0; k < 3; k++)
    O_planes_coord[0][k] = l_center_nCoord * n_v[k]
                           - layer_thickness / 2.0 * n_v[k]
                           + i_translate * i_v[k] + j_translate * j_v[k];

  for (int k = 0; k < 3; k++)
    target_point[0][k]
      = O_planes_coord[0][k]
        - n_v[k] * ((l_center_nCoord - 1.0 / 2.0 * layer_thickness) - dist_min)
            / 2.0;

  for (int k = 0; k < 3; k++)
    target_point[1][k]
      = O_planes_coord[0][k] + n_v[k] * layer_thickness
        + n_v[k] * (dist_max - (l_center_nCoord + 1.0 / 2.0 * layer_thickness))
            / 2.0;

  /* Set an array of selected cell center coordinates */
  cs_real_t *point_coord;
  BFT_MALLOC(point_coord, 3 * n_selected_cells, cs_real_t);

  for (int ii = 0; ii < n_selected_cells; ii++) {
    cs_lnum_t c_id = selected_cells[ii];
    for (int k = 0; k < 3; k++)
      point_coord[3 * ii + k] = cell_cen[3 * c_id + k];
  }

  cs_lnum_t   point_id[2];
  cs_lnum_t   rank_id[2];
  cs_real_3_t closest_point_coord[2];

  /* Initialiaze values of closest_point_coord */
  for (int jj = 0; jj < 2; jj++) {
    for (int k = 0; k < 3; k++)
      closest_point_coord[jj][k] = -1.e20;
  }

  /* For each target point, find the closest cell center in the cell sell over
     all ranks Doing so, the seeds are in the fluid domain */
  for (int jj = 0; jj < 2; jj++) {
    cs_geom_closest_point(n_selected_cells,
                          (const cs_real_3_t *)point_coord,
                          target_point[jj],
                          &point_id[jj],
                          &rank_id[jj]);
    if (point_id[jj] > -1) {
      for (int k = 0; k < 3; k++)
        closest_point_coord[jj][k] = point_coord[3 * point_id[jj] + k];
    }

    cs_parall_max(3, CS_REAL_TYPE, closest_point_coord[jj]);
  }

  BFT_MALLOC(stl_mesh->seed_coords, 3 * (n_layers - 1), cs_real_t);
  cs_lnum_t n_seeds = 0;

  /* Associate the found closest point to stl seeds */
  if (layer_id == 0) {
    for (int k = 0; k < 3; k++) {
      stl_mesh->seed_coords[n_seeds * 3 + k] = closest_point_coord[1][k];
    }
    n_seeds++;
  }
  else if (layer_id == n_layers - 1) {
    for (int k = 0; k < 3; k++) {
      stl_mesh->seed_coords[n_seeds * 3 + k] = closest_point_coord[0][k];
    }
    n_seeds++;
  }
  else {
    for (int k = 0; k < 3; k++) {
      stl_mesh->seed_coords[n_seeds * 3 + k] = closest_point_coord[0][k];
    }
    n_seeds++;
    for (int k = 0; k < 3; k++) {
      stl_mesh->seed_coords[n_seeds * 3 + k] = closest_point_coord[1][k];
    }
    n_seeds++;
  }

  stl_mesh->n_seeds = n_seeds;

  BFT_FREE(selected_cells);

  BFT_FREE(point_coord);
}

/*----------------------------------------------------------------------------
 * Set 1 STL mesh for a given layer profile
 *
 * parameters:
 *   profile  <-> pointer to the current profile structure
 *                The folowing members of the structure are updated:
 *                mesh_list: stl mesh structure (defined in cs_stl_mesh.h)
 *   layer_id <-- layer number on which the seeds needs to be set
 *----------------------------------------------------------------------------*/

static void
_set_layers_stl_mesh(user_profile_t  *profile,
                     cs_lnum_t        layer_id)
{
  cs_stl_mesh_t *stl_mesh = profile->mesh_list[layer_id];
  char           name[10];
  sprintf(name, "%s_%d", "layer", layer_id);
  strncpy(stl_mesh->name, name, 9);
  stl_mesh->name[9] = '\0';

  // Allocate memory for stl mesh
  stl_mesh->n_faces = 6 * 2; // 2 triangles per faces
  BFT_MALLOC(stl_mesh->coords, 3 * stl_mesh->n_faces, cs_real_3_t);

  // Set vector orthogonal spatial system
  cs_real_t n_v[3];
  cs_real_t dir_norm = cs_math_3_norm(profile->dir_v);
  for (int k = 0; k < 3; k++)
    n_v[k] = profile->dir_v[k] / dir_norm;

  /* Shorter variable for new base*/
  cs_real_t *i_v = profile->i_v;
  cs_real_t *j_v = profile->j_v;

  cs_real_t dist_min = profile->min_dir;

  cs_real_t layer_thickness = profile->l_thick[layer_id];
  cs_real_t l_center_nCoord = profile->pos[layer_id];

  cs_real_t i_translate = (profile->max_i + profile->min_i) / 2.0;
  cs_real_t j_translate = (profile->max_j + profile->min_j) / 2.0;

  /* plane size is calculated to be enveloppe of cell selected size*/
  cs_real_t i_size     = (profile->max_i - profile->min_i);
  cs_real_t j_size     = (profile->max_j - profile->min_j);
  cs_real_t plane_size = CS_MAX(i_size, j_size) * 1.5;

  // Calculate center  of each planes
  cs_real_3_t O_planes_coord[6];

  for (int k = 0; k < 3; k++)
    O_planes_coord[0][k] = l_center_nCoord * n_v[k]
                           - layer_thickness / 2.0 * n_v[k] + dist_min * n_v[k]
                           + i_translate * i_v[k] + j_translate * j_v[k];

  for (int k = 0; k < 3; k++)
    O_planes_coord[1][k] = O_planes_coord[0][k] + n_v[k] * layer_thickness;

  for (int k = 0; k < 3; k++)
    O_planes_coord[2][k] = O_planes_coord[0][k] + n_v[k] * layer_thickness / 2.0
                           + i_v[k] * plane_size / 2.0;

  for (int k = 0; k < 3; k++)
    O_planes_coord[3][k] = O_planes_coord[0][k] + n_v[k] * layer_thickness / 2.0
                           - i_v[k] * plane_size / 2.0;

  for (int k = 0; k < 3; k++)
    O_planes_coord[4][k] = O_planes_coord[0][k] + n_v[k] * layer_thickness / 2.0
                           + j_v[k] * plane_size / 2.0;

  for (int k = 0; k < 3; k++)
    O_planes_coord[5][k] = O_planes_coord[0][k] + n_v[k] * layer_thickness / 2.0
                           - j_v[k] * plane_size / 2.0;

  // Get points planes (4 points per plane), 2 triangles
  cs_real_t plane_vtx_coords[6][3];
  cs_lnum_t n_faces_square = 0;

  cs_real_t n_v_plane[3];

  for (int k = 0; k < 3; k++)
    n_v_plane[k] = -n_v[k];

  _square_triangle_coords(i_v,
                          n_v_plane,
                          O_planes_coord[0],
                          plane_size,
                          plane_size,
                          plane_vtx_coords);

  for (int vtx_id = 0; vtx_id < 6; vtx_id++) {
    for (int k = 0; k < 3; k++)
      stl_mesh->coords[n_faces_square * 6 + vtx_id][k]
        = plane_vtx_coords[vtx_id][k];
  }
  n_faces_square++;

  for (int k = 0; k < 3; k++)
    n_v_plane[k] = n_v[k];
  _square_triangle_coords(i_v,
                          n_v_plane,
                          O_planes_coord[1],
                          plane_size,
                          plane_size,
                          plane_vtx_coords);

  for (int vtx_id = 0; vtx_id < 6; vtx_id++) {
    for (int k = 0; k < 3; k++)
      stl_mesh->coords[n_faces_square * 6 + vtx_id][k]
        = plane_vtx_coords[vtx_id][k];
  }
  n_faces_square++;

  for (int k = 0; k < 3; k++)
    n_v_plane[k] = i_v[k];
  _square_triangle_coords(n_v,
                          n_v_plane,
                          O_planes_coord[2],
                          layer_thickness,
                          plane_size,
                          plane_vtx_coords);

  for (int vtx_id = 0; vtx_id < 6; vtx_id++) {
    for (int k = 0; k < 3; k++)
      stl_mesh->coords[n_faces_square * 6 + vtx_id][k]
        = plane_vtx_coords[vtx_id][k];
  }
  n_faces_square++;

  for (int k = 0; k < 3; k++)
    n_v_plane[k] = -i_v[k];

  _square_triangle_coords(n_v,
                          n_v_plane,
                          O_planes_coord[3],
                          layer_thickness,
                          plane_size,
                          plane_vtx_coords);

  for (int vtx_id = 0; vtx_id < 6; vtx_id++) {
    for (int k = 0; k < 3; k++)
      stl_mesh->coords[n_faces_square * 6 + vtx_id][k]
        = plane_vtx_coords[vtx_id][k];
  }
  n_faces_square++;

  for (int k = 0; k < 3; k++)
    n_v_plane[k] = j_v[k];

  _square_triangle_coords(n_v,
                          n_v_plane,
                          O_planes_coord[4],
                          layer_thickness,
                          plane_size,
                          plane_vtx_coords);

  for (int vtx_id = 0; vtx_id < 6; vtx_id++) {
    for (int k = 0; k < 3; k++)
      stl_mesh->coords[n_faces_square * 6 + vtx_id][k]
        = plane_vtx_coords[vtx_id][k];
  }
  n_faces_square++;

  for (int k = 0; k < 3; k++)
    n_v_plane[k] = -j_v[k];

  _square_triangle_coords(n_v,
                          n_v_plane,
                          O_planes_coord[5],
                          layer_thickness,
                          plane_size,
                          plane_vtx_coords);

  for (int vtx_id = 0; vtx_id < 6; vtx_id++) {
    for (int k = 0; k < 3; k++)
      stl_mesh->coords[n_faces_square * 6 + vtx_id][k]
        = plane_vtx_coords[vtx_id][k];
  }
  n_faces_square++;
}

/*----------------------------------------------------------------------------
 * Set 1 MED mesh for a given layer profile
 *
 * parameters:
 *   profile  <-> pointer to the current profile structure
 *                The following members of the structure are updated
 *                med_mesh_struct: medUmesh add in mesh structure
 *   layer_id <-- layer number on which the seeds needs to be set
 *----------------------------------------------------------------------------*/

static void
_set_med_layer_mesh([[maybe_unused]]user_profile_t  *profile,
                    [[maybe_unused]]cs_lnum_t        layer_id)
{
#if defined(HAVE_MEDCOUPLING)
  // Get mesh quantities
  const cs_real_t  *vtx_coord  = cs_glob_mesh->vtx_coord;

  // Set vector orthogonal spatial system
  cs_real_t n_v[3];
  cs_real_t dir_norm = cs_math_3_norm(profile->dir_v);
  for (int k = 0; k < 3; k++)
    n_v[k] = profile->dir_v[k] / dir_norm;

  /* Shorter variable for new base */
  cs_real_t *i_v = profile->i_v;
  cs_real_t *j_v = profile->j_v;

  cs_real_t min_i    = profile->min_i;
  cs_real_t max_i    = profile->max_i;
  cs_real_t min_j    = profile->min_j;
  cs_real_t max_j    = profile->max_j;

  cs_real_t layer_thickness = profile->l_thick[layer_id];
  cs_real_t l_center_nCoord = profile->pos[layer_id];

  /* Set number of cells in each segment */
  int n_x = 10;
  int n_y = 10;
  int n_z = 2;

  double x_min = -(max_i - min_i) / 2.0;
  double x_max = (max_i - min_i) / 2.0;

  double y_min = -(max_j - min_j) / 2.0;
  double y_max = (max_j - min_j) / 2.0;

  double z_min = 0.;
  double z_max = 0.;

  z_min = -layer_thickness / 2.0;
  z_max = layer_thickness / 2.0;

  double *x_coords, *y_coords, *z_coords;
  BFT_MALLOC(x_coords, n_x, double);
  BFT_MALLOC(y_coords, n_y, double);
  BFT_MALLOC(z_coords, n_z, double);

  for (int x_id = 0; x_id < n_x; x_id++)
    x_coords[x_id] = (x_max - x_min) * x_id / (n_x - 1) + x_min;

  for (int y_id = 0; y_id < n_y; y_id++)
    y_coords[y_id] = (y_max - y_min) * y_id / (n_y - 1) + y_min;

  for (int z_id = 0; z_id < n_z; z_id++)
    z_coords[z_id] = (z_max - z_min) * z_id / (n_z - 1) + z_min;

  /* Generate a cartesian mesh in the main base */
  MEDCoupling::DataArrayDouble *arrX = MEDCoupling::DataArrayDouble::New();
  arrX->alloc(n_x, 1);
  std::copy(x_coords, x_coords + n_x, arrX->getPointer());
  arrX->setInfoOnComponent(0, "X [m]");
  MEDCoupling::DataArrayDouble *arrY = MEDCoupling::DataArrayDouble::New();
  arrY->alloc(n_y, 1);
  std::copy(y_coords, y_coords + n_y, arrY->getPointer());
  arrY->setInfoOnComponent(0, "Y [m]");
  MEDCoupling::DataArrayDouble *arrZ = MEDCoupling::DataArrayDouble::New();
  arrZ->alloc(n_z, 1);
  std::copy(z_coords, z_coords + n_z, arrZ->getPointer());
  arrZ->setInfoOnComponent(0, "Z [m]");

  /* Retrieve Umesh preallocated */
  MEDCoupling::MEDCouplingUMesh *Umesh
    = profile->med_mesh_struct->layer_mesh[layer_id];

  MEDCoupling::MEDCouplingCMesh *Cmesh
    = MEDCoupling::MEDCouplingCMesh::New(Umesh->getName());
  Cmesh->setCoords(arrX, arrY, arrZ);

  /* Build unstructured mesh and associate it the pre allocated
   * MEDCouplingUmesh array */
  Umesh = Cmesh->buildUnstructured();

  /* Rotate and translate the mesh */

  cs_real_t i_translate = (profile->max_i + profile->min_i) / 2.0;
  cs_real_t j_translate = (profile->max_j + profile->min_j) / 2.0;

  /* Caculate translation vector */
  cs_real_3_t vec_translate;

  for (int k = 0; k < 3; k++) {
    vec_translate[k]
      = l_center_nCoord * n_v[k] + i_translate * i_v[k] + j_translate * j_v[k];
  }

  /*Rotate the mesh to make it align with i_v, j_v, n_v*/

  cs_real_t center[3] = { 0., 0., 0. };
  cs_real_t Oz[3]     = { 0., 0., 1.0 };
  cs_real_t Ox[3]     = { 1.0, 0., 0. };
  cs_real_t Oy[3]     = { 0., 1.0, 0. };

  /* Temporary base to update base after rotation */
  cs_real_t Ox_1[3] = { 0., 0., 0. };
  cs_real_t Oy_1[3] = { 0., 0., 0. };
  cs_real_t Oz_1[3] = { 0., 0., 0. };

  cs_real_t kv_plane_vec[3] = { 0., 0., 0. }; /* projection of n_v in a plane */
  cs_real_t kv_plane_norm   = 0.;
  cs_real_t rot_angle       = 0.;
  cs_real_t cos_rot_angle   = 0.;
  cs_real_t sin_rot_angle   = 0.;

  /* First make 2 rotations to have the z vectors of the mesh align with n_v */
  /* Rotation around y axis, projection of n_v in Oxz plane */

  for (int k = 0; k < 3; k++) { /*Project n_v in Oxz plane*/
    kv_plane_vec[0] += n_v[k] * Ox[k];
    kv_plane_vec[2] += n_v[k] * Oz[k];
  }

  kv_plane_norm = cs_math_3_norm(kv_plane_vec);

  if (kv_plane_norm > DBL_EPSILON) {
    cos_rot_angle = kv_plane_vec[2] / kv_plane_norm;
    sin_rot_angle = kv_plane_vec[0] / kv_plane_norm;
    rot_angle     = acos(cos_rot_angle);
  }
  else {
    rot_angle = 0.0;
  }

  /* Adapt rotation angle in case of negative sine */
  if (cos_rot_angle > 0)
    if (sin_rot_angle < 0)
      rot_angle -= M_PI + 2.0 * rot_angle;

  if (cos_rot_angle < 0)
    if (sin_rot_angle < 0)
      rot_angle += M_PI - 2.0 * rot_angle;

  Umesh->rotate(center, Oy, rot_angle); /* Rotate counter clock wise */

  /* Update the base after rotation */

  for (int k = 0; k < 3; k++) {
    Ox_1[k] = cos(rot_angle) * Ox[k] - sin(rot_angle) * Oz[k];
  }

  for (int k = 0; k < 3; k++)
    Oz_1[k] = cos(rot_angle) * Oz[k] + sin(rot_angle) * Ox[k];

  for (int k = 0; k < 3; k++)
    Oy_1[k] = Oy[k];

  for (int k = 0; k < 3; k++) {
    Ox[k] = Ox_1[k];
    Oy[k] = Oy_1[k];
    Oz[k] = Oz_1[k];
  }

  /* Rotation around x axis projection of n_v in Oyz plane */
  kv_plane_vec[0] = 0.;
  kv_plane_vec[1] = 0.;
  kv_plane_vec[2] = 0.;
  for (int k = 0; k < 3; k++) { /*Project n_v in Oyz plane*/
    kv_plane_vec[1] += n_v[k] * Oy[k];
    kv_plane_vec[2] += n_v[k] * Oz[k];
  }

  kv_plane_norm = cs_math_3_norm(kv_plane_vec);

  if (kv_plane_norm > DBL_EPSILON) {
    cos_rot_angle = kv_plane_vec[2] / kv_plane_norm;
    sin_rot_angle = kv_plane_vec[1] / kv_plane_norm;
    rot_angle     = acos(cos_rot_angle);
  }
  else
    rot_angle = 0.0;

  /* Adapt rotation angle in case of negative sine */
  if (cos_rot_angle > 0)
    if (sin_rot_angle < 0)
      rot_angle -= M_PI + 2.0 * rot_angle;

  if (cos_rot_angle < 0)
    if (sin_rot_angle < 0)
      rot_angle += M_PI - 2.0 * rot_angle;

  Umesh->rotate(center, Ox, -rot_angle); /* Rotate clock wise */

  /* Update the base after rotation */

  for (int k = 0; k < 3; k++) {
    Oy_1[k] = cos(rot_angle) * Oy[k] - sin(rot_angle) * Oz[k];
  }

  for (int k = 0; k < 3; k++)
    Oz_1[k] = cos(rot_angle) * Oz[k] + sin(rot_angle) * Oy[k];

  for (int k = 0; k < 3; k++)
    Ox_1[k] = Ox[k];

  for (int k = 0; k < 3; k++) {
    Ox[k] = Ox_1[k];
    Oy[k] = Oy_1[k];
    Oz[k] = Oz_1[k];
  }

  /* Rotation around z axis: make Ox aligne with i_v */
  kv_plane_vec[0] = 0.;
  kv_plane_vec[1] = 0.;
  kv_plane_vec[2] = 0.;
  for (int k = 0; k < 3; k++) { /* Project i_v in Oxy plane */
    kv_plane_vec[0] += i_v[k] * Ox[k];
    kv_plane_vec[1] += i_v[k] * Oy[k];
  }

  kv_plane_norm = cs_math_3_norm(kv_plane_vec);

  if (kv_plane_norm > DBL_EPSILON) {
    cos_rot_angle = kv_plane_vec[0] / kv_plane_norm;
    sin_rot_angle = kv_plane_vec[1] / kv_plane_norm;
    rot_angle     = acos(cos_rot_angle);
  }
  else
    rot_angle = 0.0;

  /* Adapt rotation angle in case of negative sine */
  if (cos_rot_angle > 0)
    if (sin_rot_angle < 0)
      rot_angle -= M_PI + 2.0 * rot_angle;

  if (cos_rot_angle < 0)
    if (sin_rot_angle < 0)
      rot_angle += M_PI - 2.0 * rot_angle;

  Umesh->rotate(center, Oz, rot_angle); /* Rotate counter clock wise */

  Umesh->translate(vec_translate); /* Finally translate the mesh */

  profile->med_mesh_struct->layer_mesh[layer_id] = Umesh;

  /* Deallocate data array and Cmesh */
  arrX->decrRef();
  arrY->decrRef();
  arrZ->decrRef();
  Cmesh->decrRef();

  BFT_FREE(x_coords);
  BFT_FREE(y_coords);
  BFT_FREE(z_coords);

#endif
}

/*----------------------------------------------------------------------------
 * Calculate percent of each cells within selection criteria within each layer
 * starting with layer 0, the cells is associated to first layer met in which
 * the cell center coordinates lies.
 *
 * parameters:
 *   profile  <-> pointer to an initalized profile structure.
 *                  The folowing members of the structure are updated
 *                  cells_layer_vol: for each layer, 1 or 0 for each cells
 *                  n_cells: number of cells associated to each layer
 *----------------------------------------------------------------------------*/

static void
_compute_cell_volume_per_layer_basic(user_profile_t  *profile)
{
  // Get mesh quantities
  const cs_mesh_quantities_t *mq       = cs_glob_mesh_quantities;
  const cs_real_t            *cell_cen = mq->cell_cen;
  const cs_lnum_t             n_cells  = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_with_ghosts  = cs_glob_mesh->n_cells_with_ghosts;

  // Profile shorter variables
  cs_lnum_t n_layers = profile->n_layers;

  // Vector profile quantities
  cs_real_t dir_norm = cs_math_3_norm(profile->dir_v);

  cs_real_t dir_normalized[3];
  dir_normalized[0] = profile->dir_v[0] / dir_norm;
  dir_normalized[1] = profile->dir_v[1] / dir_norm;
  dir_normalized[2] = profile->dir_v[2] / dir_norm;
  dir_norm          = 1.0;

  // Define pointer and variable for cs_selector
  cs_lnum_t  n_selected_cells = 0;
  cs_lnum_t *selected_cells   = NULL;
  // Allocate memory for the cells list which will be populated by cs_selector
  BFT_MALLOC(selected_cells, n_cells_with_ghosts, cs_lnum_t);

  cs_selector_get_cell_list(profile->criteria,
                            &n_selected_cells,
                            selected_cells);

  // Reset layer volume arrays

  for (int layer_id = 0; layer_id < n_layers; layer_id++) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      profile->cells_layer_vol[layer_id][c_id] = 0.0;
    }
  }

  for (cs_lnum_t ii = 0; ii < n_selected_cells; ii++) {
    cs_lnum_t c_id = selected_cells[ii];

    cs_real_t x = cell_cen[3 * c_id + 0];
    cs_real_t y = cell_cen[3 * c_id + 1];
    cs_real_t z = cell_cen[3 * c_id + 2];

    // project cell cen vector of dir vector
    cs_real_t dist
      = (x * dir_normalized[0] + y * dir_normalized[1] + z * dir_normalized[2]);

    for (int layer_id = 0; layer_id < n_layers; layer_id++) {

      cs_real_t layer_thickness = profile->l_thick[layer_id];
      cs_real_t l_center_nCoord = profile->pos[layer_id];

      cs_real_t lower_bound = l_center_nCoord - 1.0 / 2.0 * layer_thickness;
      cs_real_t upper_bound = l_center_nCoord + 1.0 / 2.0 * layer_thickness;

      if (dist <= upper_bound && dist >= lower_bound) {
        profile->n_cells[layer_id] += 1;
        profile->cells_layer_vol[layer_id][c_id] = 1.0;

        break; // stop the loop, cells cannot be associated to another layer
      }
    }

  } // end of for loop for layer cell vol / weigth calculation

  BFT_FREE(selected_cells);
}

/*----------------------------------------------------------------------------
 * Calculate percent of each cells within selection criteria within each layer
 * based on cs_stl_mesh porosity calculation
 *
 * parameters:
 *   profile  <-> pointer to an initalized profile structure.
 *                  The  folowing members of the structure are updated
 *                  cells_layer_vol: for each layer, percent of cell in
 *                  layer_id is updated.
 *----------------------------------------------------------------------------*/

static void
_compute_cell_volume_per_layer_stl(user_profile_t  *profile)
{
  /* Get Mesh quantities */
  const cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_cells             = cs_glob_mesh->n_cells;

  // Profile shorter variables
  cs_lnum_t n_layers = profile->n_layers;

  cs_real_t *cells_l_id_vol = NULL;
  BFT_MALLOC(cells_l_id_vol, n_cells_with_ghosts, cs_real_t);

  // Define pointer and variable for cs_selector
  cs_lnum_t  n_selected_cells = 0;
  cs_lnum_t *selected_cells   = NULL;
  // Allocate memory for the cells list which will be populated by cs_selector
  BFT_MALLOC(selected_cells, n_cells_with_ghosts, cs_lnum_t);

  cs_selector_get_cell_list(
    profile->criteria, &n_selected_cells, selected_cells);

  // Reset layer volume arrays

  for (int layer_id = 0; layer_id < n_layers; layer_id++) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      profile->cells_layer_vol[layer_id][c_id] = 0.0;
    }
  }

  for (int s_id = 0; s_id < profile->n_layers; s_id++) {

    for (cs_lnum_t c_id = 0; c_id < n_cells_with_ghosts; c_id++)
      cells_l_id_vol[c_id] = 0.0;

    cs_stl_mesh_t *stl_mesh = profile->mesh_list[s_id];
    /* Compute porisity associated to each layer using cs_stl_mesh features */
    cs_stl_compute_porosity(stl_mesh, cells_l_id_vol, NULL);

    for (cs_lnum_t ii = 0; ii < n_selected_cells; ii++) {
      cs_lnum_t c_id = selected_cells[ii];
      profile->cells_layer_vol[s_id][c_id] = 1.0 - cells_l_id_vol[c_id];
    }
  }

  BFT_FREE(cells_l_id_vol);
  BFT_FREE(selected_cells);
}

/*----------------------------------------------------------------------------
 * Compute intersection between a layer and med mesh
 * (method derived from cs_medcoupling_intersector.cxx)
 *
 * parameters:
 *   med_t         <-- pointer to a med mesh structure
 *   vol_intersect <-> pre allocated array to store result of intersection
 *                     (size: n_cells on local MPI rank)
 *   layer_id      <-- number of the layer to be intersected with cs local mesh
 *----------------------------------------------------------------------------*/

static void
_compute_intersection_volume_med
(
 [[maybe_unused]]user_profile_med_t *med_t,
 [[maybe_unused]]cs_real_t          *vol_intersect,
 [[maybe_unused]]cs_lnum_t           layer_id
)
{
#if defined(HAVE_MEDCOUPLING)
  cs_lnum_t n_elts = med_t->local_mesh->n_elts;

  /* initialize the pointer */
  for (cs_lnum_t c_id = 0; c_id < cs_glob_mesh->n_cells; c_id++)
    vol_intersect[c_id] = 0.;

  /* Matrix for the target mesh */
  MEDCouplingNormalizedUnstructuredMesh<3, 3> tMesh_wrapper
    (med_t->local_mesh->med_mesh);

  /* Matrix for the layer mesh, based on the bbox of the target mesh */
  const cs_real_t *bbox = med_t->local_mesh->bbox;

  const DataArrayIdType *subcells
    = med_t->layer_mesh[layer_id]->getCellsInBoundingBox(bbox, 1.05);

  MEDCouplingNormalizedUnstructuredMesh<3, 3> sMesh_wrapper
    (med_t->layer_mesh[layer_id]->buildPartOfMySelf
                                   (subcells->begin(), subcells->end(), true));
  /* Compute the intersection matrix between source and target meshes */
  std::vector<std::map<mcIdType, double> > mat;
  INTERP_KERNEL::Interpolation3D           interpolator;

  interpolator.interpolateMeshes(sMesh_wrapper, tMesh_wrapper, mat, "P0P0");

  /* Loop on the different elements of the target mesh.
   * For each element, we sum all intersected volumes to retrieve the total
   * intersected volume per cell.
   * The iterator map contains two elements:
   * -> first : which is the index of the intersected cell in source mesh
   * -> second: which the intersection volume
   */
  const cs_lnum_t *connec = med_t->local_mesh->new_to_old;
  for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
    cs_lnum_t c_id = connec[e_id];
    for (std::map<mcIdType, double>::iterator it = mat[e_id].begin();
         it != mat[e_id].end();
         ++it) {
      vol_intersect[c_id] += it->second;
    }
  }
#endif
}

/*----------------------------------------------------------------------------
 * Calculate percent of each cells within selection criteria within each layer
 * based on med coupling intersection (method derived from
 * cs_medcoupling_intersector.cxx
 *
 * parameters:
 *   profile <-> pointer to an initialized profile structure
 *                 The folowing members of the structure are updated
 *                 cells_layer_vol: for each layer, percent of cell in
 *                 layer_id is updated
 *----------------------------------------------------------------------------*/

static void
_compute_cell_vol_per_layer_med([[maybe_unused]] user_profile_t *profile)
{
#if defined(HAVE_MEDCOUPLING)

  /* Get Mesh quantities */
  const cs_lnum_t  n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t  n_cells             = cs_glob_mesh->n_cells;
  const cs_real_t *cell_vol            = cs_glob_mesh_quantities->cell_vol;

  // profile shorter variables
  cs_lnum_t n_layers = profile->n_layers;

  cs_real_t *cells_l_id_vol = NULL;
  BFT_MALLOC(cells_l_id_vol, n_cells_with_ghosts, cs_real_t);

  // define pointer and variable for cs_selector
  cs_lnum_t  n_selected_cells = 0;
  cs_lnum_t *selected_cells   = NULL;
  // Allocate memory for the cells list which will be populated by cs_selector
  BFT_MALLOC(selected_cells, n_cells_with_ghosts, cs_lnum_t);

  cs_selector_get_cell_list(
    profile->criteria, &n_selected_cells, selected_cells);

  // reset layer volume arrays

  for (int layer_id = 0; layer_id < n_layers; layer_id++) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      profile->cells_layer_vol[layer_id][c_id] = 0.0;
    }
  }

  for (int s_id = 0; s_id < profile->n_layers; s_id++) {

    for (cs_lnum_t c_id = 0; c_id < n_cells_with_ghosts; c_id++)
      cells_l_id_vol[c_id] = 0.0;

    user_profile_med_t *med_t = profile->med_mesh_struct;
    /* Compute cells intersection */
    _compute_intersection_volume_med(med_t, cells_l_id_vol, s_id);
    /* Update intersection for profile */
    for (cs_lnum_t ii = 0; ii < n_selected_cells; ii++) {
      cs_lnum_t c_id = selected_cells[ii];
      profile->cells_layer_vol[s_id][c_id]
        = cells_l_id_vol[c_id] / cell_vol[c_id];
    }
  }

  BFT_FREE(cells_l_id_vol);
  BFT_FREE(selected_cells);

#endif
}

/*----------------------------------------------------------------------------
 * Free memory of all allocated memory for a profile variables
 *
 * parameters:
 *   profile <-> pointer to the current profile structure
 *----------------------------------------------------------------------------*/

static void
_free_profile_all(user_profile_t *profile)
{
  // free stl meshes
  for (int s_id = 0; s_id < profile->n_layers; s_id++) {
    cs_stl_mesh_t *stl_mesh = profile->mesh_list[s_id];
    BFT_FREE(stl_mesh->coords);
    BFT_FREE(stl_mesh->seed_coords);
    BFT_FREE(stl_mesh->ext_mesh);
    BFT_FREE(stl_mesh);
  }

  BFT_FREE(profile->mesh_list);

#if defined(HAVE_MEDCOUPLING)
  /* Free med mesh if created */
  int test_med = strcmp(profile->intersect_method, "MEDCOUPLING");
  user_profile_med_t *med_t = profile->med_mesh_struct;
  if (test_med == 0) {
    for (int m_id = 0; m_id < profile->n_layers; m_id++) {
      med_t->layer_mesh[m_id]->decrRef();
    }
    BFT_FREE(med_t->layer_mesh);
    // Mesh will deallocated afterwards since it can be shared
    med_t->local_mesh = NULL;
    BFT_FREE(profile->med_mesh_struct);
  }

#endif

  /* Free layer histogram */
  for (int l_id = 0; l_id < profile->n_layers; l_id++) {
    user_histogram_t *histogram = profile->histogram_list[l_id];
    user_destroy_histogram(histogram);
  }
  BFT_FREE(profile->histogram_list);

  BFT_FREE(profile->l_thick);
  BFT_FREE(profile->pos);
  BFT_FREE(profile->pos_n);
  BFT_FREE(profile->n_cells);
  BFT_FREE(profile->weigth);
  BFT_FREE(profile->mean_f);
  BFT_FREE(profile->mean_f_n);
  BFT_FREE(profile->sd_f);
  BFT_FREE(profile->sd_f_n);

  for (int layer_id = 0; layer_id < profile->n_layers; layer_id++)
    BFT_FREE(profile->cells_layer_vol[layer_id]);
  BFT_FREE(profile->cells_layer_vol);
}

/*----------------------------------------------------------------------------
 * Output the setup of the profile into a log file stored in
 * ./profiles/<profile name>
 *
 * parameters:
 *   profile      <-- pointer to the current profile structure
 *   dirname      <-- pointer to directory output path
 *----------------------------------------------------------------------------*/

static void
_output_profile_setup_log(user_profile_t  *profile,
                          const char      *dirname)
{
  /* Profile shorter variables */

  char filename[80];

  sprintf(filename, "%s/profile.log", dirname);
  FILE *ptrSetupLog = NULL;

  ptrSetupLog = fopen(filename,
                      "w+"); /* create the file (should not exist
                                when this function is called) */

  fprintf(ptrSetupLog,
          "\n"
          "              ----------------------------------------------------\n"
          "                             %s \n"
          "              ----------------------------------------------------\n",
          profile->name);

  fprintf(ptrSetupLog,
          "\n\n"
          " Global informations\n"
          " ------------\n\n"
          " field: %s\n"
          " cells selection: %s\n"
          " number of layers: %d\n"
          " normal profile direction: {%.2f,%.2f,%.2f}\n"
          " minimum mesh distance along profile direction: %.2f\n"
          " maximum mesh distance along profile direction: %.2f\n"
          " Weigthed: %s\n"
          " layer progression law: %s\n"
          " progression coefficient: %.3f\n"
          " Mesh Intersection Method: %s\n",
          profile->field,
          profile->criteria,
          profile->n_layers,
          profile->dir_v[0],
          profile->dir_v[1],
          profile->dir_v[2],
          profile->min_dir,
          profile->max_dir,
          profile->weighted,
          profile->progression_law,
          profile->progression,
          profile->intersect_method);

  fprintf(ptrSetupLog, " layer thickness: [%.3f", profile->l_thick[0]);

  for (int l_id = 1; l_id < profile->n_layers; l_id++)
    fprintf(ptrSetupLog, ",%.3f", profile->l_thick[l_id]);

  fprintf(ptrSetupLog, "]\n");

  if (strcmp(profile->intersect_method, "STL") == 0) {
    fprintf(ptrSetupLog,
            "\n\n"
            " STL Layers mesh informations\n"
            " --------------\n");

    for (int layer_id = 0; layer_id < profile->n_layers; layer_id++) {
      cs_stl_mesh_t *stl_mesh = profile->mesh_list[layer_id];
      fprintf(ptrSetupLog,
              "\n"
              "    %s\n",
              stl_mesh->name);

      for (int seed_id = 0; seed_id < stl_mesh->n_seeds; seed_id++) {
        cs_real_t x = stl_mesh->seed_coords[3 * seed_id + 0];
        cs_real_t y = stl_mesh->seed_coords[3 * seed_id + 1];
        cs_real_t z = stl_mesh->seed_coords[3 * seed_id + 2];

        fprintf(ptrSetupLog,
                "      seed %d coords :{%.2f,%.2f,%.2f}\n",
                seed_id,
                x,
                y,
                z);
      }
    }
  }

  fclose(ptrSetupLog);
}

/*----------------------------------------------------------------------------
 * Output the current profile into a log file stored in ./profiles/<profile name>
 *
 * parameters:
 *   profile      <-- pointer to the current profile structure
 *   dirname      <-- pointer to directory output path
 *----------------------------------------------------------------------------*/

static void
_output_profile_log(user_profile_t  *profile,
                    const char      *dirname)
{
  /* Profile shorter variables*/
  cs_lnum_t   n_layers = profile->n_layers;
  const char *field    = profile->field;

  char filename[80];
  sprintf(filename, "%s/profile.log", dirname);
  FILE *ptrLog = NULL;

  ptrLog = fopen(filename,
                 "a+"); // File already created when this function is called

  cs_time_step_t *ts = cs_get_glob_time_step();

  fprintf(ptrLog,
          "\n\n"
          " Time step: %d - Simulation time: %.2f\n"
          " --------------\n"
          "      min %s is: %.2f - max %s is :%.2f wihtin cells selection\n\n",
          ts->nt_cur,
          ts->t_cur,
          profile->field,
          profile->min_field,
          profile->field,
          profile->max_field);

  cs_real_t total_weigth = 0.0;
  /* Output absolute profile values*/
  for (int layer_id = 0; layer_id < n_layers; layer_id++) {
    total_weigth += profile->weigth[layer_id];
    fprintf(ptrLog,
            "      mean %s of layer %d (pos = %.2f) is: %.2f - sd %.2f - Q1: "
            "%.2f - Q2: %.2f - Q3: %.2f - weigth: %.2f \n",
            field,
            layer_id,
            profile->pos[layer_id],
            profile->mean_f[layer_id],
            profile->sd_f[layer_id],
            profile->histogram_list[layer_id]->Q1,
            profile->histogram_list[layer_id]->Q2,
            profile->histogram_list[layer_id]->Q3,
            profile->weigth[layer_id]);
  }

  fprintf(ptrLog, "\n");

  /* Output normalized profile value*/

  for (int layer_id = 0; layer_id < n_layers; layer_id++) {
    fprintf(ptrLog,
            "      mean_n %s of layer %d (pos_n = %.2f) is: %.2f - sd_n %.2f "
            "- weigth: %.2f \n",
            field,
            layer_id,
            profile->pos_n[layer_id],
            profile->mean_f_n[layer_id],
            profile->sd_f_n[layer_id],
            profile->weigth[layer_id]);
  }

  fprintf(ptrLog, "\n");

  /* Compare total layer weigth and cell sel weigth
     Useful in Case stl mesh use for cell layer vol calculation */
  fprintf(ptrLog,
          "      total layer weigth is: %.2f \n"
          "      total selected cells weigth: %.2f\n",
          total_weigth,
          profile->sel_cells_weigth);

  fclose(ptrLog);
}

/*----------------------------------------------------------------------------
 * Output the current profile into a csv file stored in ./profiles/<profile name>
 *
 * parameters:
 *   profile      <-- pointer to the current profile structure
 *   dirname      <-- pointer to directory output path
 *----------------------------------------------------------------------------*/

static void
_output_profile_values_csv(user_profile_t  *profile,
                           const char      *dirname)
{
  /* Profile shorter variables */
  cs_lnum_t       n_layers = profile->n_layers;
  cs_time_step_t *ts       = cs_get_glob_time_step();

  cs_lnum_t   n_l_header = 7;
  const char *layer_header[]
    = { "pos", "weigth", "mean", "sd", "pos_n", "mean_n", "sd_n" };

  cs_real_t *profile_var[]
    = { profile->pos,   profile->weigth, profile->mean_f,
        profile->sd_f,  profile->pos_n,  profile->mean_f_n,
        profile->sd_f_n };

  char filename[80];

  sprintf(filename, "%s/results_profile.csv", dirname);
  FILE *ptrResults = NULL;

  /* Check if the file exist*/
  ptrResults = fopen(filename, "r+"); /* Create the file (should not exist
                                         when this function is called) */

  if (ptrResults == NULL) {
    /* Create a file and write header*/
    ptrResults = fopen(filename, "w+");
    fprintf(ptrResults, "time_step, time");

    for (int layer_id = 0; layer_id < n_layers; layer_id++) {
      for (int h_id = 0; h_id < n_l_header; h_id++)
        fprintf(ptrResults, ",l_%d_%s", layer_id, layer_header[h_id]);
    }

    fprintf(ptrResults, ", total_l_w, cells_sel_w\n");
  }

  fclose(ptrResults);

  ptrResults = fopen(filename, "a+"); /* Add current result in the
                                         file re using same pattern as above */

  fprintf(ptrResults, "%d,%f", ts->nt_cur, ts->t_cur);

  for (int layer_id = 0; layer_id < n_layers; layer_id++) {
    for (int h_id = 0; h_id < n_l_header; h_id++)
      fprintf(ptrResults, ",%f", profile_var[h_id][layer_id]);
  }

  cs_real_t total_l_weigth = 0.0;
  for (int layer_id = 0; layer_id < n_layers; layer_id++)
    total_l_weigth += profile->weigth[layer_id];

  fprintf(ptrResults, ",%f,%f\n", total_l_weigth, profile->sel_cells_weigth);

  fclose(ptrResults);
}

/*----------------------------------------------------------------------------
 * Output the stl mesh into stl binaries stored in ./profiles/<profile name>
 *
 * parameters:
 *   profile      <-- pointer to the current profile structure
 *   dirname      <-- pointer to directory output path
 *----------------------------------------------------------------------------*/

static void
_user_output_profile_stl_mesh(user_profile_t  *profile,
                              const char      *dirname)
{
  if (cs_glob_rank_id <= 0) {
    for (int s_id = 0; s_id < profile->n_layers; s_id++) {
      char           outfile[200];
      cs_stl_mesh_t *stl_mesh = profile->mesh_list[s_id];
      sprintf(outfile, "%s/%s.stl", dirname, stl_mesh->name);
      cs_stl_file_write(stl_mesh, outfile);
    }
  }
}

/*----------------------------------------------------------------------------
 * Output the med mesh into .med files stored in ./profiles/<profile name>
 *
 * parameters:
 *   profile      <-- pointer to the current profile structure
 *   dirname      <-- pointer to directory output path
 *----------------------------------------------------------------------------*/

static void
_output_profile_med_mesh([[maybe_unused]] user_profile_t  *profile,
                         [[maybe_unused]] const char      *dirname)
{
#if defined(HAVE_MEDCOUPLING_LOADER)
  if (cs_glob_rank_id <= 0) {
    for (int m_id = 0; m_id < profile->n_layers; m_id++) {

      char name[10];
      sprintf(name, "%s_%d", "layer", m_id); /*Generete the same name */
      char outfile[200];
      sprintf(outfile, "%s/%s.med", dirname, name);

      MEDCoupling::WriteUMesh(outfile, profile->med_mesh_struct->layer_mesh[m_id],
                              true);
    }
  }
#endif
}

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return a pointer to a profile structure, null if the profile does not exist
 *
 * parameters:
 *   name <-- pointer to profile name
 *----------------------------------------------------------------------------*/

user_profile_t *
user_profile_get_by_name(const char  *name)
{
  user_profile_t *ptr = NULL;
  for (int p_id = 0; p_id < _profile_list.n_profiles; p_id++) {
    user_profile_t *profile = _profile_list.profile_list[p_id];
    int             test      = strcmp(profile->name, name);
    if (test == 0)
      ptr = profile;
  }

  return ptr;
}

/*----------------------------------------------------------------------------
 * Return a pointer to a new and initialiazed profile structure. If profile
 * already exist, return pointer of profile and skip initialiazation
 *
 * parameters:
 *   name               <-- pointer to profile name
 *   field              <-- pointer to field name to be studied in profile
 *   criteria           <-- pointer to mesh cells sel criteria (on which
 *                          profile will be computed)
 *   dir                <-- pointer to profile direction vector
 *                          (does not need to be normalized)
 *   n_layers           <-- layer number of the profile
 *   mass_weight        <-- if true mean profile will be weigthed by cell
 *                          mass, volume otherwise
 *   intersect_method   <-- method to be used for intersect volumes
 *                          (Basic, STL or MEDCOUPLING)
 *----------------------------------------------------------------------------*/

user_profile_t *
user_create_profile(const char  *name,
                    const char  *field,
                    const char  *criteria,
                    cs_real_t    dir[3],
                    cs_lnum_t    n_layers,
                    const char  *progression_law,
                    cs_real_t    progression,
                    const char  *weighted,
                    const char  *intersect_method)
{
  /* First check if profile exists */

  user_profile_t *profile = user_profile_get_by_name(name);

  if (profile != NULL) {
    bft_printf("profile %s already exist - profile not created \n", name);
    return profile;
  }

  _profile_list.n_profiles++;
  BFT_REALLOC(_profile_list.profile_list,
              _profile_list.n_profiles,
              user_profile_t *);

  const cs_lnum_t   n_cells  = cs_glob_mesh->n_cells;

  // Initialize and allocate memory for profile
  BFT_MALLOC(profile, 1, user_profile_t);

  char *_dummy = const_cast<char *>(profile->name);
  strncpy(_dummy, name, 511);
  _dummy[511] = '\0';

  _dummy = const_cast<char *>(profile->field);
  strncpy(_dummy, field, 511);
  _dummy[511] = '\0';

  _dummy = const_cast<char *>(profile->criteria);
  strncpy(_dummy, criteria, 511);
  _dummy[511] = '\0';

  profile->n_layers        = n_layers;

  _dummy = const_cast<char *>(profile->progression_law);
  strncpy(_dummy, progression_law, 511);
  _dummy[511] = '\0';

  profile->progression     = progression;

  for (int k = 0; k < 3; k++)
    profile->dir_v[k] = dir[k];

  if (   strcmp(weighted, "MASS") != 0 && strcmp(weighted, "VOLUME") != 0
      && strcmp(weighted, "NO") != 0)
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: Weigthed param must be MASS, VOLUME or NO for '%s'\n"),
              __func__);


  _dummy = const_cast<char *>(profile->weighted);
  strncpy(_dummy, weighted, 511);
  _dummy[511] = '\0';

  /* In case of Medcoupling method, check CS has been compiled with
   * MEDCouppling*/
  int test_med = strcmp(intersect_method, "MEDCOUPLING");
  if (test_med == 0) {
#if !defined(HAVE_MEDCOUPLING)
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: This intersection method cannot be called without "
                "MEDCoupling support.\n"));
#endif
  }
  _dummy = const_cast<char *>(profile->intersect_method);
  strncpy(_dummy, intersect_method, 511);
  _dummy[511] = '\0';

  BFT_MALLOC(profile->l_thick, n_layers, cs_real_t);
  BFT_MALLOC(profile->pos, n_layers, cs_real_t);
  BFT_MALLOC(profile->pos_n, n_layers, cs_real_t);
  BFT_MALLOC(profile->n_cells, n_layers, cs_lnum_t);
  BFT_MALLOC(profile->weigth, n_layers, cs_real_t);
  BFT_MALLOC(profile->mean_f, n_layers, cs_real_t);
  BFT_MALLOC(profile->mean_f_n, n_layers, cs_real_t);
  BFT_MALLOC(profile->sd_f, n_layers, cs_real_t);
  BFT_MALLOC(profile->sd_f_n, n_layers, cs_real_t);
  BFT_MALLOC(profile->cells_layer_vol, n_layers, cs_real_t *);
  for (int layer_id = 0; layer_id < n_layers; layer_id++)
    BFT_MALLOC(profile->cells_layer_vol[layer_id], n_cells, cs_real_t);

  BFT_MALLOC(profile->mesh_list, n_layers, cs_stl_mesh_t *);
  for (int layer_id = 0; layer_id < n_layers; layer_id++)
    BFT_MALLOC(profile->mesh_list[layer_id], 1, cs_stl_mesh_t);

  /* Allocate mesh med struct*/
  if (test_med == 0)
    profile->med_mesh_struct = _allocate_med_mesh_struct(n_layers);

  /* Allocate array of histogram */
  BFT_MALLOC(profile->histogram_list, n_layers, user_histogram_t *);
  for (int layer_id = 0; layer_id < n_layers; layer_id++) {
    char pname[100];
    sprintf(pname, "layer_%d", layer_id);
    user_histogram_t *histogram
      = user_create_histogram(pname, profile->field, 300); /* n_bins_max */

    profile->histogram_list[layer_id] = histogram;
  }

  // Initialiaze value of struct profile
  profile->sel_cells_weigth = 0.0;
  profile->min_field        = 0.0;
  profile->max_field        = 0.0;
  for (int layer_id = 0; layer_id < n_layers; layer_id++) {
    profile->pos[layer_id]      = 0.0;
    profile->pos_n[layer_id]    = 0.0;
    profile->n_cells[layer_id]  = 0;
    profile->weigth[layer_id]   = 0.0;
    profile->mean_f[layer_id]   = 0.0;
    profile->mean_f_n[layer_id] = 0.0;
    profile->sd_f[layer_id]     = 0.0;
    profile->sd_f_n[layer_id]   = 0.0;
    profile->l_thick[layer_id]  = 0.0;

    cs_stl_mesh_t *stl_mesh = profile->mesh_list[layer_id];

    memset(stl_mesh->header, 0, 80);
    stl_mesh->n_faces     = 0;
    stl_mesh->coords      = NULL;
    stl_mesh->n_seeds     = 0;
    stl_mesh->seed_coords = NULL;
    stl_mesh->is_porous   = false;
    stl_mesh->ext_mesh    = NULL;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      profile->cells_layer_vol[layer_id][c_id] = 0.0;
  }

  /* Compute min and max selected mesh part for dir_v noramlized */
  _calculate_min_max_dir(profile);

  /* Compute layer thickness*/
  _compute_layer_thickness(profile);

  /* Compute position of the center of the different layer along dir_v
   * noramlized */
  _calculate_pos_dir(profile);

  // Generate layer stl mesh
  for (int layer_id = 0; layer_id < n_layers; layer_id++) {
    _set_layers_stl_mesh(profile, layer_id);
    _set_stl_layers_seeds(profile, layer_id);
  }

  // Generate Med Mesh if MEDCOUPLING method is chosen
  if (test_med == 0) {
    for (int layer_id = 0; layer_id < n_layers; layer_id++)
      _set_med_layer_mesh(profile, layer_id);
  }

  _profile_list.profile_list[_profile_list.n_profiles - 1] = profile;

  return profile;
}

/*----------------------------------------------------------------------------
 * Return an allocated histogram structure. Function is public for use histogram
 * oustide of mean profile context
 *
 * parameters :
 *   name        <-- name of the histogram
 *   field       <-- field on which the histogram is compute
 *   n_bins_max  <-- maximum number of bin of the bins
 *----------------------------------------------------------------------------*/

user_histogram_t *
user_create_histogram(char        *name,
                      const char  *field,
                      cs_lnum_t    n_bins_max)
{
  user_histogram_t *histogram;

  BFT_MALLOC(histogram, 1, user_histogram_t);
  BFT_MALLOC(histogram->h_i, n_bins_max, cs_real_t);
  BFT_MALLOC(histogram->l_i, n_bins_max, cs_real_t);
  BFT_MALLOC(histogram->c_i, n_bins_max, cs_real_t);
  BFT_MALLOC(histogram->name, 105, char);

  /* Populate histogram general informations */
  sprintf(histogram->name, "%s", name);
  histogram->field      = field;
  histogram->n_bins_max = n_bins_max;

  /* Initalize historgam values */

  for (int b_id = 0; b_id < n_bins_max; b_id++) {
    histogram->h_i[b_id] = 0.0;
    histogram->l_i[b_id] = 0.0;
    histogram->c_i[b_id] = 0.0;
  }

  histogram->n_bins = 0;
  histogram->bandwidth = 0.0;
  histogram->mean      = 0.0;
  histogram->sd        = 0.0;
  histogram->min       = 0.0;
  histogram->max       = 0.0;
  histogram->Q1        = 0.0;
  histogram->Q2        = 0.0;
  histogram->Q3        = 0.0;
  histogram->ptot      = 0.0;

  return histogram;
}

/*----------------------------------------------------------------------------
 * Compute a histogram with a given sample
 *
 * parameters :
 *   histogram    <-> name of the histogram
 *      h_id        --> compute h_i
 *      l_i         --> compute l_i
 *      c_i         --> compute c_i
 *      mean        --> compute sample mean
 *      sd          --> compute sample standar deviation
 *      Q1,Q2,Q3    --> compute sample quartiles
 *      probatot    --> compute sum of h_i (has to 1.0)
 *
 *   sample         <-- 1d sample on which the histogram is build
 *   weights        <-- weight associated to each sample element
 *   n_elts_sample  <-- number of elements in the provided sample
 *----------------------------------------------------------------------------*/

void
user_histogram_compute(user_histogram_t  *histogram,
                       cs_real_t         *sample,
                       cs_real_t         *weights,
                       cs_lnum_t          n_elts_sample)
{
  cs_real_t moment[]  = {0.0, 0.0};
  cs_real_t min_max[] = {0.0, 0.0};

  _compute_sample_moment(sample, weights, n_elts_sample, moment, min_max);

  histogram->mean = moment[0];
  histogram->sd   = moment[1];
  histogram->min  = min_max[0];
  histogram->max  = min_max[1];

  /*weight are assumed normalized*/
  _compute_histogram(histogram, sample, weights, n_elts_sample);
}

/*----------------------------------------------------------------------------
 * Free histogram structure
 *
 * parameters :
 *   histogram    <-> name of the histogram
 *----------------------------------------------------------------------------*/

void
user_destroy_histogram(user_histogram_t  *histogram)
{
  BFT_FREE(histogram->name);
  BFT_FREE(histogram->h_i);
  BFT_FREE(histogram->l_i);
  BFT_FREE(histogram->c_i);

  BFT_FREE(histogram);
}

/*----------------------------------------------------------------------------
 * Calculate percent of each cells within selection criteria within each layer
 * based on chose method (stl, basic or medcoupling)
 *
 * parameters:
 *   profile <-> pointer to an initialized profile structure
 *               The folowing members of the structure are updated:
 *               cells_layer_vol: for each layer, percent of cell
 *               in layer_id is updated
 *----------------------------------------------------------------------------*/

void
user_compute_cell_volume_per_layer(user_profile_t *profile)
{
  if (strcmp(profile->intersect_method, "MEDCOUPLING") == 0) {
#if !defined(HAVE_MEDCOUPLING)
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: This intersection method cannot be called without "
                "MEDCoupling support.\n"));
#else
    _compute_cell_vol_per_layer_med(profile);
#endif
  }
  else if (strcmp(profile->intersect_method, "STL") == 0) {
    _compute_cell_volume_per_layer_stl(profile);
  }
  else if (strcmp(profile->intersect_method, "BASIC") == 0) {
    _compute_cell_volume_per_layer_basic(profile);
  }
  else {
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: Method must be BASIC, STL or MEDCOUPLING for '%s'\n"),
              __func__);
  }
}

/*----------------------------------------------------------------------------
 * Calculate mean and standard deviation for the selected field for for the
 * cells lying in the different layers. A normlized mean and normalized
 * standard deviation is also computed. 0: minimum of the field met in the
 * selection criteria 1: maximum of the field met in the selection criteria
 * if field dimension is higher than 1, its norm norm is used for field value
 *
 * Function to be used once cells_layer_vol has been filled by:
 *   user_compute_cell_volume_per_layer()
 *
 * parameters:
 *   profile <-> pointer to a profile structure
 *               The folowing members of the structure are updated:
 *               mean_f: for each layer, mean value of the field
 *               sd_f: for each layer, field standard deviation
 *               mean_f_n: for each layer, normalized field mean
 *               sd_f_n: for each layer, standard deviation of normalized field
 *               histogram_list: each histogram is updated
 *----------------------------------------------------------------------------*/

void
user_profile_compute(user_profile_t  *profile)
{
  // Get mesh quantities

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  // profile shorter variables
  cs_lnum_t n_layers = profile->n_layers;

  /*Histogram: part to be be merged with upstream to avoid double computation*/

  /*generate sample and weights*/
  cs_real_t *sample, *weights;

  cs_real_t min_field = DBL_MAX;
  cs_real_t max_field = DBL_MIN;

  BFT_MALLOC(sample, n_cells, cs_real_t);
  BFT_MALLOC(weights, n_cells, cs_real_t);
  cs_lnum_t n_elts_sample = 0;

  for (int l_id = 0; l_id < n_layers; l_id++) {
    user_histogram_t *histogram = profile->histogram_list[l_id];

    /*Create a sample based on elements belonging to each layer*/
    _create_1d_sample_(profile, &n_elts_sample, sample, weights, l_id);

    /*Compute each layer weigth - useful for checking relevant intersection
     * calculation*/
    profile->weigth[l_id] = 0.0;

    for (cs_lnum_t iel = 0; iel < n_elts_sample; iel++)
      profile->weigth[l_id] += weights[iel];

    cs_parall_sum(1, CS_REAL_TYPE, &profile->weigth[l_id]);

    user_histogram_compute(histogram, sample, weights, n_elts_sample);

    /*Update profile values with calculated for each histogram layer*/
    profile->mean_f[l_id] = histogram->mean;
    profile->sd_f[l_id]   = histogram->sd;

    min_field = CS_MIN(min_field, histogram->min);
    max_field = CS_MAX(max_field, histogram->max);
  }

  /*update min and max field over all layer and all MPI ranks*/
  cs_parall_min(1, CS_REAL_TYPE, &min_field);
  cs_parall_max(1, CS_REAL_TYPE, &max_field);

  profile->min_field = min_field;
  profile->max_field = max_field;

  /*Compute normalized mean and sd*/

  for (int l_id = 0; l_id < n_layers; l_id++) {
    profile->mean_f_n[l_id] = 0.0;
    profile->sd_f_n[l_id]   = 0.0;

    if ((max_field - min_field) > DBL_EPSILON) {
      /*Create a sample based on elements belonging to each layer*/
      _create_1d_sample_(profile, &n_elts_sample, sample, weights, l_id);

      /*normalize sample*/
      for (cs_lnum_t iel = 0; iel < n_elts_sample; iel++)
        sample[iel] = (sample[iel] - min_field) / (max_field - min_field);

      cs_real_t moment[]  = { 0.0, 0.0 };
      cs_real_t min_max[] = { 0.0, 0.0 };

      _compute_sample_moment(sample, weights, n_elts_sample, moment, min_max);
      /*update profile value*/
      profile->mean_f_n[l_id] = moment[0];
      profile->sd_f_n[l_id]   = moment[1];
    }
  }

  BFT_FREE(sample);
  BFT_FREE(weights);
}

/*----------------------------------------------------------------------------
 * Function to compute all created profiles
 *----------------------------------------------------------------------------*/

void
user_profiles_compute_all(void)
{
  cs_lnum_t n_profiles = _profile_list.n_profiles;

  for (int p_id = 0; p_id < n_profiles; p_id++) {
    user_profile_t *profile = _profile_list.profile_list[p_id];
    user_profile_compute(profile);
  }
}

/*----------------------------------------------------------------------------
 * Output current profile state (log file and csv). directory is created in
 * ./profiles if not already existing
 *
 * parameters:
 *   profile   <-- pointer to an initialized profile structure
 *   interval  <-- output time step interval
 *----------------------------------------------------------------------------*/

void
user_profile_output(user_profile_t  *profile,
                    int              interval)
{
  cs_time_step_t *ts = cs_get_glob_time_step();

  if (ts->nt_cur % interval != 0)
    return;

  /* Profile shorter variables */
  cs_lnum_t n_layers = profile->n_layers;

  if (cs_glob_rank_id < 1) {

    const char prefix[] = "./profiles";

    char dirname[80];
    char dirname_histograms[150];

    snprintf(dirname, 80, "./profiles/%s", profile->name);
    dirname[79] = '\0';

    /* Create main directory if not present */

    if (cs_file_isdir(prefix) != 1)
      cs_file_mkdir_default(prefix);

    if (cs_file_isdir(dirname) != 1) {
      cs_file_mkdir_default(dirname);

      /* Output STL layer mesh for possible visualisation of the layer */

      _user_output_profile_stl_mesh(profile, dirname);
      _output_profile_setup_log(profile, dirname);

      /* Output med mesh if generated */
#if defined(HAVE_MEDCOUPLING_LOADER)
      if (strcmp(profile->intersect_method, "MEDCOUPLING") == 0)
        _output_profile_med_mesh(profile, dirname);
#endif
    }

    /* Create histograms directory if not present */

    snprintf(dirname_histograms, 150, "%s/histograms", dirname);
    dirname_histograms[149] = '\0';

    if (cs_file_isdir(dirname_histograms) != 1)
      cs_file_mkdir_default(dirname_histograms);

    /* Output current state of profile */

    if (ts->nt_cur % interval == 0) {
      _output_profile_values_csv(profile, dirname);
      _output_profile_log(profile, dirname);
      for (int l_id = 0; l_id < n_layers; l_id++) {
        user_histogram_t *histogram = profile->histogram_list[l_id];
        user_output_histogram_csv(histogram, dirname_histograms);
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Function to output all created profiles
 *
 * parameters:
 *   interval <-- output time step interval
 *----------------------------------------------------------------------------*/

void
user_profiles_output_all(int  interval)
{
  cs_lnum_t n_profiles = _profile_list.n_profiles;

  for (int p_id = 0; p_id < n_profiles; p_id++) {
    user_profile_t *profile = _profile_list.profile_list[p_id];
    user_profile_output(profile, interval);
  }
}

/*----------------------------------------------------------------------------
 * Output current histogram in a csv file in the given directory. File created
 * if not already existing
 *
 * parameters:
 *   histogram    <-- pointer to a histogram strucutre
 *   dirname      <-- directory in which to output the profile
 *----------------------------------------------------------------------------*/

void
user_output_histogram_csv(user_histogram_t  *histogram,
                          const char        *dirname)
{
  /*shorter cs variables*/
  cs_time_step_t *ts = cs_get_glob_time_step();

  char filename[200];

  sprintf(filename, "%s/%s.csv", dirname, histogram->name);

  /*shorter Histogram variables*/
  cs_lnum_t  n_bins = histogram->n_bins;
  cs_real_t  min    = histogram->min;
  cs_real_t *h_i    = histogram->h_i;
  cs_real_t *l_i    = histogram->l_i;
  cs_real_t  mu     = histogram->mean;
  cs_real_t  sigma  = histogram->sd;
  cs_lnum_t  n_iter = histogram->n_iter;

  cs_lnum_t   n_l_header         = 0;
  const char *histogram_header[] = {
    "time_step", "time", "n_bins", "n_iter", "min", "mu",
    "sigma",     "l_0",  "h_0",       "...",    "l_n", "h_n",
  };

  n_l_header = sizeof(histogram_header) / sizeof(char *);

  FILE *ptrResults = NULL;

  /* Check if the file exist*/
  ptrResults = fopen(filename, "r+");

  if (ptrResults == NULL) {
    /* Create a file and write header*/
    ptrResults = fopen(filename, "w+");

    fprintf(ptrResults, "%s", histogram_header[0]);

    for (int h_id = 1; h_id < n_l_header; h_id++)
      fprintf(ptrResults, ",%s", histogram_header[h_id]);

    fprintf(ptrResults, "\n");
  }

  fclose(ptrResults);

  ptrResults = fopen(
    filename,
    "a+"); /* Add current result in the file re using same pattern than above*/

  fprintf(ptrResults,
          "%d,%f,%d,%d,%f,%f,%f",
          ts->nt_cur,
          ts->t_cur,
          n_bins,
          n_iter,
          min,
          mu,
          sigma);

  /* output width and height of each histogram class */
  for (int b_id = 0; b_id < n_bins; b_id++) {
    fprintf(ptrResults, ",%f,%f", l_i[b_id], h_i[b_id]);
  }

  fprintf(ptrResults, "\n");

  fclose(ptrResults);
}

/*----------------------------------------------------------------------------
 * Output histogram graphs using OpenTurns library (need to be linked to CS)
 *
 * parameters:
 *   histogram <-- pointer to the current histogram structure
 *   dirname   <-- pointer to directory output path
 *----------------------------------------------------------------------------*/

void
user_histogram_ot_output(user_histogram_t  *histogram,
                         const char        *dirname)
{
  /*wrap c++ function*/
  _output_histogram_ot(histogram, dirname);
}

/*----------------------------------------------------------------------------
 * Output layer histogram layer of a profile.
 *
 * Warning: + OpenTurns library needs to be linked
 *          + Graph creation is relatively slow, choose relevant interval
 *
 * parameters:
 *   profile   <-- pointer to an initialized profile structure
 *   interval  <-- output time step interval
 *----------------------------------------------------------------------------*/

void
user_profile_histogram_ot_output([[maybe_unused]] user_profile_t  *profile,
                                 [[maybe_unused]] int              interval)
{
#if HAVE_OT == 1
  /*Output the PDF historgam thks to OT lib features*/
  /*Assume histograms directory already created*/
  cs_time_step_t *ts = cs_get_glob_time_step();

  cs_lnum_t n_layers = profile->n_layers;

  char dirname[80];
  char dirname_histograms[150];

  sprintf(dirname, "./profiles/%s", profile->name);
  sprintf(dirname_histograms, "%s/histograms/%d", dirname, ts->nt_cur);

  if (cs_glob_rank_id <= 0 && (fabs(ts->nt_cur % interval) < 1.E-10)) {
    mkdir(dirname_histograms, S_IRWXU | S_IRGRP | S_IROTH);
    for (int l_id = 0; l_id < n_layers; l_id++) {
      user_histogram_t *histogram = profile->histogram_list[l_id];
      _output_histogram_ot(histogram, dirname_histograms);
    }
  }

#endif
}

/*----------------------------------------------------------------------------
 * Function to output all layer histogram layer of created profiles
 * Warning: OT lib is required
 *
 * parameters:
 *   interval    <-- output time step interval
 *----------------------------------------------------------------------------*/

void
user_profiles_histogram_ot_output_all(int  interval)
{
  cs_lnum_t n_profiles = _profile_list.n_profiles;

  for (int p_id = 0; p_id < n_profiles; p_id++) {
    user_profile_t *profile = _profile_list.profile_list[p_id];
    user_profile_histogram_ot_output(profile, interval);
  }
}

/*----------------------------------------------------------------------------
 * Free memory of all profiles created
 *----------------------------------------------------------------------------*/

void
user_free_profiles(void)
{
  for (int p_id = 0; p_id < _profile_list.n_profiles; p_id++) {
    user_profile_t *profile = _profile_list.profile_list[p_id];
    _free_profile_all(profile);
  }
  _profile_list.n_profiles = 0;
  BFT_FREE(_profile_list.profile_list);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
