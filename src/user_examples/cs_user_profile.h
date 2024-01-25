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

/*----------------------------------------------------------------------------*/

#ifndef CS_USER_PROFILE_H
#define CS_USER_PROFILE_H

#include "cs_headers.h"
#include "cs_stl.h"

/* Set this define to 0 if OpenTunrs not linked to code saturne
   set this define to 1 if OpenTurns linked to CS (histogramp output) */

#define HAVE_OT 0

BEGIN_C_DECLS

typedef struct _user_profile_med_t user_profile_med_t;

/*----------------------------------------------------------------------------
 * Layer histogram structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char       *name;       /* name of the histogram*/
  const char *field;      /* field name*/
  cs_lnum_t   n_bins;     /* histogramm class numbers */
  cs_lnum_t   n_bins_max; /* maximum class number*/
  cs_real_t   bandwidth;  /* histogram bandwidth*/
  cs_real_t  *h_i;        /* density of the class i */
  cs_real_t  *l_i;        /* with of the class i*/
  cs_real_t  *c_i;        /* center of the class i*/
  cs_real_t   mean;       /* mean of the sample (weighted or not)*/
  cs_real_t   sd;         /* standard deviation of the sample*/
  cs_real_t   min;        /* class min value*/
  cs_real_t   max;        /* class max value*/
  cs_real_t   Q1;         /* First quarter quantile*/
  cs_real_t   Q2;         /* Mediane*/
  cs_real_t   Q3;         /* Third quarter quantile*/
  cs_real_t   ptot;       /* sum of h_i : has to be 1.0*/
  cs_real_t   n_iter;     /* iteration to optimize bandWidth*/

} user_histogram_t;

/*----------------------------------------------------------------------------
 * Profile main structure
 *----------------------------------------------------------------------------*/

typedef struct {

  const char name[512];             // name of the profile
  const char field[512];            // field name
  const char criteria[512];         // cell selection criteria for profile
  const char intersect_method[512]; // Method used for intersaction : "Basic",
                                // "STL", "MEDCOUPLING"
  const char progression_law[512];  // progression law for layer thickness :
                                // CONSTANT, GEOMETRIC, PARABOLIC
  const char weighted[512]; // MASS, VOLUME, NO
  cs_real_t   dir_v[3]; // profile direction
  cs_real_t   min_dir;  // minimum distance along direction (from origin)
  cs_real_t   max_dir;  // maximum distance along direction (from origin)
  cs_real_t   i_v[3];   // assuming dir_v is the k vector, i_v j_v and k_v are an
                        // orthonormal direct base
  cs_real_t min_i;  // minimum distance along direction (from origin)
  cs_real_t max_i;  // maximu distance along direction (from origin)
  cs_real_t j_v[3]; // assuming dir_v is the k vector, i_v j_v and k_v are an
                    // orthonormal direct base
  cs_real_t  min_j; // minimum distance along direction (from origin)
  cs_real_t  max_j; // maximu distance along direction (from origin)
  cs_real_t *pos; // position of the center of each layer (n_v coord, n_v: dir_v
                  // normalized)
  cs_real_t *pos_n;       // position of the center of each layer (normalized)
  cs_lnum_t  n_layers;    // number of layer of the profile
  cs_real_t  progression; // progression parameter : only usefull for geometric
                         // and parabolic laws
  cs_real_t *l_thick; // array of thickness (one per layer)
  cs_lnum_t *n_cells; // number of cells per layers without stl immersion

  cs_real_t **cells_layer_vol; // [layer_id][cell_id] percent of vol in layer_id
  cs_real_t  *weigth;          // can be either cells mass or cells volume
  cs_real_t   sel_cells_weigth;  // weigth of selected cells (with criteria)
  cs_real_t   min_field;         // minimum field value in cells selection
  cs_real_t   max_field;         // maximum field value in cells selection
  cs_real_t  *mean_f;            // mean field value (norm if fdim >1)
  cs_real_t  *mean_f_n;          // normalized mean field value
  cs_real_t  *sd_f;              // field standard deviation
  cs_real_t  *sd_f_n;            // normalized field standard deviation
  cs_stl_mesh_t     **mesh_list; // Array of stl mesh (1 per layer)
  user_profile_med_t *med_mesh_struct; // pointer to structure housing med mesh
                                       // needed to interface cxx and c code
  user_histogram_t **histogram_list; // Array of histogram (1 per layer)

} user_profile_t;

typedef struct {

  user_profile_t **profile_list; // Array of user profiles
  cs_lnum_t        n_profiles;   // number of profiles

} user_profile_info_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return a pointer to a profile structure, null if the profile does not exist
 *
 * parameters:
 *   name <-- pointer to profile name
 *----------------------------------------------------------------------------*/

user_profile_t *
user_profile_get_by_name(const char  *name);

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
                    const char  *intersect_method);

/*----------------------------------------------------------------------------
 * Return an allocated histogram structure. Function is public for use histogram
 * oustide of mean profile context
 *
 * parameters :
 *   name       <-- name of the histogram
 *   field      <-- field on which the histogram is compute
 *   n_bin_max  <-- maximum number of bin of the bins
 *----------------------------------------------------------------------------*/

user_histogram_t *
user_create_histogram(char         *name,
                      const char  *field,
                      cs_lnum_t    n_bin_max);

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
                       cs_lnum_t          n_elts_sample);

/*----------------------------------------------------------------------------
 * Free histogram structure
 *
 * parameters :
 *   histogram    <-> name of the histogram
 *----------------------------------------------------------------------------*/

void
user_destroy_histogram(user_histogram_t  *histogram);

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
user_compute_cell_volume_per_layer(user_profile_t  *profile);

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
user_profile_compute(user_profile_t  *profile);

/*----------------------------------------------------------------------------
 * Function to compute all created profiles
 *----------------------------------------------------------------------------*/

void
user_profiles_compute_all(void);

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
                    int              interval);

/*----------------------------------------------------------------------------
 * Function to output all created profiles
 *
 * parameters:
 *   interval <-- output time step interval
 *----------------------------------------------------------------------------*/

void
user_profiles_output_all(int  interval);

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
                          const char        *dirname);

/*----------------------------------------------------------------------------
 * Output histogram graphs using OpenTurns library (need to be linked to CS)
 *
 * parameters:
 *   histogram <-- pointer to the current histogram structure
 *   dirname   <-- pointer to directory output path
 *----------------------------------------------------------------------------*/

void
user_histogram_ot_output(user_histogram_t  *histogram,
                         const char        *dirname);

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
user_profile_histogram_ot_output(user_profile_t  *profile,
                                 int              interval);

/*----------------------------------------------------------------------------
 * Function to output all layer histogram layer of created profiles
 * Warning: OT lib is required
 *
 * parameters:
 *   interval    <-- output time step interval
 *----------------------------------------------------------------------------*/

void
user_profiles_histogram_ot_output_all(int  interval);

/*----------------------------------------------------------------------------
 * Free memory of all profiles created
 *----------------------------------------------------------------------------*/

void
user_free_profiles(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* defined (CS_USER_PROFILE_H) */
