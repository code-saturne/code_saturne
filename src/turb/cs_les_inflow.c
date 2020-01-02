/*============================================================================
 * Turbulent inflow generation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include <stdio.h>
#include <stdlib.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_random.h"
#include "cs_timer.h"
#include "cs_mesh_location.h"
#include "cs_restart.h"
#include "cs_restart_default.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_les_inflow.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define EPZERO  1.E-12
#define RINFIN  1.E+30

#if !defined(HUGE_VAL)
#define HUGE_VAL  1.E+12
#endif

/*
 * Definition of some mathematical functions
 *   - module                 : || U ||
 *   - dot product            :   U.V
 *   - cross product          :  U x V
 *   - tensor trace           :  Tr(T)
 *   - (double) inner product :  U.T.V
 *
 *  T must be symmetric and order as : T11, T22, T33, T12, T13, T23
 */

#undef _MODULE_3D
#undef _DOT_PRODUCT_3D
#undef _CROSS_PRODUCT_3D
#undef _TENSOR_TRACE_3D
#undef _INNER_PRODUCT_3D

#define _MODULE_3D(module_u, u) \
  ( module_u = sqrt((u)[0]*(u)[0] + (u)[1]*(u)[1] + (u)[2]*(u)[2]) )

#define _DOT_PRODUCT_3D(dot_u_v, u, v) \
  ( dot_u_v = (u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2] )

#define _CROSS_PRODUCT_3D(cross_u_v, u, v) \
  ( (cross_u_v)[0] = (u)[1]*(v)[2] - (u)[2]*(v)[1], \
    (cross_u_v)[1] = (u)[2]*(v)[0] - (u)[0]*(v)[2], \
    (cross_u_v)[2] = (u)[0]*(v)[1] - (u)[1]*(v)[0]  )

#define _TENSOR_TRACE_3D(trace_t, t) \
  ( trace_t = (t)[0] + (t)[1] + (t)[2] )

#define _INNER_PRODUCT_3D(inner_u_t_v, u, t, v)        \
  ( inner_u_t_v = \
      (u)[0]*(t)[0]*(v)[0] + (u)[0]*(t)[3]*(v)[1] + (u)[0]*(t)[4]*(v)[2] \
    + (u)[1]*(t)[3]*(v)[0] + (u)[1]*(t)[1]*(v)[1] + (u)[1]*(t)[5]*(v)[2] \
    + (u)[2]*(t)[4]*(v)[0] + (u)[2]*(t)[5]*(v)[1] + (u)[2]*(t)[2]*(v)[2] )

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Inlet definition */
/*------------------*/

struct _cs_inlet_t {

  /* Synthetic inflow method */

  cs_inflow_type_t      type;
  void                 *inflow;

  int                   initialize;            /* Indicator of initialization */
  int                   verbosity;             /* Indicator of verbosity level*/

  /* Geometric informations */

  cs_int_t              n_faces;
  cs_int_t             *parent_num;

  cs_real_t            *face_centre;
  cs_real_t            *face_surface;

  /* Mean flow information */

  double                mean_velocity[3];      /* Mean velocity               */
  double                kinetic_energy;        /* Level of energy             */
  double                dissipation_rate;      /* Level of dissipation rate   */

  /* Counters */

  double                wt_tot;                /* Total wall-clock time used  */
  double                cpu_tot;               /* Total (local) CPU used      */

};

/* Batten method */
/*---------------*/

typedef struct _cs_inflow_batten_t {

  int            n_modes;                /* Number of modes                   */

  double        *frequency;              /* Frequency   : normal law N(1,1)   */
  double        *wave_vector;            /* Wave vector : normal law N(0,1/2) */
  double        *amplitude_cos;          /* Amplitudes of the cosines         */
  double        *amplitude_sin;          /* Amplitudes of the sines           */

} cs_inflow_batten_t;

/* Synthetic Eddy Method (SEM) */
/*-----------------------------*/

typedef struct _cs_inflow_sem_t {

  int            n_structures;        /* Number of coherent structures        */

  double        *position;            /* Position of the structures           */
  double        *energy;              /* Anisotropic energy of the structures */

} cs_inflow_sem_t;

typedef struct
{
  double  val;
  int     rank;
} _mpi_double_int_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Names for synthetic turbulence generation method */

const char *cs_inflow_type_name[] = {"Laminar",
                                     "Random",
                                     "Batten",
                                     "SEM"};

/* Structures associated to the inlets */

static int            cs_glob_inflow_n_inlets    = 0;
static cs_inlet_t   **cs_glob_inflow_inlet_array = NULL;
static cs_restart_t  *cs_glob_inflow_suite = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Generation of synthetic turbulence via a Gaussian random method.
 *
 * parameters:
 *   n_points          --> Local number of points where turbulence is generated
 *   fluctuations      <-- Velocity fluctuations generated
 *----------------------------------------------------------------------------*/

static void
_random_method(const cs_int_t   n_points,
               cs_real_t       *fluctuations)
{
  cs_int_t  point_id;

  int       coo_id;

  double    random[3];

  for (point_id = 0; point_id < n_points; point_id++) {
    cs_random_normal(3, random);
    for (coo_id = 0; coo_id < 3; coo_id++)
      fluctuations[point_id*3 + coo_id] = random[coo_id];
  }
}

/*----------------------------------------------------------------------------
 * Generation of synthetic turbulence via the Batten method.
 *
 * parameters:
 *   n_points          --> Local number of points where turbulence is generated
 *   point_coordinates --> Coordinates of the points
 *   initialize        --> Indicator of initialization
 *   inflow            --> Specific structure for Batten method
 *   time              --> Current time at the present iteration
 *   reynolds_stresses --> Reynolds stresses at each point
 *   dissipation_rate  --> Dissipation rate at each point
 *   fluctuations      <-- Velocity fluctuations generated at each point
 *----------------------------------------------------------------------------*/

static void
_batten_method(const cs_int_t       n_points,
               const cs_real_t     *point_coordinates,
               const int            initialize,
               cs_inflow_batten_t  *inflow,
               const cs_real_t      time,
               const cs_real_t     *reynolds_stresses,
               const cs_real_t     *dissipation_rate,
               cs_real_t           *fluctuations)
{
  cs_int_t  point_id;

  int       coo_id;
  int       mode_id;

  const double  two_pi              = 2.*acos(-1.);
  const double  sqrt_three_half     = sqrt(1.5);
  const double  sqrt_two_by_n_modes = sqrt(2./inflow->n_modes);

  const double *frequency     = inflow->frequency;
  const double *wave_vector   = inflow->wave_vector;
  const double *amplitude_cos = inflow->amplitude_cos;
  const double *amplitude_sin = inflow->amplitude_sin;

  if (initialize == 1) {

    const int     three_n_modes   = 3*inflow->n_modes;
    const double  one_by_sqrt_two = sqrt(0.5);

    if (cs_glob_rank_id <= 0) {

      /* Random generation of the n_modes frequencies following a normal
         law with a mean of 1 and a variance of 1 (i.e. N(1,1)). */

      cs_random_normal(inflow->n_modes, inflow->frequency);

      for (mode_id = 0; mode_id < inflow->n_modes; mode_id++)
        inflow->frequency[mode_id] += 1.;

      /* Random generation of the n_modes wave vectors following a normal
         law with a mean of 0 and a variance of 0.5 (i.e. N(0,1/2)). */

      cs_random_normal(three_n_modes, inflow->wave_vector);

      for (mode_id = 0; mode_id < inflow->n_modes; mode_id++)
        for (coo_id = 0; coo_id < 3; coo_id++)
          inflow->wave_vector[mode_id*3 + coo_id] *= one_by_sqrt_two;

      /* Generation of the n_modes amplitude vector for both the sines and
         the cosines. */

      for (mode_id = 0; mode_id < inflow->n_modes; mode_id++) {

        double  rcos[3];
        double  rsin[3];

        /* Temporary random vectors following a normal law N(0,1) necessary
           to compute the random amplitudes of the sines and cosines */

        cs_random_normal(3, rcos);
        cs_random_normal(3, rsin);

        _CROSS_PRODUCT_3D(inflow->amplitude_cos + mode_id*3,
                          rcos,
                          inflow->wave_vector   + mode_id*3);

        _CROSS_PRODUCT_3D(inflow->amplitude_sin + mode_id*3,
                          rsin,
                          inflow->wave_vector   + mode_id*3);

      }

    }

#if defined(HAVE_MPI)

    if (cs_glob_rank_id >= 0) {

      MPI_Bcast(inflow->frequency,   inflow->n_modes, CS_MPI_REAL, 0,
                cs_glob_mpi_comm);
      MPI_Bcast(inflow->wave_vector,   three_n_modes, CS_MPI_REAL, 0,
                cs_glob_mpi_comm);
      MPI_Bcast(inflow->amplitude_cos, three_n_modes, CS_MPI_REAL, 0,
                cs_glob_mpi_comm);
      MPI_Bcast(inflow->amplitude_sin, three_n_modes, CS_MPI_REAL, 0,
                cs_glob_mpi_comm);

    }

#endif

  }

  for (point_id = 0; point_id < n_points; point_id++) {

    double kinetic_energy;

    double time_scale;
    double velocity_scale;
    double lenght_scale;

    double spectral_time;
    double spectral_coordinates[3];

    /*
      Compute integral scales of turbulence :
      -  Tb = k / epsilon
      -  Vb = sqrt(k)
      -  Lb = Tb * Vb     ( = k^(3/2) / epsilon )
    */

    _TENSOR_TRACE_3D(kinetic_energy, reynolds_stresses + point_id*6);
    kinetic_energy *= 0.5;

    time_scale     = kinetic_energy / dissipation_rate[point_id];
    velocity_scale = sqrt(kinetic_energy);
    lenght_scale   = time_scale * velocity_scale;

    /* Spectral position of the point in space and time */

    spectral_time = two_pi * time / time_scale;

    for (coo_id = 0; coo_id < 3; coo_id++)
      spectral_coordinates[coo_id] =
        two_pi * point_coordinates[point_id*3 + coo_id] / lenght_scale;

    /* Compute the velocity fluctuations */

    for (mode_id = 0; mode_id < inflow->n_modes; mode_id++) {

      double spectral_velocity_scale = 0.;
      double norm_wave_vector        = 0.;
      double dxpot = 0.;
      double mod_wave_vector[3];

      _MODULE_3D(norm_wave_vector, wave_vector + mode_id*3);

      _INNER_PRODUCT_3D(spectral_velocity_scale,
                        wave_vector + mode_id*3,
                        reynolds_stresses + point_id*6,
                        wave_vector + mode_id*3);

      spectral_velocity_scale = sqrt_three_half*sqrt(spectral_velocity_scale)
        /norm_wave_vector;


      for (coo_id = 0; coo_id < 3; coo_id++)
        mod_wave_vector[coo_id] = wave_vector[mode_id*3 + coo_id]
          * velocity_scale / spectral_velocity_scale;

      _DOT_PRODUCT_3D(dxpot, mod_wave_vector, spectral_coordinates);

       dxpot += frequency[mode_id]*spectral_time;

      for (coo_id = 0; coo_id < 3; coo_id++)
        fluctuations[point_id*3 + coo_id] +=
          amplitude_cos[mode_id*3 + coo_id]*cos(dxpot)
          + amplitude_sin[mode_id*3 + coo_id]*sin(dxpot);

    }

    for (coo_id = 0; coo_id < 3; coo_id++)
      fluctuations[point_id*3 + coo_id] *= sqrt_two_by_n_modes;

  }
}

/*----------------------------------------------------------------------------
 * Generation of synthetic turbulence via the Synthetic Eddy Method (SEM).
 *
 * parameters:
 *   n_points          --> Local number of points where turbulence is generated
 *   num_face          --> Local id of inlet boundary faces
 *   point_coordinates --> Coordinates of the points
 *   point_ponderation --> Ponderation of the points (surface, volume or NULL)
 *   initialize        --> Indicator of initialization
 *   verbosity         --> Indicator of verbosity level
 *   inflow            --> Specific structure for Batten method
 *   time_step         --> Time step at the present iteration
 *   velocity          --> Velocity at each point
 *   reynolds_stresses --> Reynolds stresses at each point
 *   dissipation_rate  --> Dissipation rate at each point
 *   fluctuations      <-- Velocity fluctuations generated at each point
 *----------------------------------------------------------------------------*/

static void
_synthetic_eddy_method(const cs_int_t    n_points,
                       const cs_int_t   *num_face,
                       const cs_real_t  *point_coordinates,
                       const cs_real_t  *point_ponderation,
                       const int         initialize,
                       const int         verbosity,
                       cs_inflow_sem_t  *inflow,
                       const cs_real_t   time_step,
                       const cs_real_t  *velocity,
                       const cs_real_t  *reynolds_stresses,
                       const cs_real_t  *dissipation_rate,
                       cs_real_t        *fluctuations)
{
  cs_lnum_t   point_id;

  int        coo_id;
  int        struct_id;

  cs_lnum_t  j;

  double     alpha;
  double     random = -1.;

  double     pond_tot;

  double     average_velocity[3];

  double     box_volume;

  double     box_length[3];
  double     box_min_coord[3];
  double     box_max_coord[3];

  double    *length_scale;
  double    *ponderation;

  int compt[3] = {0, 0, 0};

  /* Computation of the characteristic scale of the synthetic eddies */
  /*-----------------------------------------------------------------*/

  BFT_MALLOC(length_scale, 3*n_points, double);

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *mesh_q = cs_glob_mesh_quantities;

  for (point_id = 0; point_id < n_points; point_id++) {

    cs_int_t  b_face_id = num_face[point_id] - 1;
    cs_lnum_t cell_id = mesh->b_face_cells[b_face_id];

    for (coo_id = 0; coo_id < 3; coo_id++) {

      double length_scale_min = -HUGE_VAL;

      for (j = mesh->b_face_vtx_idx[b_face_id];
           j < mesh->b_face_vtx_idx[b_face_id + 1]; j++) {
              cs_lnum_t vtx_id = mesh->b_face_vtx_lst[j];
              length_scale_min = CS_MAX(length_scale_min,
                              2.*CS_ABS(mesh_q->cell_cen[3*cell_id + coo_id]
                                        - mesh->vtx_coord[3*vtx_id + coo_id]));
      }

      length_scale[3*point_id + coo_id] =
              pow(1.5*reynolds_stresses[6*point_id + coo_id],1.5)
              /dissipation_rate[point_id];

      length_scale[3*point_id + coo_id] = 0.5*length_scale[3*point_id + coo_id];

      length_scale[3*point_id + coo_id] =
        CS_MAX(length_scale[3*point_id + coo_id],length_scale_min);

      if (CS_ABS(length_scale[3*point_id + coo_id]-length_scale_min) < EPZERO)
              compt[coo_id]++;

    }

  }

  if (verbosity > 0) {

    char     direction[3] = "xyz";
    cs_int_t sum[3] = {compt[0], compt[1], compt[2]};

    bft_printf(_("Max. size of synthetic eddies:\n"));

    for (coo_id = 0; coo_id < 3; coo_id++) {

      cs_real_t ls_max = -HUGE_VAL;
      cs_real_t xyzmax[3] = {0., 0., 0.};

      for (point_id = 0; point_id < n_points; point_id++) {

        ls_max = CS_MAX(length_scale[3*point_id + coo_id],ls_max);

        if (CS_ABS(ls_max - length_scale[3*point_id + coo_id]) < EPZERO)
          for (j = 0; j < 3; j++)
            xyzmax[j] = point_coordinates[3*point_id + j];

      }

#if defined(HAVE_MPI)

      if (cs_glob_rank_id >= 0) {

        _mpi_double_int_t  val_in, val_max;

        val_in.val  = ls_max;
        val_in.rank = cs_glob_rank_id;

        MPI_Allreduce(&val_in, &val_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
                      cs_glob_mpi_comm);

        ls_max = val_max.val;

        MPI_Bcast(xyzmax, 3, CS_MPI_REAL, val_max.rank, cs_glob_mpi_comm);

      }

#endif

      bft_printf(_("   max(sigma_%c) = %f, at coordinates (%f,%f,%f)\n"),
                 direction[coo_id], ls_max, xyzmax[0], xyzmax[1], xyzmax[2]);

    }
    bft_printf(_("\n"));

    bft_printf(_("Number of min. clippings (eddy size equals grid size):\n"));

#if defined(HAVE_MPI)

    if (cs_glob_rank_id >= 0)
      MPI_Allreduce(&compt, &sum, 3, CS_MPI_INT, MPI_SUM,
                    cs_glob_mpi_comm);

#endif

    for (coo_id = 0; coo_id < 3; coo_id++)
      bft_printf(_("   sigma_%c clipped %d times\n"),
                   direction[coo_id], sum[coo_id]);

    bft_printf(_("\n"));
  }

  /* Definition of the box on which eddies are generated */
  /*-----------------------------------------------------*/

  for (coo_id = 0; coo_id < 3; coo_id++) {
    box_min_coord[coo_id] =  HUGE_VAL;
    box_max_coord[coo_id] = -HUGE_VAL;
  }

  for (point_id = 0; point_id < n_points; point_id++)

    for (coo_id = 0; coo_id < 3; coo_id++) {

      box_min_coord[coo_id] =
        CS_MIN(box_min_coord[coo_id],
               point_coordinates[point_id*3 + coo_id]
               - length_scale[3*point_id + coo_id]);

      box_max_coord[coo_id] =
        CS_MAX(box_max_coord[coo_id],
               point_coordinates[point_id*3 + coo_id]
               + length_scale[3*point_id + coo_id]);

    }

#if defined(HAVE_MPI)

  if (cs_glob_rank_id >= 0) {

    double min_glob[3], max_glob[3];

    MPI_Allreduce(box_min_coord, &min_glob, 3, CS_MPI_REAL, MPI_MIN,
                  cs_glob_mpi_comm);

    MPI_Allreduce(box_max_coord, &max_glob, 3, CS_MPI_REAL, MPI_MAX,
                  cs_glob_mpi_comm);

    for (coo_id = 0; coo_id < 3; coo_id++) {
      box_min_coord[coo_id] = min_glob[coo_id];
      box_max_coord[coo_id] = max_glob[coo_id];
    }

  }

#endif

  for (coo_id = 0; coo_id < 3; coo_id++)
    box_length[coo_id] = box_max_coord[coo_id] - box_min_coord[coo_id];

  box_volume = box_length[0]*box_length[1]*box_length[2];

  if (verbosity > 0)
    bft_printf(_("Dimensions of the virtual box: \n"
                 "   Lx = %f, coo_min : %f, coo_max : %f\n"
                 "   Ly = %f, coo_min : %f, coo_max : %f\n"
                 "   Lz = %f, coo_min : %f, coo_max : %f\n\n"),
               box_length[0], box_min_coord[0], box_max_coord[0],
               box_length[1], box_min_coord[1], box_max_coord[1],
               box_length[2], box_min_coord[2], box_max_coord[2]);

  /* Initialization of the eddy field */
  /*----------------------------------*/

  if (initialize == 1) {

    if (cs_glob_rank_id <= 0) {

      for (struct_id = 0; struct_id < inflow->n_structures; struct_id++) {

        /* Random intensities */

        for (coo_id = 0; coo_id < 3; coo_id++) {

          cs_random_uniform(1, &random);
          inflow->energy[struct_id*3 + coo_id]
            = (random < 0.5) ? -1. : 1.;

        }

        /* Position of the eddies in the box */

        for (coo_id = 0; coo_id < 3; coo_id++) {

          cs_random_uniform(1, &random);
          inflow->position[struct_id*3 + coo_id] =
            box_min_coord[coo_id] + random*box_length[coo_id];

        }

      }

    }

#if defined(HAVE_MPI)

    if (cs_glob_rank_id >= 0) {

      MPI_Bcast(inflow->energy,   3*inflow->n_structures,
                CS_MPI_REAL, 0, cs_glob_mpi_comm);
      MPI_Bcast(inflow->position, 3*inflow->n_structures,
                CS_MPI_REAL, 0, cs_glob_mpi_comm);

    }

#endif

  }

  /* Estimation of the convection speed (with ponderation by surface) */
  /*------------------------------------------------------------------*/

  BFT_MALLOC(ponderation, n_points, double);

  if (point_ponderation == NULL)
    for (point_id = 0; point_id < n_points; point_id++)
      ponderation[point_id] = 1.;
  else
    for (point_id = 0; point_id < n_points; point_id++)
      ponderation[point_id] = point_ponderation[point_id];


  pond_tot = 0.;

  for (coo_id = 0; coo_id < 3; coo_id++)
    average_velocity[coo_id] = 0.;


  for (point_id = 0; point_id < n_points; point_id++) {

    pond_tot += ponderation[point_id];

    for (coo_id = 0; coo_id < 3; coo_id++)
      average_velocity[coo_id] += velocity[point_id*3 + coo_id]
        *ponderation[point_id];

  }

  BFT_FREE(ponderation);

  if (cs_glob_rank_id < 0)

    for (coo_id = 0; coo_id < 3; coo_id++)
      average_velocity[coo_id] /= pond_tot;

#if defined(HAVE_MPI)

  else {

    double pond_glob;
    double sum[3];

    MPI_Allreduce(&pond_tot, &pond_glob, 1, CS_MPI_REAL, MPI_SUM,
                  cs_glob_mpi_comm);

    MPI_Allreduce(average_velocity, &sum, 3, CS_MPI_REAL, MPI_SUM,
                  cs_glob_mpi_comm);

    for (coo_id = 0; coo_id < 3; coo_id++)
      average_velocity[coo_id] = sum[coo_id] / pond_glob;

  }

#endif

  /* Time evolution of the eddies */
  /*------------------------------*/

  if (cs_glob_rank_id <= 0) {

    /* Time advancement of the eddies */

    for (struct_id = 0; struct_id < inflow->n_structures; struct_id++) {

      for (coo_id = 0; coo_id < 3; coo_id++)
        inflow->position[struct_id*3 + coo_id] +=
          average_velocity[coo_id]*time_step;

    }

    /* Checking if the structures are still in the box */

    int compt_born = 0;

    for (struct_id = 0; struct_id < inflow->n_structures; struct_id++) {

      int new_struct = 0;
      int randomize[3] = {1, 1, 1};

      /* If the eddy leaves the box by one side, one convects it */

      for (coo_id = 0; coo_id < 3; coo_id++) {

        if (inflow->position[struct_id*3 + coo_id]
            < box_min_coord[coo_id]) {
          new_struct = 1;
          randomize[coo_id] = 0;
          inflow->position[struct_id*3 + coo_id] += box_length[coo_id];
        }
        else if (inflow->position[struct_id*3 + coo_id]
                 > box_max_coord[coo_id]) {
          new_struct = 1;
          randomize[coo_id] = 0;
          inflow->position[struct_id*3 + coo_id] -= box_length[coo_id];
        }

      }

      if (new_struct == 1) {

        /* The other directions are randomized */

        for (coo_id = 0; coo_id < 3; coo_id++) {

          if (randomize[coo_id] == 1) {
            cs_random_uniform(1, &random);
            inflow->position[struct_id*3 + coo_id]
              = box_min_coord[coo_id] + random*box_length[coo_id];
          }

        }

        /* New randomization of the energy */

        for (coo_id = 0; coo_id < 3; coo_id++) {

          cs_random_uniform(1, &random);
          inflow->energy[struct_id*3 + coo_id]
            = (random < 0.5) ? -1. : 1.;
        }

      }

      compt_born += new_struct;

    }

    if (verbosity > 0)
      bft_printf(_("Number of eddies leaving the box (regenerated): %i\n\n"),
                 compt_born);

  }

#if defined(HAVE_MPI)

  if (cs_glob_rank_id >= 0) {

    MPI_Bcast(inflow->energy,   3*inflow->n_structures,
              CS_MPI_REAL, 0, cs_glob_mpi_comm);
    MPI_Bcast(inflow->position, 3*inflow->n_structures,
              CS_MPI_REAL, 0, cs_glob_mpi_comm);

  }

#endif

  /* Computation of the eddy signal */
  /*--------------------------------*/

  alpha = sqrt(box_volume / (double) inflow->n_structures);

  for (point_id = 0; point_id < n_points; point_id++) {

    double distance[3];

    for (struct_id = 0; struct_id < inflow->n_structures; struct_id++) {

      for (coo_id = 0; coo_id < 3; coo_id++)
        distance[coo_id] =
          CS_ABS(point_coordinates[point_id*3 + coo_id]
                 - inflow->position[struct_id*3 + coo_id]);

      if (distance[0] < length_scale[point_id*3 + 0] &&
          distance[1] < length_scale[point_id*3 + 1] &&
          distance[2] < length_scale[point_id*3 + 2]) {

        double form_function = 1.;
        for (coo_id = 0; coo_id < 3; coo_id++)
          form_function *=
            (1.-distance[coo_id]/length_scale[point_id*3 + coo_id])
            /sqrt(2./3.*length_scale[point_id*3 + coo_id]);

        for (coo_id = 0; coo_id < 3; coo_id++)
          fluctuations[point_id*3 + coo_id] +=
            inflow->energy[struct_id*3 + coo_id]*form_function;

      }

    }

    for (coo_id = 0; coo_id < 3; coo_id++)
      fluctuations[point_id*3 + coo_id] *= alpha;

  }

  BFT_FREE(length_scale);
}

/*----------------------------------------------------------------------------
 * Rescaling of the fluctuations by the statistics following the Lund method.
 *
 * One assumes that the statistics are interlaced and ordered as follow :
 *   <u'u'>  <v'v'>  <w'w'>  <u'v'>  <u'w'>  <v'w'>
 *
 * parameters:
 *   n_points          --> Local number of points where turbulence is generated
 *   statistics        --> Statistics (i.e. Reynolds stresses)
 *   fluctuations      <-- Velocity fluctuations generated
 *----------------------------------------------------------------------------*/

static void
_rescale_fluctuations(cs_int_t     n_points,
                      cs_real_t   *statistics,
                      cs_real_t   *fluctuations)
{
  cs_int_t   point_id;

  for (point_id = 0; point_id < n_points; point_id++) {

    /* Reynolds stresses */

    double R11 = statistics[point_id*6 + 0];
    double R22 = statistics[point_id*6 + 1];
    double R33 = statistics[point_id*6 + 2];
    double R12 = statistics[point_id*6 + 3];
    double R13 = statistics[point_id*6 + 4];
    double R23 = statistics[point_id*6 + 5];

    /* Lund's coefficients */

    double a11 = sqrt(R11);
    double a21 = R12 / a11;
    double a22 = sqrt(R22 - a21*a21);
    double a31 = R13 / a11;
    double a32 = (R23 - a21*a31) / a22;
    double a33 = sqrt(R33 - a31*a31 - a32*a32);

    /* Rescaling of velocity fluctuations */

    double up_corr = a11*fluctuations[point_id*3];
    double vp_corr = a21*fluctuations[point_id*3]
      + a22*fluctuations[point_id*3 + 1];
    double wp_corr = a31*fluctuations[point_id*3]
      + a32*fluctuations[point_id*3 + 1]
      + a33*fluctuations[point_id*3 + 2];

    fluctuations[point_id*3    ] = up_corr;
    fluctuations[point_id*3 + 1] = vp_corr;
    fluctuations[point_id*3 + 2] = wp_corr;

  }

}

/*----------------------------------------------------------------------------
 * Modify the normal component of the fluctuations such that the mass flowrate
 * of the fluctuating field is zero.
 *
 * parameters:
 *   n_points          --> Local number of points where turbulence is generated
 *   num_face          --> Local id of inlet boundary faces
 *   fluctuations      <-> Velocity fluctuations
 *----------------------------------------------------------------------------*/

static void
_rescale_flowrate(cs_int_t     n_points,
                  cs_int_t    *num_face,
                  cs_real_t   *fluctuations)
{
  /* Compute the mass flow rate of the fluctuating field */
  /* and the area of the inlet */

  cs_lnum_t point_id;

  double mass_flow_rate = 0., mass_flow_rate_g = 0.;
  double area = 0., area_g = 0.;
  cs_real_t *density = CS_F_(rho)->val;
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *mesh_q = cs_glob_mesh_quantities;

  for (point_id = 0; point_id < n_points; point_id++) {

    double dot_product = 0.;

    cs_lnum_t b_face_id = num_face[point_id] - 1;
    cs_lnum_t cell_id = mesh->b_face_cells[b_face_id];

    double fluct[3] = { fluctuations[point_id*3],
                        fluctuations[point_id*3 + 1],
                        fluctuations[point_id*3 + 2] };
    double normal[3] = { mesh_q->b_face_normal[b_face_id*3],
                         mesh_q->b_face_normal[b_face_id*3 + 1],
                         mesh_q->b_face_normal[b_face_id*3 + 2] };

    _DOT_PRODUCT_3D(dot_product, fluct, normal);

    mass_flow_rate += density[cell_id]*dot_product;
    area = area + mesh_q->b_face_surf[b_face_id];

  }
  mass_flow_rate_g = mass_flow_rate;
  area_g = area;

#if defined(HAVE_MPI)
  if (cs_glob_rank_id >= 0) {
    MPI_Allreduce(&mass_flow_rate, &mass_flow_rate_g, 1, CS_MPI_REAL, MPI_SUM,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&area, &area_g, 1, CS_MPI_REAL, MPI_SUM, cs_glob_mpi_comm);
  }
#endif

  for (point_id = 0; point_id < n_points; point_id++) {

    /* Decompose the fluctuation in a local coordinate system */
    /* (not valid for warped boundary faces) */

    int coo_id;
    cs_lnum_t b_face_id = num_face[point_id] - 1;
    cs_lnum_t cell_id = mesh->b_face_cells[b_face_id];

    cs_lnum_t idx = mesh->b_face_vtx_idx[b_face_id];
    cs_lnum_t vtx_id1 = mesh->b_face_vtx_lst[idx];
    cs_lnum_t vtx_id2 = mesh->b_face_vtx_lst[idx+1];

    double norm = 0.;
    double normal_comp = 0., tangent_comp1 = 0., tangent_comp2 = 0.;
    double normal_unit[3], tangent_unit1[3], tangent_unit2[3];

    double fluct[3] = { fluctuations[point_id*3],
                        fluctuations[point_id*3 + 1],
                        fluctuations[point_id*3 + 2] };

    for (coo_id = 0; coo_id < 3; coo_id++) {
      normal_unit[coo_id] = mesh_q->b_face_normal[b_face_id*3 + coo_id];
      normal_unit[coo_id] /= mesh_q->b_face_surf[b_face_id];
    }

    for (coo_id = 0; coo_id < 3; coo_id++) {
      tangent_unit1[coo_id] = mesh->vtx_coord[3*vtx_id1 + coo_id]
                            - mesh->vtx_coord[3*vtx_id2 + coo_id];
    }

    _CROSS_PRODUCT_3D(tangent_unit2, normal_unit, tangent_unit1);

    _MODULE_3D(norm, tangent_unit1);
    for (coo_id = 0; coo_id < 3; coo_id++)
      tangent_unit1[coo_id] /= norm;

    _MODULE_3D(norm, tangent_unit2);
    for (coo_id = 0; coo_id < 3; coo_id++)
      tangent_unit2[coo_id] /= norm;

    _DOT_PRODUCT_3D(normal_comp, fluct, normal_unit);
    _DOT_PRODUCT_3D(tangent_comp1, fluct, tangent_unit1);
    _DOT_PRODUCT_3D(tangent_comp2, fluct, tangent_unit2);

    /* Rescale the normal component and return in cartesian coordinates*/

    normal_comp -= mass_flow_rate_g/(density[cell_id]*area_g);

    for (coo_id = 0; coo_id < 3; coo_id++)
      fluctuations[point_id*3 + coo_id] = normal_comp*normal_unit[coo_id]
                                        + tangent_comp1*tangent_unit1[coo_id]
                                        + tangent_comp2*tangent_unit2[coo_id];

  }
}

/*----------------------------------------------------------------------------
 * Add a new inlet for synthetic turbulence inflow generation
 *----------------------------------------------------------------------------*/

static void
_cs_inflow_add_inlet(cs_inflow_type_t   type,
                     const cs_int_t     n_faces,
                     const cs_int_t    *num_face,
                     const int          n_entities,
                     const int          verbosity,
                     const cs_real_t   *mean_velocity,
                     const cs_real_t    kinetic_energy,
                     const cs_real_t    dissipation_rate)
{
  cs_int_t  face_id;

  int       coo_id;

  cs_inlet_t   *inlet = NULL;

  const cs_mesh_quantities_t  *mesh_q = cs_glob_mesh_quantities;

  /* Allocating inlet structures */
  /*-----------------------------*/

  BFT_REALLOC(cs_glob_inflow_inlet_array,
              cs_glob_inflow_n_inlets + 1, cs_inlet_t *);

  BFT_MALLOC(inlet, 1, cs_inlet_t);

  /* Mesh */

  inlet->n_faces = n_faces;

  inlet->parent_num  = NULL;
  inlet->face_centre = NULL;
  inlet->face_surface = NULL;

  if (inlet->n_faces > 0) {

    BFT_MALLOC(inlet->parent_num, inlet->n_faces, cs_int_t);
    for (face_id = 0; face_id < n_faces; face_id++)
      inlet->parent_num[face_id] = num_face[face_id];

    BFT_MALLOC(inlet->face_centre, 3*inlet->n_faces, cs_real_t);
    for (face_id = 0; face_id < inlet->n_faces; face_id++)
      for (coo_id = 0; coo_id < 3; coo_id++)
        inlet->face_centre[face_id*3 + coo_id]
          = mesh_q->b_face_cog[(num_face[face_id]-1)*3 + coo_id];

    BFT_MALLOC(inlet->face_surface, inlet->n_faces, cs_real_t);
    for (face_id = 0; face_id < inlet->n_faces; face_id++)
      inlet->face_surface[face_id]
        = sqrt(  pow(mesh_q->b_face_normal[(num_face[face_id]-1)*3 + 0],2)
               + pow(mesh_q->b_face_normal[(num_face[face_id]-1)*3 + 1],2)
               + pow(mesh_q->b_face_normal[(num_face[face_id]-1)*3 + 2],2));

  }

  /* Turbulence level */

  for (coo_id = 0; coo_id < 3; coo_id++)
    inlet->mean_velocity[coo_id] = mean_velocity[coo_id];

  inlet->kinetic_energy   = kinetic_energy;
  inlet->dissipation_rate = dissipation_rate;

  /* Generation method of synthetic turbulence */
  /*-------------------------------------------*/

  if (type > 3)
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid choice of synthetic turbulence generation method (%d).\n"
                "Valid choices are:\n"
                "\t0 -> laminar\n\t1 -> random\n\t2 -> batten\n\t3 -> SEM\n"),
              type);
  else
    inlet->type = type;

  switch(inlet->type) {

  case CS_INFLOW_LAMINAR:

    inlet->inflow = NULL;
    break;

  case CS_INFLOW_RANDOM:

    inlet->inflow = NULL;
    break;

  case CS_INFLOW_BATTEN:

    {
      cs_inflow_batten_t *inflow;

      if (n_entities <= 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("The number of modes for the Batten method "
                    "must be strictly positive. %d is given here.\n"),
                  n_entities);

      BFT_MALLOC(inflow, 1, cs_inflow_batten_t);

      inflow->n_modes = n_entities;

      BFT_MALLOC(inflow->frequency,       inflow->n_modes, double);
      BFT_MALLOC(inflow->wave_vector,   3*inflow->n_modes, double);
      BFT_MALLOC(inflow->amplitude_cos, 3*inflow->n_modes, double);
      BFT_MALLOC(inflow->amplitude_sin, 3*inflow->n_modes, double);

      inlet->inflow = inflow;
    }
    break;

  case CS_INFLOW_SEM:

    {
      cs_inflow_sem_t *inflow;

      if (n_entities <= 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("The number of eddies for the SEM "
                    "must be strictly positive. %d is given here.\n"),
                  n_entities);

      BFT_MALLOC(inflow, 1, cs_inflow_sem_t);

      inflow->n_structures = n_entities;

      BFT_MALLOC(inflow->position, 3*inflow->n_structures, double);
      BFT_MALLOC(inflow->energy,   3*inflow->n_structures, double);

      inlet->inflow = inflow;
    }
    break;

  }

  /* Others */
  /*--------*/

  inlet->initialize = 1;

  inlet->verbosity = verbosity;

  inlet->wt_tot  = 0.;
  inlet->cpu_tot = 0.;

  /* Global array of inlets */
  /*------------------------*/

  cs_glob_inflow_inlet_array[cs_glob_inflow_n_inlets] = inlet;
  cs_glob_inflow_n_inlets++;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a structure for the inlets
 *----------------------------------------------------------------------------*/

void CS_PROCF(defsyn, DEFSYN)
(
 cs_int_t            *n_inlets     /* --> number of inlets                    */
)
{
  cs_int_t inlet_id;

  /* Definition of the global parameters of the inlets */

  CS_PROCF(cs_user_les_inflow_init, CS_USER_LES_INFLOW_INIT)(n_inlets);

  for (inlet_id = 0; inlet_id < *n_inlets; inlet_id++) {

    const cs_mesh_t  *mesh = cs_glob_mesh;
    const cs_int_t    nument = inlet_id + 1;
    int      j;
    cs_int_t n_faces_g;

    /* Initializations */

    int        type = 0;
    cs_int_t   n_faces = 0;
    int        n_entities = 0;
    int        verbosity = 0;

    cs_int_t  *index_face = NULL;

    cs_real_t  mean_velocity[3] = {0., 0., 0.};

    cs_real_t  kinetic_energy = 0.;
    cs_real_t  dissipation_rate = 0.;

    BFT_MALLOC(index_face, mesh->n_b_faces, cs_int_t);

    for (j = 0; j < mesh->n_b_faces; j++)
      index_face[j] = 0;

    bft_printf(_(" Definition of the LES inflow boundary \"%d\" \n"),
               cs_glob_inflow_n_inlets + 1);

    /* User definition of the parameters of the inlet */

    CS_PROCF(cs_user_les_inflow_define, CS_USER_LES_INFLOW_DEFINE)(
                             &nument,
                             &type,
                             &n_entities,
                             &verbosity,
                             &n_faces,
                             index_face,
                             mean_velocity,
                             &kinetic_energy,
                             &dissipation_rate);

    n_faces_g = n_faces;

#if defined(HAVE_MPI)

    if (cs_glob_rank_id >= 0)
      MPI_Allreduce(&n_faces, &n_faces_g, 1, CS_MPI_INT, MPI_SUM,
                    cs_glob_mpi_comm);

#endif

    if (n_faces_g == 0)
      bft_error(__FILE__, __LINE__, 0,
               _("Abort while defing the LES inlets.\n"
                 "The LES inlet \"%d\" does not contain any boundary face.\n"
                 "Verify the definition of the LES inlets "
                 "(cs_user_les_inflow.f90 file).\n"),nument);

    /* Add an inlet to the global inlets array */

    _cs_inflow_add_inlet((cs_inflow_type_t) type,
                         n_faces,
                         index_face,
                         n_entities,
                         verbosity,
                         mean_velocity,
                         kinetic_energy,
                         dissipation_rate);

    BFT_FREE(index_face);

    bft_printf(_("   Method: %d (%s)\n"
                 "   Number of boundary faces (global): %d\n"),
               type,cs_inflow_type_name[type],n_faces_g);

    if (type == 2)
      bft_printf(_("   Number of modes: %d\n\n"),n_entities);
    else if (type == 3)
      bft_printf(_("   Number of structures: %d\n\n"),n_entities);
    else
      bft_printf(_("   \n"));

  }

  bft_printf(" ------------------------------------------------------------- \n"
             "\n");

}

/*----------------------------------------------------------------------------
 * General synthetic turbulence generation
 *----------------------------------------------------------------------------*/

void CS_PROCF(synthe, SYNTHE)
(
 const cs_int_t  *const nvar,      /* --> number of variables                 */
 const cs_int_t  *const nscal,     /* --> number of scalars                   */
 const cs_int_t  *const iu,        /* --> index of velocity component         */
 const cs_int_t  *const iv,        /* --> index of velocity component         */
 const cs_int_t  *const iw,        /* --> index of velocity component         */
 const cs_real_t *const ttcabs,    /* --> current physical time               */
 const cs_real_t        dt[],      /* --> time step                           */
       cs_real_t        rcodcl[]   /* <-> boundary conditions array           */
)
{
  cs_int_t  face_id;

  int       inlet_id;
  int       coo_id;

  const double two_third = 2./3.;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  const cs_lnum_t  n_b_faces = mesh->n_b_faces;

  if (cs_glob_inflow_n_inlets == 0)
    return;

  for (inlet_id = 0; inlet_id < cs_glob_inflow_n_inlets; inlet_id++) {

    int nument = inlet_id + 1;

    cs_inlet_t *inlet = cs_glob_inflow_inlet_array[inlet_id];

    cs_real_t   *mean_velocity = NULL;
    cs_real_t   *reynolds_stresses = NULL;
    cs_real_t   *dissipation_rate = NULL;

    cs_real_t   *fluctuations = NULL;

    double wt_start, wt_stop, cpu_start, cpu_stop;

    wt_start  = cs_timer_wtime();
    cpu_start = cs_timer_cpu_time();

    /* Mean velocity profile, one-point statistics and dissipation rate */
    /*------------------------------------------------------------------*/

    BFT_MALLOC(mean_velocity,     3*inlet->n_faces, cs_real_t);
    BFT_MALLOC(reynolds_stresses, 6*inlet->n_faces, cs_real_t);
    BFT_MALLOC(dissipation_rate,    inlet->n_faces, cs_real_t);

    /* Initialization by the turbulence scales given by the user */

    for (face_id = 0; face_id < inlet->n_faces; face_id++) {

      for (coo_id = 0; coo_id < 3; coo_id++)
        mean_velocity[face_id*3 + coo_id] = inlet->mean_velocity[coo_id];

      for (coo_id = 0; coo_id < 3; coo_id++)
        reynolds_stresses[face_id*6 + coo_id] = two_third*inlet->kinetic_energy;

      for (coo_id = 0; coo_id < 3; coo_id++)
        reynolds_stresses[face_id*6 + 3 + coo_id] = 0.;

      dissipation_rate[face_id] = inlet->dissipation_rate;

    }

    /* Modification by the user */

    CS_PROCF(cs_user_les_inflow_advanced, CS_USER_LES_INFLOW_ADVANCED)
      (&nument, &inlet->n_faces,
       nvar, nscal,
       inlet->parent_num,
       dt,
       mean_velocity, reynolds_stresses, dissipation_rate);

    /* Generation of the synthetic turbulence */
    /*----------------------------------------*/

    BFT_MALLOC(fluctuations, 3*inlet->n_faces, cs_real_t);

    for (face_id = 0; face_id < inlet->n_faces; face_id++)
      for (coo_id = 0; coo_id < 3; coo_id++)
        fluctuations[face_id*3 + coo_id] = 0.;

    switch(inlet->type) {

    case CS_INFLOW_LAMINAR:
      break;
    case CS_INFLOW_RANDOM:
      _random_method(inlet->n_faces,
                     fluctuations);
      break;
    case CS_INFLOW_BATTEN:
      _batten_method(inlet->n_faces,
                     inlet->face_centre,
                     inlet->initialize,
                     (cs_inflow_batten_t *) inlet->inflow,
                     *ttcabs,
                     reynolds_stresses,
                     dissipation_rate,
                     fluctuations);
      break;
    case CS_INFLOW_SEM:
      {
        if (inlet->verbosity > 0)
          bft_printf(_("\n------------------------------"
                       "-------------------------------\n\n"
                       "SEM INFO, inlet \"%d\" \n\n"), nument);

        _synthetic_eddy_method(inlet->n_faces,
                               inlet->parent_num,
                               inlet->face_centre,
                               inlet->face_surface,
                               inlet->initialize,
                               inlet->verbosity,
                               (cs_inflow_sem_t *) inlet->inflow,
                               dt[0],
                               mean_velocity,
                               reynolds_stresses,
                               dissipation_rate,
                               fluctuations);

        if (inlet->verbosity > 0)
          bft_printf("------------------------------"
                     "-------------------------------\n");
      }
      break;
    }

    inlet->initialize = 0;

    BFT_FREE(dissipation_rate);

    /* Rescaling of the synthetic fluctuations by the statistics */
    /*-----------------------------------------------------------*/

    if (inlet->type == CS_INFLOW_RANDOM || inlet->type == CS_INFLOW_BATTEN
        || inlet->type == CS_INFLOW_SEM)
      _rescale_fluctuations(inlet->n_faces,
                            reynolds_stresses,
                            fluctuations);

    BFT_FREE(reynolds_stresses);

    /* Rescaling of the mass flow rate */
    /*---------------------------------*/

    if (inlet->type == CS_INFLOW_RANDOM || inlet->type == CS_INFLOW_BATTEN
        || inlet->type == CS_INFLOW_SEM)
      _rescale_flowrate(inlet->n_faces,
                        inlet->parent_num,
                        fluctuations);

    /* Boundary conditions */
    /*---------------------*/

    for (face_id = 0; face_id < inlet->n_faces; face_id++) {

      cs_int_t index_face = inlet->parent_num[face_id]-1;

      rcodcl[0*n_b_faces*(*nvar) + (*iu - 1)*n_b_faces + index_face] =
        mean_velocity[face_id*3 + 0]
        + fluctuations[face_id*3 + 0];

      rcodcl[0*n_b_faces*(*nvar) + (*iv - 1)*n_b_faces + index_face] =
        mean_velocity[face_id*3 + 1]
        + fluctuations[face_id*3 + 1];

      rcodcl[0*n_b_faces*(*nvar) + (*iw - 1)*n_b_faces + index_face] =
        mean_velocity[face_id*3 + 2]
        + fluctuations[face_id*3 + 2];

    }

    BFT_FREE(mean_velocity);
    BFT_FREE(fluctuations);

    wt_stop  = cs_timer_wtime();
    cpu_stop = cs_timer_cpu_time();

    inlet->wt_tot  += (wt_stop - wt_start);
    inlet->cpu_tot += (cpu_stop - cpu_start);
  }

}

/*----------------------------------------------------------------------------
 * Read the restart file of the LES inflow module
 *
 * Fortran interface:
 *
 * SUBROUTINE LECSYN
 * *****************
 *
 * character(kind=c_char)  filename : <-- : Name of the restart file
 *----------------------------------------------------------------------------*/

void CS_PROCF(lecsyn, LECSYN)
(
 const char  *filename
)
{
  bool                corresp_cel, corresp_fac, corresp_fbr, corresp_som;
  cs_int_t            nbvent;
  cs_int_t            indfac, ierror;

  cs_restart_t             *suite;
  cs_mesh_location_type_t   support;
  cs_restart_val_type_t     typ_val;

  bft_printf(_(" Reading the LES inflow module restart file...\n"));

  ierror = CS_RESTART_SUCCESS;

  /* Open the restart file */
  cs_glob_inflow_suite
    = cs_restart_create(filename, NULL, CS_RESTART_MODE_READ);

  if (cs_glob_inflow_suite == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while opening the LES inflow module restart file "
                "in read mode.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              filename);

  /* Pointer to the global restart structure */
  suite = cs_glob_inflow_suite;

  /* Verification of the associated "support" to the restart file */
  cs_restart_check_base_location(suite, &corresp_cel, &corresp_fac,
                                 &corresp_fbr, &corresp_som);

  /* Only boundary faces are of interest */
  indfac = (corresp_fbr == true ? 1 : 0);
  if (indfac == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while reading the LES inflow module restart file.\n"
                "The number of boundary faces has been modified\n"
                "Verify that the restart file corresponds to "
                "the present study.\n"));


  { /* Read the header */
    char       nomrub[] = "version_fichier_suite_turbulence_synthetique";
    cs_int_t   *tabvar;

    BFT_MALLOC(tabvar, 1, cs_int_t);

    nbvent  = 1;
    support = CS_MESH_LOCATION_NONE;
    typ_val = CS_TYPE_cs_int_t;

    ierror = cs_restart_read_section(suite,
                                     nomrub,
                                     support,
                                     nbvent,
                                     typ_val,
                                     tabvar);

    if (ierror < CS_RESTART_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("Abort while reading the LES inflow module restart file.\n"
                  "\n"
                  "The file %s does not seem to be a restart file\n"
                  "for the LES inflow module.\n"
                  "The calculation will not be run.\n"
                  "\n"
                  "Verify that the restart file corresponds to a\n"
                  "restart file for the LES inflow module.\n"),
                filename);

    BFT_FREE(tabvar);
  }

  { /* Read the number of inlets */
    char       nomrub[] = "nb_inlets";
    cs_int_t   *tabvar;
    cs_int_t   n_inlets;

    BFT_MALLOC(tabvar, 1, cs_int_t);

    nbvent  = 1;
    support = CS_MESH_LOCATION_NONE;
    typ_val = CS_TYPE_cs_int_t;

    ierror = cs_restart_read_section(suite,
                                     nomrub,
                                     support,
                                     nbvent,
                                     typ_val,
                                     tabvar);

    if (ierror < CS_RESTART_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("Problem while reading section in the restart file\n"
                  "for the LES inflow module:\n"
                  "<%s>\n"
                  "The calculation will not be run.\n"), nomrub);

    n_inlets = *tabvar;

    /* Coherence check */
    if (cs_glob_inflow_n_inlets != n_inlets)
      bft_error(__FILE__, __LINE__, 0,
                _("Stop reading the LES inflow module restart file.\n"
                  "The calculation is defined with %d LES inlets "
                  "while the restart file contains %d.\n"),
                cs_glob_inflow_n_inlets, n_inlets);

    BFT_FREE(tabvar);
  }

  { /* Read the structure of each inlet */
    int inlet_id;
    cs_int_t   *tabvar;

    for (inlet_id = 0; inlet_id < cs_glob_inflow_n_inlets; inlet_id++) {

      cs_inlet_t *inlet = cs_glob_inflow_inlet_array[inlet_id];

      { /* type of inlet */
        char             nomrub[] = "type_inlet";
        cs_inflow_type_t type;

        BFT_MALLOC(tabvar, 1, cs_int_t);

        nbvent  = 1;
        support = CS_MESH_LOCATION_NONE;
        typ_val = CS_TYPE_cs_int_t;

        ierror = cs_restart_read_section(suite,
                                         nomrub,
                                         support,
                                         nbvent,
                                         typ_val,
                                         tabvar);

        if (ierror < CS_RESTART_SUCCESS)
          bft_error(__FILE__, __LINE__, 0,
                    _("Problem while reading section in the restart file\n"
                      "for the LES inflow module:\n"
                      "<%s>\n"
                      "The calculation will not be run.\n"), nomrub);

        type = (cs_inflow_type_t) *tabvar;

        /* Coherence check */
        if (inlet->type != type)
          bft_error(__FILE__, __LINE__, 0,
                    _("Stop reading the LES inflow module restart file.\n"
                      "The inlet %d uses the method %d (%s) instead of "
                      "%d (%s) in the restart file.\n"),
                    inlet_id + 1,
                    inlet->type, cs_inflow_type_name[inlet->type],
                    type, cs_inflow_type_name[type]);

        BFT_FREE(tabvar);
      }

      switch(inlet->type) {

      case CS_INFLOW_LAMINAR:
        break;

      case CS_INFLOW_RANDOM:
        break;

      case CS_INFLOW_BATTEN:

        {
          cs_inflow_batten_t *inflow = (cs_inflow_batten_t *) inlet->inflow;

          { /* number of modes */
            char nomrub[] = "batten_number_modes";
            cs_int_t n_modes;

            BFT_MALLOC(tabvar, 1, cs_int_t);

            nbvent  = 1;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_int_t;

            ierror = cs_restart_read_section(suite,
                                             nomrub,
                                             support,
                                             nbvent,
                                             typ_val,
                                             tabvar);

            if (ierror < CS_RESTART_SUCCESS)
              bft_error(__FILE__, __LINE__, 0,
                        _("Problem while reading section in the restart file\n"
                          "for the LES inflow module:\n"
                          "<%s>\n"
                          "The calculation will not be run.\n"), nomrub);

            n_modes = *tabvar;

            /* Coherence check */
            if (inflow->n_modes != n_modes)
              bft_error(__FILE__, __LINE__, 0,
                        _("Stop reading the LES inflow module restart file.\n"
                          "%d modes are given for the Batten method "
                          "while the restart file contains %d.\n"),
                        inflow->n_modes, n_modes);

            BFT_FREE(tabvar);
          }

          { /* frequencies */
            char nomrub1[] = "batten_frequencies";

            nbvent  = inflow->n_modes;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            ierror = cs_restart_read_section(suite,
                                             nomrub1,
                                             support,
                                             nbvent,
                                             typ_val,
                                             inflow->frequency);

            if (ierror < CS_RESTART_SUCCESS)
              bft_error(__FILE__, __LINE__, 0,
                        _("Problem while reading section in the restart file\n"
                          "for the LES inflow module:\n"
                          "<%s>\n"
                          "The calculation will not be run.\n"), nomrub1);

            /* wave vector */
            char nomrub2[] = "batten_wave_vector";

            nbvent  = 3*inflow->n_modes;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            ierror = cs_restart_read_section(suite,
                                             nomrub2,
                                             support,
                                             nbvent,
                                             typ_val,
                                             inflow->wave_vector);

            if (ierror < CS_RESTART_SUCCESS)
              bft_error(__FILE__, __LINE__, 0,
                        _("Problem while reading section in the restart file\n"
                          "for the LES inflow module:\n"
                          "<%s>\n"
                          "The calculation will not be run.\n"), nomrub2);

            /* amplitude cos */
            char nomrub3[] = "batten_amplitude_cos";

            nbvent  = 3*inflow->n_modes;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            ierror = cs_restart_read_section(suite,
                                             nomrub3,
                                             support,
                                             nbvent,
                                             typ_val,
                                             inflow->amplitude_cos);

            if (ierror < CS_RESTART_SUCCESS)
              bft_error(__FILE__, __LINE__, 0,
                        _("Problem while reading section in the restart file\n"
                          "for the LES inflow module:\n"
                          "<%s>\n"
                          "The calculation will not be run.\n"), nomrub3);

            /* amplitude sin */
            char nomrub4[] = "batten_amplitude_sin";

            nbvent  = 3*inflow->n_modes;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            ierror = cs_restart_read_section(suite,
                                             nomrub4,
                                             support,
                                             nbvent,
                                             typ_val,
                                             inflow->amplitude_sin);

            if (ierror < CS_RESTART_SUCCESS)
              bft_error(__FILE__, __LINE__, 0,
                        _("Problem while reading section in the restart file\n"
                          "for the LES inflow module:\n"
                          "<%s>\n"
                          "The calculation will not be run.\n"), nomrub4);

          }

        }
        break;

      case CS_INFLOW_SEM:

        {
          cs_inflow_sem_t *inflow = (cs_inflow_sem_t *) inlet->inflow;

          { /* number of structures */
            char nomrub[] = "sem_number_structures";
            cs_int_t n_structures;

            BFT_MALLOC(tabvar, 1, cs_int_t);

            nbvent  = 1;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_int_t;

            ierror = cs_restart_read_section(suite,
                                             nomrub,
                                             support,
                                             nbvent,
                                             typ_val,
                                             tabvar);

            if (ierror < CS_RESTART_SUCCESS)
              bft_error(__FILE__, __LINE__, 0,
                        _("Problem while reading section in the restart file\n"
                          "for the LES inflow module:\n"
                          "<%s>\n"
                          "The calculation will not be run.\n"), nomrub);

            n_structures = *tabvar;

            /* Coherence check */
            if (inflow->n_structures != n_structures)
              bft_error(__FILE__, __LINE__, 0,
                        _("Stop reading the LES inflow module restart file.\n"
                          "%d eddies are given for the SEM "
                          "while the restart file contains %d.\n"),
                        inflow->n_structures, n_structures);

            BFT_FREE(tabvar);
          }

          { /* positions */
            char nomrub1[] = "sem_positions";

            nbvent  = 3*inflow->n_structures;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            ierror = cs_restart_read_section(suite,
                                             nomrub1,
                                             support,
                                             nbvent,
                                             typ_val,
                                             inflow->position);

            if (ierror < CS_RESTART_SUCCESS)
              bft_error(__FILE__, __LINE__, 0,
                        _("Problem while reading section in the restart file\n"
                          "for the LES inflow module:\n"
                          "<%s>\n"
                          "The calculation will not be run.\n"), nomrub1);

            /* energies */
            char nomrub2[] = "sem_energies";

            nbvent  = 3*inflow->n_structures;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            ierror = cs_restart_read_section(suite,
                                             nomrub2,
                                             support,
                                             nbvent,
                                             typ_val,
                                             inflow->energy);

            if (ierror < CS_RESTART_SUCCESS)
              bft_error(__FILE__, __LINE__, 0,
                        _("Problem while reading section in the restart file\n"
                          "for the LES inflow module:\n"
                          "<%s>\n"
                          "The calculation will not be run.\n"), nomrub2);

          }

        }
        break;

      }

      inlet->initialize = 0;

    }

  }

  cs_restart_read_fields(suite, CS_RESTART_LES_INFLOW);

  /* Close the restart file and free structures */
  cs_restart_destroy(&cs_glob_inflow_suite);

  bft_printf(_(" ...completed\n"));
}

/*----------------------------------------------------------------------------
 * Write the restart file of les inflow module
 *
 * Fortran interface:
 *
 * SUBROUTINE ECRSYN
 * *****************
 *
 * character(kind=c_char)  filename : <-- : Name of the restart file
 *----------------------------------------------------------------------------*/

void CS_PROCF(ecrsyn, ECRSYN)
(
 const char  *filename
)
{
  int   inlet_id;

  cs_int_t            nbvent;

  cs_restart_t             *suite;
  cs_mesh_location_type_t   support;
  cs_restart_val_type_t     typ_val;

  if (cs_glob_inflow_n_inlets == 0)
    return;

  bft_printf(_("\n Writing the LES inflow module restart file...\n"));

  /* Open the restart file */
  cs_glob_inflow_suite
    = cs_restart_create(filename, NULL, CS_RESTART_MODE_WRITE);

  if (cs_glob_inflow_suite == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while opening the LES inflow module restart "
                "file in write mode.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              filename);

  /* Pointer to the global restart structure */
  suite = cs_glob_inflow_suite;


  { /* Write the header */
    char       nomrub[] = "version_fichier_suite_turbulence_synthetique";
    cs_int_t   *tabvar;

    BFT_MALLOC(tabvar, 1, cs_int_t);

    nbvent  = 1;
    support = CS_MESH_LOCATION_NONE;
    typ_val = CS_TYPE_cs_int_t;

    *tabvar = 120;

    cs_restart_write_section(suite,
                             nomrub,
                             support,
                             nbvent,
                             typ_val,
                             tabvar);

    BFT_FREE(tabvar);
  }

  { /* Write the number of inlets */
    char       nomrub[] = "nb_inlets";
    cs_int_t   *tabvar;

    BFT_MALLOC(tabvar, 1, cs_int_t);

    *tabvar = cs_glob_inflow_n_inlets;

    nbvent  = 1;
    support = CS_MESH_LOCATION_NONE;
    typ_val = CS_TYPE_cs_int_t;

    cs_restart_write_section(suite,
                             nomrub,
                             support,
                             nbvent,
                             typ_val,
                             tabvar);

    BFT_FREE(tabvar);
  }

  { /* Write the structure of each inlet */
    cs_int_t   *tabvar;

    for (inlet_id = 0; inlet_id < cs_glob_inflow_n_inlets; inlet_id++) {

      cs_inlet_t *inlet = cs_glob_inflow_inlet_array[inlet_id];

      { /* type of inlet */
        char       nomrub[] = "type_inlet";

        BFT_MALLOC(tabvar, 1, cs_int_t);

        *tabvar = (cs_int_t)inlet->type;

        nbvent  = 1;
        support = CS_MESH_LOCATION_NONE;
        typ_val = CS_TYPE_cs_int_t;

        cs_restart_write_section(suite,
                                 nomrub,
                                 support,
                                 nbvent,
                                 typ_val,
                                 tabvar);

        BFT_FREE(tabvar);
      }

      switch(inlet->type) {

      case CS_INFLOW_LAMINAR:
        break;

      case CS_INFLOW_RANDOM:
        break;

      case CS_INFLOW_BATTEN:
        {
          cs_inflow_batten_t *inflow = (cs_inflow_batten_t *) inlet->inflow;

          { /* number of modes */
            char nomrub[] = "batten_number_modes";

            BFT_MALLOC(tabvar, 1, cs_int_t);

            *tabvar = inflow->n_modes;

            nbvent  = 1;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_int_t;

            cs_restart_write_section(suite,
                                     nomrub,
                                     support,
                                     nbvent,
                                     typ_val,
                                     tabvar);

            BFT_FREE(tabvar);
          }

          { /* frequencies */
            char nomrub1[] = "batten_frequencies";

            nbvent  = inflow->n_modes;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            cs_restart_write_section(suite,
                                     nomrub1,
                                     support,
                                     nbvent,
                                     typ_val,
                                     inflow->frequency);

            /* wave vector */
            char nomrub2[] = "batten_wave_vector";

            nbvent  = 3*inflow->n_modes;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            cs_restart_write_section(suite,
                                     nomrub2,
                                     support,
                                     nbvent,
                                     typ_val,
                                     inflow->wave_vector);

            /* amplitude cos */
            char nomrub3[] = "batten_amplitude_cos";

            nbvent  = 3*inflow->n_modes;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            cs_restart_write_section(suite,
                                     nomrub3,
                                     support,
                                     nbvent,
                                     typ_val,
                                     inflow->amplitude_cos);

            /* amplitude sin */
            char nomrub4[] = "batten_amplitude_sin";

            nbvent  = 3*inflow->n_modes;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            cs_restart_write_section(suite,
                                     nomrub4,
                                     support,
                                     nbvent,
                                     typ_val,
                                     inflow->amplitude_sin);
          }

        }
        break;

      case CS_INFLOW_SEM:

        {
          cs_inflow_sem_t *inflow = (cs_inflow_sem_t *) inlet->inflow;

          { /* number of structures */
            char nomrub[] = "sem_number_structures";

            BFT_MALLOC(tabvar, 1, cs_int_t);

            *tabvar = inflow->n_structures;

            nbvent  = 1;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_int_t;

            cs_restart_write_section(suite,
                                     nomrub,
                                     support,
                                     nbvent,
                                     typ_val,
                                     tabvar);

            BFT_FREE(tabvar);
          }

          { /* positions */
            char nomrub1[] = "sem_positions";

            nbvent  = 3*inflow->n_structures;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            cs_restart_write_section(suite,
                                     nomrub1,
                                     support,
                                     nbvent,
                                     typ_val,
                                     inflow->position);

            /* energies */
            char nomrub2[] = "sem_energies";

            nbvent  = 3*inflow->n_structures;
            support = CS_MESH_LOCATION_NONE;
            typ_val = CS_TYPE_cs_real_t;

            cs_restart_write_section(suite,
                                     nomrub2,
                                     support,
                                     nbvent,
                                     typ_val,
                                     inflow->energy);
          }

        }
        break;

      }

    }

  }

  cs_restart_write_fields(suite, CS_RESTART_LES_INFLOW);

  /* Close the restart file and free structures */
  cs_restart_destroy(&cs_glob_inflow_suite);

  bft_printf(_(" ...completed\n"));

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Destroy cs_inlet_t structures
 *----------------------------------------------------------------------------*/

void
cs_inflow_finalize(void)
{
  int coo_id;
  int inlet_id;

  if (cs_glob_inflow_n_inlets == 0)
    return;

  /* Destruction of each inlet structure */
  /*-------------------------------------*/

  for (inlet_id = 0; inlet_id < cs_glob_inflow_n_inlets; inlet_id++) {

    cs_inlet_t *inlet = cs_glob_inflow_inlet_array[inlet_id];

    bft_printf(_("\n"
                 "Summary of synthetic turbulence generation for inlet \"%d\""
                 " (%s) :\n\n"
                 "  Accumulated wall-clock time:      %12.3f\n"),
               inlet_id + 1, cs_inflow_type_name[inlet->type],
               inlet->wt_tot);

    if (cs_glob_rank_id < 0)
      bft_printf(_("  Accumulated CPU time:             %12.3f\n"),
                 inlet->cpu_tot);

#if defined(HAVE_MPI)

    else {

      double cpu_min, cpu_max, cpu_tot;
      double cpu_loc = inlet->cpu_tot;

      MPI_Allreduce(&cpu_loc, &cpu_min, 1, MPI_DOUBLE, MPI_MIN,
                    cs_glob_mpi_comm);
      MPI_Allreduce(&cpu_loc, &cpu_max, 1, MPI_DOUBLE, MPI_MAX,
                    cs_glob_mpi_comm);
      MPI_Allreduce(&cpu_loc, &cpu_tot, 1, MPI_DOUBLE, MPI_SUM,
                    cs_glob_mpi_comm);

      bft_printf(_("  Accumulated CPU time:\n"
                   "    local min:                      %12.3f\n"
                   "    local max:                      %12.3f\n"
                   "    mean:                           %12.3f\n"),
                 cpu_min, cpu_max, cpu_tot/cs_glob_n_ranks);

    }

#endif

    /* Mesh */

    if (inlet->n_faces > 0) {
    BFT_FREE(inlet->parent_num);
    BFT_FREE(inlet->face_centre);
    BFT_FREE(inlet->face_surface);
    }

    inlet->n_faces   = 0;

    /* Turbulence level */

    for (coo_id = 0; coo_id < 3; coo_id++)
      inlet->mean_velocity[coo_id] = 0.;

    inlet->kinetic_energy   = 0.;
    inlet->dissipation_rate = 0.;

    /* Generation method of synthetic turbulence */

    inlet->initialize = 0;

    switch(inlet->type) {

    case CS_INFLOW_LAMINAR:

      inlet->inflow = NULL;
      break;

    case CS_INFLOW_RANDOM:

      inlet->inflow = NULL;
      break;

    case CS_INFLOW_BATTEN:

      {
        cs_inflow_batten_t *inflow = (cs_inflow_batten_t *) inlet->inflow;

        inflow->n_modes = 0;

        BFT_FREE(inflow->frequency);
        BFT_FREE(inflow->wave_vector);
        BFT_FREE(inflow->amplitude_cos);
        BFT_FREE(inflow->amplitude_sin);

        BFT_FREE(inflow);

        inlet->inflow = NULL;
      }
      break;

    case CS_INFLOW_SEM:

      {
        cs_inflow_sem_t *inflow = (cs_inflow_sem_t *) inlet->inflow;;

        inflow->n_structures = 0;

        BFT_FREE(inflow->position);
        BFT_FREE(inflow->energy);

        BFT_FREE(inflow);

        inlet->inflow = NULL;
      }
      break;

    default:
      break;

    }

    inlet->wt_tot  = 0.;
    inlet->cpu_tot = 0.;

    BFT_FREE(inlet);
  }

  /* Global array of inlets */
  /*------------------------*/

  cs_glob_inflow_n_inlets = 0;
  BFT_FREE(cs_glob_inflow_inlet_array);
}

#undef _MODULE_3D
#undef _DOT_PRODUCT_3D
#undef _CROSS_PRODUCT_3D
#undef _TENSOR_TRACE_3D
#undef _INNER_PRODUCT_3D

/*----------------------------------------------------------------------------*/

END_C_DECLS
