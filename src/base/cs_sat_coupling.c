/*============================================================================
 * Functions associated with code coupling.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>
#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_nodal.h"
#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_coupling.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_connect.h"
#include "cs_prototypes.h"
#include "cs_selector.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sat_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Structure associated with Code_Saturne coupling */

typedef struct {

  int      match_id;        /* Id of matched application, -1 initially */
  char    *app_name;        /* Application name */
  char    *face_cpl_sel_c;  /* Face selection criteria */
  char    *cell_cpl_sel_c;  /* Cell selection criteria */
  char    *face_sup_sel_c;  /* Face selection criteria */
  char    *cell_sup_sel_c;  /* Cell selection criteria */
  int      verbosity;       /* Verbosity level */

} _cs_sat_coupling_builder_t;


struct _cs_sat_coupling_t {

  char            *sat_name;     /* Application name */

  char            *face_cpl_sel; /* Face selection criteria */
  char            *cell_cpl_sel; /* Face selection criteria */
  char            *face_sup_sel; /* Face selection criteria */
  char            *cell_sup_sel; /* Face selection criteria */

  ple_locator_t   *localis_cel;  /* Locator associated with cells */
  ple_locator_t   *localis_fbr;  /* Locator associated with boundary faces */

  cs_int_t         nbr_cel_sup;  /* Number of associated cell locations */
  cs_int_t         nbr_fbr_sup;  /* Number of associated face locations */

  fvm_nodal_t     *cells_sup;    /* Local cells at which distant values are
                                    interpolated*/
  fvm_nodal_t     *faces_sup;    /* Local faces at which distant values are
                                    interpolated*/

  cs_real_t       *distant_dist_fbr; /* Distant vectors (distance JJ') */
  cs_real_t       *distant_of;
  cs_real_t       *local_of;
  cs_real_t       *distant_pond_fbr; /* Distant weighting coefficient */
  cs_real_t       *local_pond_fbr;   /* Local weighting coefficient */

  int              verbosity; /* Verbosity level */

  /* Communication-related members */

#if defined(HAVE_MPI)

  MPI_Comm         comm;           /* Associated MPI communicator */

  int              n_sat_ranks;    /* Number of associated Code_Saturne ranks */
  int              sat_root_rank;  /* First associated Code_Saturne rank */

#endif

};

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Array of couplings */

static int _cs_glob_n_sat_cp = -1;

static int                         _sat_coupling_builder_size = 0;
static _cs_sat_coupling_builder_t *_sat_coupling_builder = NULL;

static int                  cs_glob_sat_n_couplings = 0;
static cs_sat_coupling_t  **cs_glob_sat_couplings = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Remove matched builder entries from the coupling builder.
 *----------------------------------------------------------------------------*/

static void
_remove_matched_builder_entries(void)
{
  int i;
  int n_unmatched_entries = 0;

  /* First, free arrays associated with marked entries */

  for (i = 0; i < _sat_coupling_builder_size; i++) {

    _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + i;

    if (scb->match_id > -1) {
      if (scb->face_cpl_sel_c != NULL) BFT_FREE(scb->face_cpl_sel_c);
      if (scb->cell_cpl_sel_c != NULL) BFT_FREE(scb->cell_cpl_sel_c);
      if (scb->face_sup_sel_c != NULL) BFT_FREE(scb->face_sup_sel_c);
      if (scb->cell_sup_sel_c != NULL) BFT_FREE(scb->cell_sup_sel_c);
      if (scb->app_name != NULL) BFT_FREE(scb->app_name);
    }
  }

  /* Now, remove marked entries and resize */

  for (i = 0; i < _sat_coupling_builder_size; i++) {
    _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + i;
    if (scb->match_id < 0) {
      *(_sat_coupling_builder + n_unmatched_entries) = *scb;
      n_unmatched_entries += 1;
    }
  }

  _sat_coupling_builder_size = n_unmatched_entries;

  BFT_REALLOC(_sat_coupling_builder,
              _sat_coupling_builder_size,
              _cs_sat_coupling_builder_t);
}

/*----------------------------------------------------------------------------
 * Print information on yet unmatched Code_Saturne couplings.
 *----------------------------------------------------------------------------*/

static void
_print_all_unmatched_sat(void)
{
  int i;

  const char empty_string[] = "";

  /* Loop on defined Code_Saturne instances */

  for (i = 0; i < _sat_coupling_builder_size; i++) {

    _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + i;

    if (scb->match_id < 0) {

      const char *local_name = empty_string;

      if (scb->app_name != NULL)
        local_name = scb->app_name;

      bft_printf(_(" Code_Saturne coupling:\n"
                   "   coupling id:              %d\n"
                   "   local name:               \"%s\"\n\n"),
                 i, local_name);
    }
  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Initialize communicator for Code_Saturne coupling
 *
 * parameters:
 *   sat_coupling  <-> Code_Saturne coupling structure
 *   coupling_id   <-- id of this coupling (for log file message)
 *----------------------------------------------------------------------------*/

static void
_init_comm(cs_sat_coupling_t *sat_coupling,
           int                coupling_id)

{
#if defined(HAVE_MPI)

  int  mpi_flag = 0;
  int local_range[2] = {-1, -1};
  int distant_range[2] = {-1, -1};

  MPI_Initialized(&mpi_flag);

  if (mpi_flag == 0)
    return;

  bft_printf(_(" Code_Saturne coupling %d: initializing MPI communication ... "),
             coupling_id);
  bft_printf_flush();

  ple_coupling_mpi_intracomm_create(MPI_COMM_WORLD,
                                    cs_glob_mpi_comm,
                                    sat_coupling->sat_root_rank,
                                    &(sat_coupling->comm),
                                    local_range,
                                    distant_range);

  bft_printf(_("[ok]\n"));
  bft_printf(_("  Local ranks = [%d..%d], distant ranks = [%d..%d].\n\n"),
             local_range[0], local_range[1] - 1,
             distant_range[0], distant_range[1] - 1);
  bft_printf_flush();

  sat_coupling->n_sat_ranks = distant_range[1] - distant_range[0];
  sat_coupling->sat_root_rank = distant_range[0];

#endif
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Add a Code_Saturne coupling using MPI.
 *
 * parameters:
 *   builder_id    <-- Code_Saturne application id in coupling builder
 *   sat_root_rank <-- root rank associated with Code_Saturne
 *   n_sat_ranks   <-- number of ranks associated with Code_Saturne
 *----------------------------------------------------------------------------*/

static void
_sat_add_mpi(int builder_id,
             int sat_root_rank,
             int n_sat_ranks)
{
  cs_sat_coupling_t *sat_coupling = NULL;
  _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + builder_id;

  /* Similarly to SYRTHES 4, we might be able to add
     Code_Saturne couplings directly (without resorting
     to a temporary builder), then match communications */

  cs_sat_coupling_add(scb->face_cpl_sel_c,
                      scb->cell_cpl_sel_c,
                      scb->face_sup_sel_c,
                      scb->cell_sup_sel_c,
                      scb->app_name,
                      scb->verbosity);

  sat_coupling = cs_sat_coupling_by_id(cs_sat_coupling_n_couplings() - 1);

  sat_coupling->sat_root_rank = sat_root_rank;
  sat_coupling->n_sat_ranks = n_sat_ranks;

  _init_comm(sat_coupling,
             builder_id);
}

/*----------------------------------------------------------------------------
 * Print information on identified Code_Saturne couplings using MPI.
 *
 * This function requires coupling_builder information, and must thus
 * be called before removing matched builder entries.
 *----------------------------------------------------------------------------*/

static void
_print_all_mpi_sat(void)
{
  int i;

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();
  const char empty_string[] = "";

  /* Loop on defined Code_Saturne instances */

  for (i = 0; i < _sat_coupling_builder_size; i++) {

    _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + i;

    if (scb->match_id > -1) {

      const char *local_name = empty_string;
      const char *distant_name = empty_string;

      const ple_coupling_mpi_set_info_t
        ai = ple_coupling_mpi_set_get_info(mpi_apps, scb->match_id);

      if (scb->app_name != NULL)
        local_name = scb->app_name;
      if (ai.app_name != NULL)
        distant_name = ai.app_name;

      bft_printf(_(" Code_Saturne coupling:\n"
                   "   coupling id:              %d\n"
                   "   local name:               \"%s\"\n"
                   "   distant application name: \"%s\"\n"
                   "   MPI application id:       %d\n"
                   "   MPI root rank:            %d\n"
                   "   number of MPI ranks:      %d\n\n"),
                 i, local_name, distant_name,
                 scb->match_id, ai.root_rank, ai.n_ranks);
    }
  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Initialize MPI Code_Saturne couplings using MPI.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

static void
_init_all_mpi_sat(void)
{
  int i;

  int n_apps = 0;
  int n_matched_apps = 0;
  int n_sat_apps = 0;

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == NULL)
    return;

  n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);

  /* First pass to count available Code_Saturne couplings */

  for (i = 0; i < n_apps; i++) {
    const ple_coupling_mpi_set_info_t
      ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
    if (strncmp(ai.app_type, "Code_Saturne", 12) == 0)
      n_sat_apps += 1;
  }

  /* In single-coupling mode, no identification necessary */

  if (n_sat_apps == 2 && _sat_coupling_builder_size == 1) {

    const int local_app_id = ple_coupling_mpi_set_get_app_id(mpi_apps);

    for (i = 0; i < n_apps; i++) {
      const ple_coupling_mpi_set_info_t
        ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
      if (   strncmp(ai.app_type, "Code_Saturne", 12) == 0
          && i != local_app_id) {
        _sat_coupling_builder->match_id = i;
      }
    }

    n_matched_apps += 1;

  }

  /* In multiple-coupling mode, identification is necessary */

  else {

    int j;
    ple_coupling_mpi_set_info_t ai;

    int *sat_appinfo = NULL;

    /* First, build an array of matched/unmatched Code_Saturne applications,
       with 2 entries per instance: matched indicator, app_id */

    BFT_MALLOC(sat_appinfo, n_sat_apps*2, int);

    n_sat_apps = 0;

    for (i = 0; i < n_apps; i++) {
      ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
      if (strncmp(ai.app_type, "Code_Saturne", 12) == 0) {
        sat_appinfo[n_sat_apps*2] = 0;
        sat_appinfo[n_sat_apps*2 + 1] = i;
        n_sat_apps += 1;
      }
    }

    /* Loop on defined Code_Saturne instances */

    for (i = 0; i < _sat_coupling_builder_size; i++) {

      _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + i;

      /* Loop on available Code_Saturne instances to match app_names */

      if (scb->app_name != NULL) {

        for (j = 0; j < n_sat_apps; j++) {

          if (sat_appinfo[j*2] != 0) /* Consider only unmatched applications */
            continue;

          ai = ple_coupling_mpi_set_get_info(mpi_apps, sat_appinfo[j*2 + 1]);
          if (ai.app_name == NULL)
            continue;

          if (strcmp(ai.app_name, scb->app_name) == 0) {
            scb->match_id = sat_appinfo[j*2 + 1];
            sat_appinfo[j*2] = i;
            n_matched_apps += 1;
            break;
          }

        }

      }

    } /* End of loop on defined Code_Saturne instances */

    BFT_FREE(sat_appinfo);

  } /* End of test on single or multiple Code_Saturne matching algorithm */

  /* Print matching info */

  _print_all_mpi_sat();

  /* Now initialize matched couplings */
  /*----------------------------------*/

  for (i = 0; i < _sat_coupling_builder_size; i++) {

    _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + i;

    if (scb->match_id > -1) {
      const ple_coupling_mpi_set_info_t
        ai = ple_coupling_mpi_set_get_info(mpi_apps, scb->match_id);

      if (strncmp(ai.app_type, "Code_Saturne", 12) == 0)
        _sat_add_mpi(i, ai.root_rank, ai.n_ranks);
    }

  }

  /* Cleanup */

  _remove_matched_builder_entries();
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Destroy a coupling structure
 *
 * parameters:
 *   couplage <-> pointer to coupling structure to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

static cs_sat_coupling_t *
_sat_coupling_destroy(cs_sat_coupling_t  *couplage)
{
  BFT_FREE(couplage->sat_name);

  BFT_FREE(couplage->face_cpl_sel);
  BFT_FREE(couplage->cell_cpl_sel);
  BFT_FREE(couplage->face_sup_sel);
  BFT_FREE(couplage->cell_sup_sel);

  ple_locator_destroy(couplage->localis_cel);
  ple_locator_destroy(couplage->localis_fbr);

  if (couplage->cells_sup != NULL)
    fvm_nodal_destroy(couplage->cells_sup);
  if (couplage->faces_sup != NULL)
    fvm_nodal_destroy(couplage->faces_sup);

  BFT_FREE(couplage->distant_dist_fbr);
  BFT_FREE(couplage->distant_of);
  BFT_FREE(couplage->local_of);
  BFT_FREE(couplage->distant_pond_fbr);
  BFT_FREE(couplage->local_pond_fbr);

#if defined(HAVE_MPI)
  if (   couplage->comm != MPI_COMM_WORLD
      && couplage->comm != cs_glob_mpi_comm)
    MPI_Comm_free(&(couplage->comm));
#endif

  BFT_FREE(couplage);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Computed some quantities needed for a centered-like interpolation
 *  - distance JJ' for distant boundary faces
 *  - local weighting coefficients
 *----------------------------------------------------------------------------*/

static void
_sat_coupling_interpolate(cs_sat_coupling_t  *couplage)
{
  int    icoo;
  int    reverse;

  cs_int_t    ind;
  cs_int_t    iel;
  cs_int_t    ifac;

  cs_int_t    n_fbr_loc  = 0;
  cs_int_t    n_fbr_dist = 0;

  cs_real_t   pdt_scal;
  cs_real_t   surface;

  cs_real_t   distance_fbr_cel;
  cs_real_t   distance_cel_cel;

  cs_real_t   dist_cel_fbr[3];
  cs_real_t   vect_surf_norm[3];

  cs_real_t  *local_surf     = NULL;
  cs_real_t  *local_xyzcen   = NULL;
  cs_real_t  *distant_surf   = NULL;
  cs_real_t  *distant_xyzcen = NULL;

  const cs_lnum_t   *lstfbr        = NULL;
  const cs_lnum_t   *element       = NULL;
  const cs_coord_t  *distant_coord = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *mesh_quantities = cs_glob_mesh_quantities;

  /* Removing the connectivity and localization informations in case of
     coupling update */

  if (couplage->distant_dist_fbr != NULL)
    BFT_FREE(couplage->distant_dist_fbr);
  if (couplage->distant_of != NULL)
    BFT_FREE(couplage->distant_of);
  if (couplage->local_of != NULL)
    BFT_FREE(couplage->local_of);
  if (couplage->distant_pond_fbr != NULL)
    BFT_FREE(couplage->distant_pond_fbr);
  if (couplage->local_pond_fbr != NULL)
    BFT_FREE(couplage->local_pond_fbr);


  /* Interpolation structure */

  n_fbr_loc  = ple_locator_get_n_interior(couplage->localis_fbr);
  lstfbr     = ple_locator_get_interior_list(couplage->localis_fbr);

  n_fbr_dist    = ple_locator_get_n_dist_points(couplage->localis_fbr);
  element       = ple_locator_get_dist_locations(couplage->localis_fbr);
  distant_coord = ple_locator_get_dist_coords(couplage->localis_fbr);


  /* Calculation of the distance DJJPB defining the distance from */
  /* the local cell center to the distant boundary face norm      */
  /*--------------------------------------------------------------*/

  BFT_MALLOC(couplage->distant_dist_fbr, 3*n_fbr_dist, cs_real_t);

  /* Store the local surface vector of the coupled boundary faces */

  BFT_MALLOC(local_surf, 3*n_fbr_loc, cs_real_t);

  for (ind = 0 ; ind < n_fbr_loc ; ind++) {

    ifac = lstfbr[ind] - 1;

    for (icoo = 0 ; icoo < 3 ; icoo++)
      local_surf[ind*3 + icoo] = mesh_quantities->b_face_normal[ifac*3 + icoo];

 }

  /* Get the distant faces surface vector (reverse = 1) */

  reverse = 1;

  BFT_MALLOC(distant_surf, 3*n_fbr_dist, cs_real_t);

  ple_locator_exchange_point_var(couplage->localis_fbr,
                                 distant_surf,
                                 local_surf,
                                 NULL,
                                 sizeof(cs_real_t),
                                 3,
                                 reverse);

  BFT_FREE(local_surf);

  /* Calculation of the JJ' vectors */

  BFT_MALLOC(distant_xyzcen, 3*n_fbr_dist, cs_real_t);

  for (ind = 0; ind < n_fbr_dist; ind++) {

    iel = element[ind] - 1;

    surface = 0.;
    for (icoo = 0; icoo < 3; icoo++)
      surface += distant_surf[ind*3 + icoo]*distant_surf[ind*3 + icoo];
    surface = sqrt(surface);

    pdt_scal = 0.;
    for (icoo = 0; icoo < 3; icoo++) {

      dist_cel_fbr[icoo] =
        distant_coord[ind*3 + icoo] - mesh_quantities->cell_cen[iel*3 + icoo];

      /* Store the distant coordinates to compute the weighting coefficients */
      distant_xyzcen[ind*3 + icoo] = mesh_quantities->cell_cen[iel*3 + icoo];

      vect_surf_norm[icoo] =
        distant_surf[ind*3 + icoo] / surface;

      pdt_scal += dist_cel_fbr[icoo]*vect_surf_norm[icoo];

    }

    for (icoo = 0; icoo < 3; icoo++)
      couplage->distant_dist_fbr[ind*3 + icoo] =
        dist_cel_fbr[icoo] - pdt_scal*vect_surf_norm[icoo];

  }

  BFT_FREE(distant_surf);


  /* Calculation of the local weighting coefficient */
  /*------------------------------------------------*/

  BFT_MALLOC(couplage->distant_pond_fbr, n_fbr_dist, cs_real_t);
  BFT_MALLOC(couplage->local_pond_fbr, n_fbr_loc, cs_real_t);

  /* Get the cell center coordinates (reverse = 0) */

  reverse = 0;

  BFT_MALLOC(local_xyzcen, 3*n_fbr_loc, cs_real_t);

  ple_locator_exchange_point_var(couplage->localis_fbr,
                                 distant_xyzcen,
                                 local_xyzcen,
                                 NULL,
                                 sizeof(cs_real_t),
                                 3,
                                 reverse);

  BFT_FREE(distant_xyzcen);

  /* Calculation of the local weighting coefficients */

  for (ind = 0 ; ind < n_fbr_loc ; ind++) {

    ifac = lstfbr[ind] - 1;
    iel  = mesh->b_face_cells[ifac];

    surface = 0.;

    distance_fbr_cel = 0.;
    distance_cel_cel = 0.;

    for (icoo = 0 ; icoo < 3 ; icoo++) {

      surface += mesh_quantities->b_face_normal[ifac*3 + icoo]
        * mesh_quantities->b_face_normal[ifac*3 + icoo];

      distance_fbr_cel += mesh_quantities->b_face_normal[ifac*3 + icoo] *
        (local_xyzcen[ind*3 + icoo] - mesh_quantities->b_face_cog[ifac*3 + icoo]);

      distance_cel_cel += mesh_quantities->b_face_normal[ifac*3 + icoo] *
        (local_xyzcen[ind*3 + icoo] - mesh_quantities->cell_cen[iel*3 + icoo]);

    }

    surface = sqrt(surface);

    distance_fbr_cel /= surface;
    distance_cel_cel /= surface;

    if (fabs(distance_cel_cel) > 1.e-12)
      couplage->local_pond_fbr[ind] = distance_fbr_cel / distance_cel_cel;
    else
      couplage->local_pond_fbr[ind] = 0.5;

  }


  /* Get the distant weighting coefficients (reverse = 1) */

  reverse = 1;

  ple_locator_exchange_point_var(couplage->localis_fbr,
                                 couplage->distant_pond_fbr,
                                 couplage->local_pond_fbr,
                                 NULL,
                                 sizeof(cs_real_t),
                                 1,
                                 reverse);



  /* Calculation of the OF distance */
  /*--------------------------------*/

  BFT_MALLOC(couplage->distant_of, 3*n_fbr_dist, cs_real_t);
  BFT_MALLOC(couplage->local_of, 3*n_fbr_loc, cs_real_t);

  for (ind = 0 ; ind < n_fbr_loc ; ind++) {

    ifac = lstfbr[ind] - 1;
    iel  = mesh->b_face_cells[ifac];

    surface = 0.;

    distance_fbr_cel = 0.;
    distance_cel_cel = 0.;

    for (icoo = 0 ; icoo < 3 ; icoo++) {

      surface += mesh_quantities->b_face_normal[ifac*3 + icoo]
        * mesh_quantities->b_face_normal[ifac*3 + icoo];

      distance_fbr_cel += mesh_quantities->b_face_normal[ifac*3 + icoo] *
        (local_xyzcen[ind*3+icoo] - mesh_quantities->b_face_cog[ifac*3+icoo]);

      distance_cel_cel += mesh_quantities->b_face_normal[ifac*3 + icoo] *
        (local_xyzcen[ind*3+icoo] - mesh_quantities->cell_cen[iel*3+icoo]);

    }

    surface = sqrt(surface);

    distance_fbr_cel /= surface;
    distance_cel_cel /= surface;

    for (icoo = 0 ; icoo < 3 ; icoo++)

      couplage->local_of[ind*3 + icoo] =
        mesh_quantities->b_face_cog[ifac*3 + icoo]
        -  (mesh_quantities->b_face_cog[ifac*3 + icoo] /*  O'  */
            + mesh_quantities->b_face_normal[ifac*3 + icoo] *distance_fbr_cel/surface   /*J'=F+n*FJ'*/
            - 0.5*mesh_quantities->b_face_normal[ifac*3 + icoo] *distance_cel_cel/surface );  /*-n*I'J'/2*/

  }

  reverse = 1;

  ple_locator_exchange_point_var(couplage->localis_fbr,
                                 couplage->distant_of,
                                 couplage->local_of,
                                 NULL,
                                 sizeof(cs_real_t),
                                 3,
                                 reverse);

  BFT_FREE(local_xyzcen);


}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * User function wrapper for definition of Code_Saturne couplings
 *
 * Fortran Interface:
 *
 * SUBROUTINE USSATC
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (ussatc, USSATC)
(
 void
)
{
  cs_user_saturne_coupling();
}

/*----------------------------------------------------------------------------
 * Get number of code coupling
 *
 * Fortran interface:
 *
 * SUBROUTINE NBCCPL
 * *****************
 *
 * INTEGER          NBRCPL         : <-- : number of code couplings
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbccpl, NBCCPL)
(
 cs_lnum_t   *n_couplings
)
{
  if (_cs_glob_n_sat_cp < 0) {
    if (_sat_coupling_builder_size > 0)
      _cs_glob_n_sat_cp = _sat_coupling_builder_size;
    else
      _cs_glob_n_sat_cp = cs_sat_coupling_n_couplings();
  }

  *n_couplings = _cs_glob_n_sat_cp;
}

/*----------------------------------------------------------------------------
 * Set the list of cells and boundary faces associated to a coupling
 * and a cloud of point.
 *
 * The local "support" cells and boundary faces are used to localize
 * the values in the distant "coupled" cells and faces.
 * Depending on the role of sender and/or receiver of the current process
 * in the coupling, some of these sets can be empty or not.
 *
 * The cell values are always localized and interpolated on the distant
 * "cells" support. The face values are localized and interpolated on
 * the distant "face" support if present, or on the distant "cell" support
 * if not.
 *
 * If the input arrays LCESUP and LFBSUP are not ordered, they will be
 * orderd in output.
 *
 * Fortran interface:
 *
 * SUBROUTINE DEFLOC
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 *----------------------------------------------------------------------------*/

void CS_PROCF (defloc, DEFLOC)
(
 const cs_int_t  *numcpl
)
{
  cs_int_t  ind;
  cs_int_t  nbr_fbr_cpl = 0, nbr_cel_cpl = 0;

  int  indic_glob[2] = {0, 0};
  int  indic_loc[2] = {0, 0};

  char coupled_mesh_name[64];
  cs_lnum_t *c_elt_list = NULL;
  cs_lnum_t *f_elt_list = NULL;
  cs_sat_coupling_t  *coupl = NULL;
  fvm_nodal_t  *support_fbr = NULL;
  cs_mesh_quantities_t  *mesh_quantities = cs_glob_mesh_quantities;

  int locator_options[PLE_LOCATOR_N_OPTIONS];
  locator_options[PLE_LOCATOR_NUMBERING] = 1;

  const double tolerance = 0.1;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  /* Removing the connectivity and localization informations in case of
     coupling update */

  if (coupl->cells_sup != NULL) fvm_nodal_destroy(coupl->cells_sup);
  if (coupl->faces_sup != NULL) fvm_nodal_destroy(coupl->faces_sup);

  /* Create the local lists */

  if (coupl->cell_sup_sel != NULL) {

    BFT_MALLOC(c_elt_list, cs_glob_mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_num_list(coupl->cell_sup_sel,
                                  &(coupl->nbr_cel_sup),
                                  c_elt_list);

  }

  if (coupl->face_sup_sel != NULL) {

    BFT_MALLOC(f_elt_list, cs_glob_mesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_num_list(coupl->face_sup_sel,
                                    &(coupl->nbr_fbr_sup),
                                    f_elt_list);

  }

  if (coupl->nbr_cel_sup > 0) indic_loc[0] = 1;
  if (coupl->nbr_fbr_sup > 0) indic_loc[1] = 1;

  for (ind = 0 ; ind < 2 ; ind++)
    indic_glob[ind] = indic_loc[ind];

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Allreduce (indic_loc, indic_glob, 2, MPI_INT, MPI_MAX,
                   cs_glob_mpi_comm);
#endif

  if (indic_glob[0] > 0) {

    sprintf(coupled_mesh_name, _("coupled_cells_%d"), *numcpl);

    coupl->cells_sup = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                                      coupled_mesh_name,
                                                      false,
                                                      coupl->nbr_cel_sup,
                                                      c_elt_list);

  }

  if (indic_glob[1] > 0) {

    sprintf(coupled_mesh_name, _("coupled_faces_%d"), *numcpl);

    coupl->faces_sup = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                      coupled_mesh_name,
                                                      false,
                                                      0,
                                                      coupl->nbr_fbr_sup,
                                                      NULL,
                                                      f_elt_list);

  }

  if (coupl->cell_sup_sel != NULL) BFT_FREE(c_elt_list);
  if (coupl->face_sup_sel != NULL) BFT_FREE(f_elt_list);

  /* Build and initialize associated locator */

#if defined(PLE_HAVE_MPI)

  if (coupl->localis_cel == NULL)
    coupl->localis_cel = ple_locator_create(coupl->comm,
                                            coupl->n_sat_ranks,
                                            coupl->sat_root_rank);

  if (coupl->localis_fbr == NULL)
    coupl->localis_fbr = ple_locator_create(coupl->comm,
                                            coupl->n_sat_ranks,
                                            coupl->sat_root_rank);

#else

  if (coupl->localis_cel == NULL)
    coupl->localis_cel = ple_locator_create();

  if (coupl->localis_fbr == NULL)
    coupl->localis_fbr = ple_locator_create();

#endif

  /* Initialization of the distant point localization */

  if (coupl->cell_cpl_sel != NULL) {

    BFT_MALLOC(c_elt_list, cs_glob_mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_num_list(coupl->cell_cpl_sel,
                                  &nbr_cel_cpl,
                                  c_elt_list);

  }

  ple_locator_set_mesh(coupl->localis_cel,
                       coupl->cells_sup,
                       locator_options,
                       0.,
                       tolerance,
                       3,
                       nbr_cel_cpl,
                       c_elt_list,
                       NULL,
                       mesh_quantities->cell_cen,
                       NULL,
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh_p);

  if (coupl->cell_cpl_sel != NULL) BFT_FREE(c_elt_list);


  if (coupl->face_cpl_sel != NULL) {

    BFT_MALLOC(f_elt_list, cs_glob_mesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_num_list(coupl->face_cpl_sel,
                                    &nbr_fbr_cpl,
                                    f_elt_list);

  }

  if (indic_glob[1] > 0)
    support_fbr = coupl->faces_sup;
  else
    support_fbr = coupl->cells_sup;

  ple_locator_set_mesh(coupl->localis_fbr,
                       support_fbr,
                       locator_options,
                       0.,
                       tolerance,
                       3,
                       nbr_fbr_cpl,
                       f_elt_list,
                       NULL,
                       mesh_quantities->b_face_cog,
                       NULL,
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh_p);

  if (coupl->face_cpl_sel != NULL) BFT_FREE(f_elt_list);


  /* Computed some quantities needed for a centered-like interpolation */

  if (coupl->localis_fbr != NULL)
    _sat_coupling_interpolate(coupl);


#if 0
  /* TODO: associate the FVM meshes to the post-processing,
     with a fonction giving a pointer to the associated FVM structures,
     and another enabling its compacting or removing */
  {
    fvm_writer_t *w = fvm_writer_init("coupled_mesh",
                                      NULL,
                                      "EnSight Gold",
                                      "binary",
                                      FVM_WRITER_FIXED_MESH);

    fvm_writer_export_nodal(w, coupl->cells_sup);
    fvm_writer_finalize(w);

  }
#endif

  /* Compacting the interpolation support (could be removed) */

  if (coupl->cells_sup != NULL)
    fvm_nodal_reduce(coupl->cells_sup, 1);
  if (coupl->faces_sup != NULL)
    fvm_nodal_reduce(coupl->faces_sup, 1);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  ple_locator_dump(coupl->localis_cel);
  ple_locator_dump(coupl->localis_fbr);
#endif
}

/*----------------------------------------------------------------------------
 * Get the number of cells and boundary faces, "support", coupled and not
 * localized associated to a given coupling
 *
 * Fortran interface:
 *
 * SUBROUTINE NBECPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCESUP         : <-- : number of "support" cells
 * INTEGER          NFBSUP         : <-- : number of "support" boundary faces
 * INTEGER          NCECPL         : <-- : number of coupled cells
 * INTEGER          NFBCPL         : <-- : number of coupled boundary faces
 * INTEGER          NCENCP         : <-- : number of not coupled cells
 *                                 :     : (since not localized)
 * INTEGER          NFBNCP         : <-- : number of not coupled boundary faces
 *                                 :     : (since not localized)
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbecpl, NBECPL)
(
 const cs_int_t  *numcpl,
       cs_int_t  *ncesup,
       cs_int_t  *nfbsup,
       cs_int_t  *ncecpl,
       cs_int_t  *nfbcpl,
       cs_int_t  *ncencp,
       cs_int_t  *nfbncp
)
{
  cs_sat_coupling_t  *coupl = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  *ncesup = coupl->nbr_cel_sup;
  *nfbsup = coupl->nbr_fbr_sup;

  *ncecpl = 0;
  *nfbcpl = 0;

  *ncencp = 0;
  *nfbncp = 0;

  if (coupl->localis_cel != NULL) {
    *ncecpl = ple_locator_get_n_interior(coupl->localis_cel);
    *ncencp = ple_locator_get_n_exterior(coupl->localis_cel);
  }

  if (coupl->localis_fbr != NULL) {
    *nfbcpl = ple_locator_get_n_interior(coupl->localis_fbr);
    *nfbncp = ple_locator_get_n_exterior(coupl->localis_fbr);
  }

}

/*----------------------------------------------------------------------------
 * Get the lists of coupled cells and boundary faces (i.e. receiving)
 * associated to a given coupling
 *
 * The number of cells and boundary faces, got with NBECPL(), are used
 * for arguments coherency checks.
 *
 * Fortran interface:
 *
 * SUBROUTINE LELCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCECPL         : --> : number of coupled cells
 * INTEGER          NFBCPL         : --> : number of coupled boundary faces
 * INTEGER          LCECPL(*)      : <-- : list of coupled cells
 * INTEGER          LFBCPL(*)      : <-- : list of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (lelcpl, LELCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *ncecpl,
 const cs_int_t  *nfbcpl,
       cs_int_t  *lcecpl,
       cs_int_t  *lfbcpl
)
{
  cs_int_t  ind;

  cs_int_t  _ncecpl = 0;
  cs_int_t  _nfbcpl = 0;

  cs_sat_coupling_t  *coupl = NULL;

  const cs_int_t  *lst = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  if (coupl->localis_cel != NULL)
    _ncecpl = ple_locator_get_n_interior(coupl->localis_cel);

  if (coupl->localis_fbr != NULL)
    _nfbcpl = ple_locator_get_n_interior(coupl->localis_fbr);

  if (*ncecpl != _ncecpl || *nfbcpl != _nfbcpl)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for LELCPL()\n"
                "NCECPL = %d and NFBCPL = %d are indicated.\n"
                "The values for this coupling should be %d and %d."),
              *numcpl, (int)(*ncecpl), (int)(*nfbcpl),
              (int)_ncecpl, (int)_nfbcpl);

  /* Copy lists (would be useless with a pure C API) */

  if (_ncecpl > 0) {
    lst = ple_locator_get_interior_list(coupl->localis_cel);
    for (ind = 0 ; ind < _ncecpl ; ind++)
      lcecpl[ind] = lst[ind];
  }

  if (_nfbcpl > 0) {
    lst = ple_locator_get_interior_list(coupl->localis_fbr);
    for (ind = 0 ; ind < _nfbcpl ; ind++)
      lfbcpl[ind] = lst[ind];
  }
}

/*----------------------------------------------------------------------------
 * Get the lists of not coupled cells and boundary faces (i.e. receiving but
 * not localized) associated to a given coupling
 *
 * The number of cells and boundary faces, got with NBECPL(), are used
 * for arguments coherency checks.
 *
 * Fortran interface:
 *
 * SUBROUTINE LENCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCENCP         : --> : number of not coupled cells
 * INTEGER          NFBNCP         : --> : number of not coupled boundary faces
 * INTEGER          LCENCP(*)      : <-- : list of not coupled cells
 * INTEGER          LFBNCP(*)      : <-- : list of not coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (lencpl, LENCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *ncencp,
 const cs_int_t  *nfbncp,
       cs_int_t  *lcencp,
       cs_int_t  *lfbncp
)
{
  cs_int_t  ind;

  cs_int_t  _ncencp = 0;
  cs_int_t  _nfbncp = 0;
  cs_sat_coupling_t  *coupl = NULL;

  const cs_int_t  *lst = NULL;


  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  if (coupl->localis_cel != NULL)
    _ncencp = ple_locator_get_n_exterior(coupl->localis_cel);

  if (coupl->localis_fbr != NULL)
    _nfbncp = ple_locator_get_n_exterior(coupl->localis_fbr);

  if (*ncencp != _ncencp || *nfbncp != _nfbncp)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for LELNCP()\n"
                "NCENCP = %d and NFBNCP = %d are indicated.\n"
                "The values for this coupling should be %d and %d."),
              *numcpl, (int)(*ncencp), (int)(*nfbncp),
              (int)_ncencp, (int)_nfbncp);

  /* Copy lists (would be useless with a pure C API) */

  if (_ncencp > 0) {
    lst = ple_locator_get_exterior_list(coupl->localis_cel);
    for (ind = 0 ; ind < _ncencp ; ind++)
      lcencp[ind] = lst[ind];
  }

  if (_nfbncp > 0) {
    lst = ple_locator_get_exterior_list(coupl->localis_fbr);
    for (ind = 0 ; ind < _nfbncp ; ind++)
      lfbncp[ind] = lst[ind];
  }
}

/*----------------------------------------------------------------------------
 * Get the number of distant point associated to a given coupling
 * and localized on the local domain
 *
 * Fortran interface:
 *
 * SUBROUTINE NPDCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCEDIS         : <-- : number of distant cells
 * INTEGER          NFBDIS         : <-- : numbre de distant boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (npdcpl, NPDCPL)
(
 const cs_int_t  *numcpl,
       cs_int_t  *ncedis,
       cs_int_t  *nfbdis
)
{
  cs_sat_coupling_t  *coupl = NULL;

  /* Verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  /* Get the number of points */

  *ncedis = 0;
  *nfbdis = 0;

  if (coupl->localis_cel != NULL)
    *ncedis = ple_locator_get_n_dist_points(coupl->localis_cel);

  if (coupl->localis_fbr != NULL)
    *nfbdis = ple_locator_get_n_dist_points(coupl->localis_fbr);

}

/*----------------------------------------------------------------------------
 * Get the distant points coordinates associated to a given coupling
 * and a list of points, and the elements number and type (cell or face)
 * "containing" this points.
 *
 * The number of distant points NBRPTS must be equal to one the arguments
 * NCEDIS or NFBDIS given by NPDCPL(), and is given here for coherency checks
 * between the arguments NUMCPL and ITYSUP.
 *
 * Fortran interface:
 *
 * SUBROUTINE COOCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRPTS         : --> : number of distant points
 * INTEGER          ITYDIS         : --> : 1 : access to the points associated
 *                                 :     :     to the distant cells
 *                                 :     : 2 : access to the points associated
 *                                 :     :     to the distant boundary faces
 * INTEGER          ITYLOC         : <-- : 1 : localization on the local cells
 *                                 :     : 2 : localization on the local faces
 * INTEGER          LOCPTS(*)      : <-- : "containing" number associated to
 *                                 :     :   each point
 * DOUBLE PRECISION COOPTS(3,*)    : <-- : distant point coordinates
 * DOUBLE PRECISION DJPPTS(3,*)    : <-- : distant vectors to the coupled face
 * DOUBLE PRECISION PNDPTS(*)      : <-- : distant weighting coefficients
 *----------------------------------------------------------------------------*/

void CS_PROCF (coocpl, COOCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *nbrpts,
 const cs_int_t  *itydis,
       cs_int_t  *ityloc,
       cs_int_t  *locpts,
       cs_real_t *coopts,
       cs_real_t *djppts,
       cs_real_t *dofpts,
       cs_real_t *pndpts
)
{
  cs_int_t  ind, icoo;

  cs_int_t  n_pts_dist = 0;
  cs_sat_coupling_t  *coupl = NULL;
  ple_locator_t  *localis = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  *ityloc = 0;

  if (*itydis == 1) {
    localis = coupl->localis_cel;
    *ityloc = 1;
  }
  else if (*itydis == 2) {
    localis = coupl->localis_fbr;
    if (coupl->nbr_fbr_sup > 0)
      *ityloc = 2;
    else
      *ityloc = 1;
  }

  if (localis != NULL)
    n_pts_dist = ple_locator_get_n_dist_points(localis);

  if (*nbrpts != n_pts_dist)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for COOCPL()\n"
                "ITYDIS = %d and NBRPTS = %d are indicated.\n"
                "The value for NBRPTS should be %d."),
              *numcpl, (int)(*itydis), (int)(*nbrpts), (int)n_pts_dist);

  /* Creation the local lists */

  if (localis != NULL) {

    n_pts_dist = ple_locator_get_n_dist_points(localis);

    if (n_pts_dist > 0) {

      const cs_lnum_t   *element;
      const cs_coord_t  *coord;

      element = ple_locator_get_dist_locations(localis);
      coord   = ple_locator_get_dist_coords(localis);

      for (ind = 0 ; ind < n_pts_dist ; ind++) {
        locpts[ind] = element[ind];
        for (icoo = 0 ; icoo < 3 ; icoo++)
          coopts[ind*3 + icoo] = coord[ind*3 + icoo];
      }

      if (*itydis == 2)
        for (ind = 0 ; ind < n_pts_dist ; ind++)
          for (icoo = 0 ; icoo < 3 ; icoo++) {
            djppts[ind*3 + icoo] = coupl->distant_dist_fbr[ind*3 + icoo];
            dofpts[ind*3 + icoo] = coupl->distant_of[ind*3 + icoo];
            pndpts[ind] = coupl->distant_pond_fbr[ind];
          }

    }

  }

}

/*----------------------------------------------------------------------------
 * Get the weighting coefficient needed for a centered-like interpolation
 * in the case of a coupling on boundary faces.
 *
 * Fortran interface:
 *
 * SUBROUTINE PONDCP
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRCPL         : --> : number of distant points
 * INTEGER          ITYLOC         : <-- : 1 : localization on the local cells
 *                                 :     : 2 : localization on the local faces
 * DOUBLE PRECISION PNDCPL(*)      : <-- : weighting coefficients
 *----------------------------------------------------------------------------*/

void CS_PROCF (pondcp, PONDCP)
(
 const cs_int_t  *const numcpl,
 const cs_int_t  *const nbrpts,
       cs_int_t  *const ityloc,
       cs_real_t *const pndcpl,
       cs_real_t *const distof
)
{
  int             icoo;
  cs_int_t        ind;
  cs_int_t        nfbcpl = 0;
  cs_sat_coupling_t  *coupl = NULL;
  ple_locator_t  *localis = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  if (*ityloc == 1)
    bft_error(__FILE__, __LINE__, 0,
              _("The centered interpolation scheme is not available\n"
                "when coupling cells"));
  else if (*ityloc == 2)
    localis = coupl->localis_fbr;


  if (localis != NULL)
    nfbcpl = ple_locator_get_n_interior(localis);

  if (*nbrpts != nfbcpl)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for PNDCPL().\n"
                "ITYLOC = %d and NBRPTS = %d are indicated.\n"
                "NBRPTS should be %d."),
              *numcpl, (int)(*ityloc), (int)(*nbrpts), (int)nfbcpl);

  /* Creation of the local lists */

  if (localis != NULL) {

    if (nfbcpl > 0) {

      for (ind = 0 ; ind < nfbcpl ; ind++) {
        pndcpl[ind] = coupl->local_pond_fbr[ind];
        for (icoo = 0 ; icoo < 3 ; icoo++)
          distof[ind*3 + icoo] = coupl->local_of[ind*3 + icoo];
      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Exchange a variable associated to a set of point and a coupling.
 *
 * Fortran interface:
 *
 * SUBROUTINE VARCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRDIS         : --> : number of values to send
 * INTEGER          NBRLOC         : --> : number of values to receive
 * INTEGER          ITYVAR         : --> : 1 : variables defined at cells
 *                                 :     : 2 : variables defined at faces
 * INTEGER          STRIDE         : --> : 1 : for scalars
 *                                 :     : 3 : for vectors
 * DOUBLE PRECISION VARDIS(*)      : --> : distant variable(to send)
 * DOUBLE PRECISION VARLOC(*)      : <-- : local variable (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (varcpl, VARCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *nbrdis,
 const cs_int_t  *nbrloc,
 const cs_int_t  *ityvar,
 const cs_int_t  *stride,
       cs_real_t *vardis,
       cs_real_t *varloc
)
{
  cs_int_t  n_val_dist_ref = 0;
  cs_int_t  n_val_loc_ref = 0;
  cs_real_t  *val_dist = NULL;
  cs_real_t  *val_loc = NULL;
  cs_sat_coupling_t  *coupl = NULL;
  ple_locator_t  *localis = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  if (*ityvar == 1)
    localis = coupl->localis_cel;
  else if (*ityvar == 2)
    localis = coupl->localis_fbr;

  if (localis != NULL) {
    n_val_dist_ref = ple_locator_get_n_dist_points(localis);
    n_val_loc_ref  = ple_locator_get_n_interior(localis);
  }

  if (*nbrdis > 0 && *nbrdis != n_val_dist_ref)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for VARCPL()\n"
                "ITYVAR = %d and NBRDIS = %d are indicated.\n"
                "NBRDIS should be 0 or %d."),
              *numcpl, (int)(*ityvar), (int)(*nbrdis), (int)n_val_dist_ref);

  if (*nbrloc > 0 && *nbrloc != n_val_loc_ref)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for VARCPL()\n"
                "ITYVAR = %d and NBRLOC = %d are indicated.\n"
                "NBRLOC should be 0 or %d."),
              *numcpl, (int)(*ityvar), (int)(*nbrloc), (int)n_val_loc_ref);

  /* Create the local lists */

  if (localis != NULL) {

    if (*nbrdis > 0)
      val_dist = vardis;
    if (*nbrloc > 0)
      val_loc = varloc;

    ple_locator_exchange_point_var(localis,
                                   val_dist,
                                   val_loc,
                                   NULL,
                                   sizeof(cs_real_t),
                                   *stride,
                                   0);

  }

}

/*----------------------------------------------------------------------------
 * Array of integers exchange, associated to a given coupling.
 *
 * It is assumed that the arrays have the same size and the same values on
 * each group of processus (local and distant).
 *
 * Fortran interface:
 *
 * SUBROUTINE TBICPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRDIS         : --> : number of values to send
 * INTEGER          NBRLOC         : --> : number of values to receive
 * INTEGER          TABDIS(*)      : --> : distant values (to send)
 * INTEGER          TABLOC(*)      : <-- : local values (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tbicpl, TBICPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *nbrdis,
 const cs_int_t  *nbrloc,
       cs_int_t  *vardis,
       cs_int_t  *varloc
)
{
  cs_int_t  ind;
  cs_int_t  nbr = 0;
  bool  distant = false;

#if defined(HAVE_MPI)

  MPI_Status  status;

  cs_sat_coupling_t  *coupl = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  if (coupl->comm != MPI_COMM_NULL) {

    distant = true;

    /* Exchange between the groups master node */

    if (cs_glob_rank_id < 1)
      MPI_Sendrecv(vardis, *nbrdis, CS_MPI_INT, coupl->sat_root_rank, 0,
                   varloc, *nbrloc, CS_MPI_INT, coupl->sat_root_rank, 0,
                   coupl->comm, &status);

    /* Synchronization inside a group */

    if (cs_glob_n_ranks > 1)
      MPI_Bcast (varloc, *nbrloc, CS_MPI_INT, 0, cs_glob_mpi_comm);

  }

#endif /* defined(HAVE_MPI) */

  if (distant == false) {

    nbr = CS_MIN(*nbrdis, *nbrloc);

    for (ind = 0; ind < nbr; ind++)
      varloc[ind] = vardis[ind];

  }
}

/*----------------------------------------------------------------------------
 * Array of reals exchange, associated to a given coupling.
 *
 * It is assumed that the arrays have the same size and the same values on
 * each group of processus (local and distant).
 *
 * Fortran interface:
 *
 * SUBROUTINE TBRCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRDIS         : --> : number of values to send
 * INTEGER          NBRLOC         : --> : number of values to receive
 * DOUBLE PRECISION TABDIS(*)      : --> : distant values (to send)
 * DOUBLE PRECISION TABLOC(*)      : <-- : local values (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tbrcpl, TBRCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *nbrdis,
 const cs_int_t  *nbrloc,
       cs_real_t *vardis,
       cs_real_t *varloc
)
{
  cs_int_t  ind;
  cs_int_t  nbr = 0;
  bool  distant = false;

#if defined(HAVE_MPI)

  MPI_Status  status;

  cs_sat_coupling_t  *coupl = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  if (coupl->comm != MPI_COMM_NULL) {

    distant = true;

    /* Exchange between the groups master node */

    if (cs_glob_rank_id < 1)
      MPI_Sendrecv(vardis, *nbrdis, CS_MPI_REAL, coupl->sat_root_rank, 0,
                   varloc, *nbrloc, CS_MPI_REAL, coupl->sat_root_rank, 0,
                   coupl->comm, &status);

    /* Synchronization inside a group */

    if (cs_glob_n_ranks > 1)
      MPI_Bcast(varloc, *nbrloc, CS_MPI_REAL, 0, cs_glob_mpi_comm);

  }

#endif /* defined(HAVE_MPI) */

  if (distant == false) {

    nbr = CS_MIN(*nbrdis, *nbrloc);

    for (ind = 0; ind < nbr; ind++)
      varloc[ind] = vardis[ind];

  }
}

/*----------------------------------------------------------------------------
 * Compute the maximum value of an integer variable associated to a coupling.
 *
 * It is assumed that the integer value is the same for each group of
 * processus (local and distant).
 *
 * Fortran interface:
 *
 * SUBROUTINE MXICPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          VALDIS         : --> : distant value (to send)
 * INTEGER          VALMAX         : <-- : local maximum (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (mxicpl, MXICPL)
(
 const cs_int_t  *numcpl,
       cs_int_t  *vardis,
       cs_int_t  *varmax
)
{
  bool  distant = false;

#if defined(HAVE_MPI)

  cs_sat_coupling_t  *coupl = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[*numcpl - 1];

  if (coupl->comm != MPI_COMM_NULL) {

    distant = true;

    MPI_Allreduce(vardis, varmax, 1, CS_MPI_INT, MPI_MAX, coupl->comm);

  }

#endif /* defined(HAVE_MPI) */

  if (distant == false) {

    *varmax = *vardis;

  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define new Code_Saturne coupling.
 *
 * The arguments to \ref cs_sat_coupling_define are:
 * \param[in] saturne_name          matching Code_Saturne application name
 * \param[in] boundary_cpl_criteria boundary face selection criteria for coupled
 *                                  faces, or NULL
 * \param[in] volume_cpl_criteria   cell selection criteria for coupled cells, or
                                    NULL
 * \param[in] boundary_sup_criteria boundary face selection criteria for support
 *                                  (not functional)
 * \param[in] volume_sup_criteria   cell selection criteria for support
 * \param[in] verbosity             verbosity level
 *
 * In the case of only 2 Code_Saturne instances, the 'saturne_name' argument
 * is ignored, as there is only one matching possibility.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * Code_Saturne instances based on the 'saturne_name' argument.
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_define(const char  *saturne_name,
                       const char  *boundary_cpl_criteria,
                       const char  *volume_cpl_criteria,
                       const char  *boundary_sup_criteria,
                       const char  *volume_sup_criteria,
                       int          verbosity)
{
  _cs_sat_coupling_builder_t *scb = NULL;

  /* Add corresponding coupling to temporary Code_Saturne couplings array */

  BFT_REALLOC(_sat_coupling_builder,
              _sat_coupling_builder_size + 1,
              _cs_sat_coupling_builder_t);

  scb = &(_sat_coupling_builder[_sat_coupling_builder_size]);

  scb->match_id = -1;

  scb->app_name = NULL;
  if (saturne_name != NULL) {
    BFT_MALLOC(scb->app_name, strlen(saturne_name) + 1, char);
    strcpy(scb->app_name, saturne_name);
  }

  scb->face_cpl_sel_c = NULL;
  if (boundary_cpl_criteria != NULL) {
    BFT_MALLOC(scb->face_cpl_sel_c, strlen(boundary_cpl_criteria) + 1, char);
    strcpy(scb->face_cpl_sel_c, boundary_cpl_criteria);
  }

  scb->cell_cpl_sel_c = NULL;
  if (volume_cpl_criteria != NULL) {
    BFT_MALLOC(scb->cell_cpl_sel_c, strlen(volume_cpl_criteria) + 1, char);
    strcpy(scb->cell_cpl_sel_c, volume_cpl_criteria);
  }

  scb->face_sup_sel_c = NULL;
  if (boundary_sup_criteria != NULL) {
    BFT_MALLOC(scb->face_sup_sel_c, strlen(boundary_sup_criteria) + 1, char);
    strcpy(scb->face_sup_sel_c, boundary_sup_criteria);
  }

  scb->cell_sup_sel_c = NULL;
  if (volume_sup_criteria != NULL) {
    BFT_MALLOC(scb->cell_sup_sel_c, strlen(volume_sup_criteria) + 1, char);
    strcpy(scb->cell_sup_sel_c, volume_sup_criteria);
  }

  scb->verbosity = verbosity;

  _sat_coupling_builder_size += 1;
}

/*----------------------------------------------------------------------------
 * Get number of Code_Saturne couplings.
 *
 * returns:
 *   number of Code_Saturne couplings
 *----------------------------------------------------------------------------*/

int
cs_sat_coupling_n_couplings(void)
{
  return cs_glob_sat_n_couplings;
}

/*----------------------------------------------------------------------------
 * Get pointer to Code_Saturne coupling.
 *
 * parameters:
 *   coupling_id <-- Id (0 to n-1) of Code_Saturne coupling
 *
 * returns:
 *   pointer to Code_Saturne coupling structure
 *----------------------------------------------------------------------------*/

cs_sat_coupling_t *
cs_sat_coupling_by_id(int coupling_id)
{
  cs_sat_coupling_t  *retval = NULL;

  if (   coupling_id > -1
      && coupling_id < cs_glob_sat_n_couplings)
    retval = cs_glob_sat_couplings[coupling_id];

  return retval;
}

/*----------------------------------------------------------------------------
 * Initialize Code_Saturne couplings.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_all_init(void)
{
  /* First try using MPI */

#if defined(HAVE_MPI)

  if (_sat_coupling_builder_size > 0)
    _init_all_mpi_sat();

#endif

  /* Print unmatched instances */

  if (_sat_coupling_builder_size > 0) {

    bft_printf("Unmatched Code_Saturne couplings:\n"
               "---------------------------------\n\n");

    _print_all_unmatched_sat();

    bft_error(__FILE__, __LINE__, 0,
              _("At least 1 Code_Saturne coupling was defined for which\n"
                "no communication with a Code_Saturne instance is possible."));
  }
}

/*----------------------------------------------------------------------------
 * Create a sat_coupling_t structure.
 *
 * parameters:
 *   ref_axis           <-- reference axis
 *   face_sel_criterion <-- criterion for selection of boundary faces
 *   cell_sel_criterion <-- criterion for selection of cells
 *   sat_name           <-- Code_Saturne application name
 *   verbosity          <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_add(const char  *face_cpl_sel_c,
                    const char  *cell_cpl_sel_c,
                    const char  *face_sup_sel_c,
                    const char  *cell_sup_sel_c,
                    const char  *sat_name,
                    int          verbosity)
{
  cs_sat_coupling_t *sat_coupling = NULL;

  /* Allocate _cs_sat_coupling_t structure */

  BFT_REALLOC(cs_glob_sat_couplings,
              cs_glob_sat_n_couplings + 1, cs_sat_coupling_t*);
  BFT_MALLOC(sat_coupling, 1, cs_sat_coupling_t);

  sat_coupling->sat_name = NULL;

  if (sat_name != NULL) {
    BFT_MALLOC(sat_coupling->sat_name, strlen(sat_name) + 1, char);
    strcpy(sat_coupling->sat_name, sat_name);
  }

  /* Selection criteria  */

  sat_coupling->face_cpl_sel = NULL;
  sat_coupling->cell_cpl_sel = NULL;
  sat_coupling->face_sup_sel = NULL;
  sat_coupling->cell_sup_sel = NULL;

  if (face_cpl_sel_c != NULL) {
    BFT_MALLOC(sat_coupling->face_cpl_sel, strlen(face_cpl_sel_c) + 1, char);
    strcpy(sat_coupling->face_cpl_sel, face_cpl_sel_c);
  }
  if (cell_cpl_sel_c != NULL) {
    BFT_MALLOC(sat_coupling->cell_cpl_sel, strlen(cell_cpl_sel_c) + 1, char);
    strcpy(sat_coupling->cell_cpl_sel, cell_cpl_sel_c);
  }

  if (face_sup_sel_c != NULL) {
    BFT_MALLOC(sat_coupling->face_sup_sel, strlen(face_sup_sel_c) + 1, char);
    strcpy(sat_coupling->face_sup_sel, face_sup_sel_c);
  }
  if (cell_sup_sel_c != NULL) {
    BFT_MALLOC(sat_coupling->cell_sup_sel, strlen(cell_sup_sel_c) + 1, char);
    strcpy(sat_coupling->cell_sup_sel, cell_sup_sel_c);
  }

  sat_coupling->faces_sup = NULL;
  sat_coupling->cells_sup = NULL;

  sat_coupling->localis_fbr = NULL;
  sat_coupling->localis_cel = NULL;

  sat_coupling->nbr_fbr_sup = 0;
  sat_coupling->nbr_cel_sup = 0;

  sat_coupling->verbosity = verbosity;

  /* Geometric quantities arrays for interpolation */

  sat_coupling->distant_dist_fbr = NULL;
  sat_coupling->distant_of = NULL;
  sat_coupling->local_of = NULL;
  sat_coupling->distant_pond_fbr = NULL;
  sat_coupling->local_pond_fbr = NULL;

  /* Initialize communicators */

#if defined(HAVE_MPI)

  sat_coupling->comm = MPI_COMM_NULL;
  sat_coupling->n_sat_ranks = 0;
  sat_coupling->sat_root_rank = -1;

#endif

  /* Update coupling array and return */

  cs_glob_sat_couplings[cs_glob_sat_n_couplings] = sat_coupling;
  cs_glob_sat_n_couplings++;
}

/*----------------------------------------------------------------------------
 * Destroy all couplings
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_all_finalize(void)
{
  int  i;

  for (i = 0 ; i < cs_glob_sat_n_couplings ; i++)
    _sat_coupling_destroy(cs_glob_sat_couplings[i]);

  BFT_FREE(cs_glob_sat_couplings);

  cs_glob_sat_n_couplings = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
