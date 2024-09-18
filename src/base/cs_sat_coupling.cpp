/*============================================================================
 * Functions associated with code coupling.
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

#include "cs_ale.h"
#include "cs_assert.h"
#include "cs_array.h"
#include "cs_atmo.h"
#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_coupling.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_connect.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_selector.h"
#include "cs_turbulence_model.h"
#include "cs_turbomachinery.h"

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

/* Structure associated with code_saturne coupling */

typedef struct {

  int      match_id;        /* Id of matched application, -1 initially */
  char    *app_name;        /* Application name */
  char    *face_cpl_sel_c;  /* Face selection criteria */
  char    *cell_cpl_sel_c;  /* Cell selection criteria */
  char    *face_loc_sel_c;  /* Face selection criteria */
  char    *cell_loc_sel_c;  /* Cell selection criteria */
  int      verbosity;       /* Verbosity level */
  int      reverse;         /* Reverse mode if 1 */

} _cs_sat_coupling_builder_t;

struct _cs_sat_coupling_t {

  char                   *sat_name;        /* Application name */
  cs_sat_coupling_tag_t  *tag_func;        /* Tagging function pointer */
  void                   *tag_context;     /* Tagging context */
  int                    reverse;          /* Reverse mode */

  char            *face_cpl_sel; /* Local face overlapped selection criteria */
  char            *cell_cpl_sel; /* Local cell overlapped selection criteria */
  char            *face_loc_sel; /* Distant face overlapped selection criteria */
  char            *cell_loc_sel; /* Distant cell overlapped selection criteria */

  ple_locator_t   *localis_cel;  /* Locator associated with cells */
  ple_locator_t   *localis_fbr;  /* Locator associated with boundary faces */

  cs_lnum_t        nbr_cel_sup;  /* Number of associated cell locations */
  cs_lnum_t        nbr_fbr_sup;  /* Number of associated face locations */

  fvm_nodal_t     *cells_sup;    /* Local cells at which distant values are
                                    interpolated*/
  fvm_nodal_t     *faces_sup;    /* Local faces at which distant values are
                                    interpolated*/

  cs_real_t       *distant_dist_fbr; /* Distant vectors (distance JJ') */
  cs_real_t       *distant_of;       /* Distant vector OF */
  cs_real_t       *local_of;         /* Local vector OF */
  cs_real_t       *distant_pond_fbr; /* Distant weighting coefficient */
  cs_real_t       *local_pond_fbr;   /* Local weighting coefficient */

  cs_real_t        tolerance; /* location tolerance */
  int              verbosity; /* Verbosity level */

  /* variables used for coherency checks over coupling */
  int              icorio; /* Reference frame of resolution */
  int              ale;    /* ALE model */
  int              imajcp; /* ALE/Deformable mesh */
  int              nvarcp; /* Local number of solved variables */
  int              nvarto; /* Max number of solved variables in coupling */
  int              iturb;  /* Global turbulence model */

  /* Communication-related members */

#if defined(HAVE_MPI)

  MPI_Comm         comm;           /* Associated MPI communicator */

  int              n_sat_ranks;    /* Number of associated code_saturne ranks */
  int              sat_root_rank;  /* First associated code_saturne rank */

#endif

};

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Array of couplings */

static int                         _sat_coupling_builder_size = 0;
static _cs_sat_coupling_builder_t *_sat_coupling_builder = nullptr;

static int                  cs_glob_sat_n_couplings = 0;
static cs_sat_coupling_t  **cs_glob_sat_couplings = nullptr;

static int _cs_sat_coupling_initialized = 0;

/*============================================================================
 * Global variables
 *============================================================================*/

int  cs_glob_sat_coupling_face_interpolation_type = 0;

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

    if (scb->match_id != -1) {
      if (scb->face_cpl_sel_c != nullptr) BFT_FREE(scb->face_cpl_sel_c);
      if (scb->cell_cpl_sel_c != nullptr) BFT_FREE(scb->cell_cpl_sel_c);
      if (scb->face_loc_sel_c != nullptr) BFT_FREE(scb->face_loc_sel_c);
      if (scb->cell_loc_sel_c != nullptr) BFT_FREE(scb->cell_loc_sel_c);
      if (scb->app_name != nullptr) BFT_FREE(scb->app_name);
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
 * Print information on yet unmatched code_saturne couplings.
 *----------------------------------------------------------------------------*/

static void
_print_all_unmatched_sat(void)
{
  const char empty_string[] = "";

  /* Loop on defined code_saturne instances */

  for (int i = 0; i < _sat_coupling_builder_size; i++) {

    _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + i;

    if (scb->match_id < 0) {

      const char *local_name = empty_string;

      if (scb->app_name != nullptr)
        local_name = scb->app_name;

      bft_printf(_(" code_saturne coupling:\n"
                   "   coupling id:              %d\n"
                   "   local name:               \"%s\"\n\n"),
                 i, local_name);
    }
  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Initialize communicator for code_saturne coupling
 *
 * parameters:
 *   sat_coupling  <-> code_saturne coupling structure
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

  bft_printf(_(" code_saturne coupling %d: initializing MPI communication ... "),
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
 * Add a code_saturne coupling using MPI.
 *
 * parameters:
 *   builder_id    <-- code_saturne application id in coupling builder
 *   sat_root_rank <-- root rank associated with code_saturne
 *   n_sat_ranks   <-- number of ranks associated with code_saturne
 *----------------------------------------------------------------------------*/

static void
_sat_add_mpi(int builder_id,
             int sat_root_rank,
             int n_sat_ranks)
{
  cs_sat_coupling_t *sat_coupling = nullptr;
  _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + builder_id;

  /* Similarly to SYRTHES, we might be able to add
     code_saturne couplings directly (without resorting
     to a temporary builder), then match communications */

  cs_sat_coupling_add(scb->face_cpl_sel_c,
                      scb->cell_cpl_sel_c,
                      scb->face_loc_sel_c,
                      scb->cell_loc_sel_c,
                      scb->app_name,
                      scb->reverse,
                      scb->verbosity);

  sat_coupling = cs_sat_coupling_by_id(cs_sat_coupling_n_couplings() - 1);

  sat_coupling->sat_root_rank = sat_root_rank;
  sat_coupling->n_sat_ranks = n_sat_ranks;

  _init_comm(sat_coupling,
             builder_id);
}

/*----------------------------------------------------------------------------
 * Print information on identified code_saturne couplings using MPI.
 *
 * This function requires coupling_builder information, and must thus
 * be called before removing matched builder entries.
 *----------------------------------------------------------------------------*/

static void
_print_all_mpi_sat(void)
{
  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();
  const char empty_string[] = "";

  /* Loop on defined code_saturne instances */

  for (int i = 0; i < _sat_coupling_builder_size; i++) {

    _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + i;

    if (scb->match_id > -1) {

      const char *local_name = empty_string;
      const char *distant_name = empty_string;

      const ple_coupling_mpi_set_info_t
        ai = ple_coupling_mpi_set_get_info(mpi_apps, scb->match_id);

      if (scb->app_name != nullptr)
        local_name = scb->app_name;
      if (ai.app_name != nullptr)
        distant_name = ai.app_name;

      bft_printf(_(" code_saturne coupling:\n"
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
 * Initialize MPI code_saturne couplings using MPI.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

static void
_init_all_mpi_sat(void)
{
  int n_apps = 0;
  int n_sat_apps = 0;

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == nullptr)
    return;

  n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);

  /* First pass to count available code_saturne couplings */

  for (int i = 0; i < n_apps; i++) {
    const ple_coupling_mpi_set_info_t
      ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
    if (   strncmp(ai.app_type, "code_saturne", 12) == 0
        || strncmp(ai.app_type, "Code_Saturne", 12) == 0)
      n_sat_apps += 1;
  }

  /* In single-coupling mode, no identification necessary */

  if (n_sat_apps == 2 && _sat_coupling_builder_size == 1) {

    const int local_app_id = ple_coupling_mpi_set_get_app_id(mpi_apps);

    for (int i = 0; i < n_apps; i++) {
      const ple_coupling_mpi_set_info_t
        ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
      if (   (   strncmp(ai.app_type, "code_saturne", 12) == 0
              || strncmp(ai.app_type, "Code_Saturne", 12) == 0)
          && i != local_app_id) {
        _sat_coupling_builder->match_id = i;
      }
    }

  }

  /* In multiple-coupling mode, identification is necessary */

  else {

    ple_coupling_mpi_set_info_t ai;

    int *sat_appinfo = nullptr;

    /* First, build an array of matched/unmatched code_saturne applications,
       with 2 entries per instance: matched indicator, app_id */

    BFT_MALLOC(sat_appinfo, n_sat_apps*2, int);

    n_sat_apps = 0;

    for (int i = 0; i < n_apps; i++) {
      ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
      if (   strncmp(ai.app_type, "code_saturne", 12) == 0
          || strncmp(ai.app_type, "Code_Saturne", 12) == 0) {
        sat_appinfo[n_sat_apps*2] = 0;
        sat_appinfo[n_sat_apps*2 + 1] = i;
        n_sat_apps += 1;
      }
    }

    /* Loop on defined code_saturne instances */

    for (int i = 0; i < _sat_coupling_builder_size; i++) {

      _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + i;

      /* Loop on available code_saturne instances to match app_names */

      if (scb->app_name != nullptr) {

        for (int j = 0; j < n_sat_apps; j++) {

          if (sat_appinfo[j*2] != 0) /* Consider only unmatched applications */
            continue;

          ai = ple_coupling_mpi_set_get_info(mpi_apps, sat_appinfo[j*2 + 1]);
          if (ai.app_name == nullptr)
            continue;

          if (strcmp(ai.app_name, scb->app_name) == 0) {
            scb->match_id = sat_appinfo[j*2 + 1];
            sat_appinfo[j*2] = i;
            break;
          }

        }

      }

    } /* End of loop on defined code_saturne instances */

    BFT_FREE(sat_appinfo);

  } /* End of test on single or multiple code_saturne matching algorithm */

  /* Print matching info */

  _print_all_mpi_sat();

  /* Now initialize matched couplings */
  /*----------------------------------*/

  for (int i = 0; i < _sat_coupling_builder_size; i++) {

    _cs_sat_coupling_builder_t *scb = _sat_coupling_builder + i;

    if (scb->match_id > -1) {
      const ple_coupling_mpi_set_info_t
        ai = ple_coupling_mpi_set_get_info(mpi_apps, scb->match_id);

      if (   strncmp(ai.app_type, "code_saturne", 12) == 0
          || strncmp(ai.app_type, "Code_Saturne", 12) == 0)
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
 *   cpl <-> pointer to coupling structure to destroy
 *
 * returns:
 *   nullptr pointer
 *----------------------------------------------------------------------------*/

static cs_sat_coupling_t *
_sat_coupling_destroy(cs_sat_coupling_t  *cpl)
{
  BFT_FREE(cpl->sat_name);

  BFT_FREE(cpl->face_cpl_sel);
  BFT_FREE(cpl->cell_cpl_sel);
  BFT_FREE(cpl->face_loc_sel);
  BFT_FREE(cpl->cell_loc_sel);

  ple_locator_destroy(cpl->localis_cel);
  ple_locator_destroy(cpl->localis_fbr);

  if (cpl->cells_sup != nullptr)
    fvm_nodal_destroy(cpl->cells_sup);
  if (cpl->faces_sup != nullptr)
    fvm_nodal_destroy(cpl->faces_sup);

  BFT_FREE(cpl->distant_dist_fbr);
  BFT_FREE(cpl->distant_of);
  BFT_FREE(cpl->local_of);
  BFT_FREE(cpl->distant_pond_fbr);
  BFT_FREE(cpl->local_pond_fbr);

#if defined(HAVE_MPI)
  if (   cpl->comm != MPI_COMM_WORLD
      && cpl->comm != cs_glob_mpi_comm)
    MPI_Comm_free(&(cpl->comm));
#endif

  BFT_FREE(cpl);

  return nullptr;
}

/*----------------------------------------------------------------------------
 * Computed some quantities needed for a centered-like interpolation
 *  - distance JJ' for distant boundary faces
 *  - local weighting coefficients
 *----------------------------------------------------------------------------*/

static void
_sat_coupling_interpolate(cs_sat_coupling_t  *cpl)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;
  const cs_real_3_t *b_face_cog = (const cs_real_3_t *)mq->b_face_cog;
  const cs_real_3_t *b_face_u_normal = (const cs_real_3_t *)mq->b_face_u_normal;

  /* Removing the connectivity and localization informations in case of
     coupling update */

  if (cpl->distant_dist_fbr != nullptr)
    BFT_FREE(cpl->distant_dist_fbr);
  if (cpl->distant_of != nullptr)
    BFT_FREE(cpl->distant_of);
  if (cpl->local_of != nullptr)
    BFT_FREE(cpl->local_of);
  if (cpl->distant_pond_fbr != nullptr)
    BFT_FREE(cpl->distant_pond_fbr);
  if (cpl->local_pond_fbr != nullptr)
    BFT_FREE(cpl->local_pond_fbr);

  /* Interpolation structure */

  const cs_lnum_t n_fbr_loc = ple_locator_get_n_interior(cpl->localis_fbr);
  const cs_lnum_t *lstfbr = ple_locator_get_interior_list(cpl->localis_fbr);

  const cs_lnum_t n_fbr_dist
    = ple_locator_get_n_dist_points(cpl->localis_fbr);
  const cs_lnum_t *element
    = ple_locator_get_dist_locations(cpl->localis_fbr);
  const cs_coord_t *distant_coord
    = ple_locator_get_dist_coords(cpl->localis_fbr);

  /* Calculation of the distance DJJPB defining the distance from */
  /* the local cell center to the distant boundary face norm      */
  /*--------------------------------------------------------------*/

  BFT_MALLOC(cpl->distant_dist_fbr, 3*n_fbr_dist, cs_real_t);

  /* Store the local surface vector of the coupled boundary faces */

  cs_real_t *local_u_norm;
  BFT_MALLOC(local_u_norm, 3*n_fbr_loc, cs_real_t);

  for (cs_lnum_t ind = 0; ind < n_fbr_loc; ind++) {
    cs_lnum_t ifac = lstfbr[ind];

    for (cs_lnum_t icoo = 0; icoo < 3; icoo++)
      local_u_norm[ind*3 + icoo] = b_face_u_normal[ifac][icoo];
  }

  /* Get the distant faces surface vector (reverse = 1) */

  cs_real_t *distant_u_norm;
  BFT_MALLOC(distant_u_norm, 3*n_fbr_dist, cs_real_t);

  ple_locator_exchange_point_var(cpl->localis_fbr,
                                 distant_u_norm,
                                 local_u_norm,
                                 nullptr,
                                 sizeof(cs_real_t),
                                 3,
                                 1); /* reverse */

  BFT_FREE(local_u_norm);

  /* Calculation of the JJ' vectors */

  cs_real_t *distant_xyzcen;
  BFT_MALLOC(distant_xyzcen, 3*n_fbr_dist, cs_real_t);

  for (cs_lnum_t ind = 0; ind < n_fbr_dist; ind++) {

    cs_lnum_t c_id = element[ind];

    cs_real_t pdt_scal = 0.;
    cs_real_t dist_cel_fbr[3];
    const cs_real_t *vect_surf_norm = distant_u_norm + ind*3;

    for (cs_lnum_t icoo = 0; icoo < 3; icoo++) {
      dist_cel_fbr[icoo] = distant_coord[ind*3 + icoo] - cell_cen[c_id][icoo];

      /* Store the distant coordinates to compute the weighting coefficients */
      distant_xyzcen[ind*3 + icoo] = cell_cen[c_id][icoo];
      pdt_scal += dist_cel_fbr[icoo]*vect_surf_norm[icoo];
    }

    for (cs_lnum_t icoo = 0; icoo < 3; icoo++)
      cpl->distant_dist_fbr[ind*3 + icoo]
        = dist_cel_fbr[icoo] - pdt_scal*vect_surf_norm[icoo];

  }

  BFT_FREE(distant_u_norm);

  /* Calculation of the local weighting coefficient */
  /*------------------------------------------------*/

  BFT_MALLOC(cpl->distant_pond_fbr, n_fbr_dist, cs_real_t);
  BFT_MALLOC(cpl->local_pond_fbr, n_fbr_loc, cs_real_t);

  /* Get the cell center coordinates (reverse = 0) */

  cs_real_t *local_xyzcen;
  BFT_MALLOC(local_xyzcen, 3*n_fbr_loc, cs_real_t);

  ple_locator_exchange_point_var(cpl->localis_fbr,
                                 distant_xyzcen,
                                 local_xyzcen,
                                 nullptr,
                                 sizeof(cs_real_t),
                                 3,
                                 0);

  BFT_FREE(distant_xyzcen);

  /* Calculation of the local weighting coefficients */

  for (cs_lnum_t ind = 0; ind < n_fbr_loc; ind++) {

    cs_lnum_t ifac = lstfbr[ind];
    cs_lnum_t c_id  = mesh->b_face_cells[ifac];

    cs_real_t distance_fbr_cel = 0.;
    cs_real_t distance_cel_cel = 0.;

    for (cs_lnum_t icoo = 0; icoo < 3; icoo++) {
      distance_fbr_cel
        +=   b_face_u_normal[ifac][icoo]
           * (local_xyzcen[ind*3 + icoo] - b_face_cog[ifac][icoo]);

      distance_cel_cel
        +=   b_face_u_normal[ifac][icoo]
           * (local_xyzcen[ind*3 + icoo] - cell_cen[c_id][icoo]);
    }

    if (fabs(distance_cel_cel) > 1.e-12)
      cpl->local_pond_fbr[ind] = distance_fbr_cel / distance_cel_cel;
    else
      cpl->local_pond_fbr[ind] = 0.5;

  }

  /* Get the distant weighting coefficients (reverse = 1) */

  ple_locator_exchange_point_var(cpl->localis_fbr,
                                 cpl->distant_pond_fbr,
                                 cpl->local_pond_fbr,
                                 nullptr,
                                 sizeof(cs_real_t),
                                 1,
                                 1); /* reverse */

  /* Calculation of the OF distance */
  /*--------------------------------*/

  BFT_MALLOC(cpl->distant_of, 3*n_fbr_dist, cs_real_t);
  BFT_MALLOC(cpl->local_of, 3*n_fbr_loc, cs_real_t);

  for (cs_lnum_t ind = 0; ind < n_fbr_loc; ind++) {

    cs_lnum_t ifac = lstfbr[ind];
    cs_lnum_t c_id  = mesh->b_face_cells[ifac];

    /* n/norm(n) . IF = FJ' */
    cs_real_t distance_fbr_cel = 0.;
    /* n/norm(n) . IJ = I'J'*/
    cs_real_t distance_cel_cel = 0.;

    for (cs_lnum_t icoo = 0; icoo < 3; icoo++) {

      distance_fbr_cel +=   b_face_u_normal[ifac][icoo]
                          * (local_xyzcen[ind*3+icoo] - b_face_cog[ifac][icoo]);

      distance_cel_cel +=   b_face_u_normal[ifac][icoo]
                          * (local_xyzcen[ind*3+icoo] - cell_cen[c_id][icoo]);

    }

    for (cs_lnum_t icoo = 0; icoo < 3; icoo++) {
      cpl->local_of[ind*3 + icoo]
        =     b_face_cog[ifac][icoo]
          -  (      b_face_cog[ifac][icoo]        /*  O'  */
              +       b_face_u_normal[ifac][icoo]
                    * distance_fbr_cel            /*J'=F+n*FJ'*/
              - 0.5 * b_face_u_normal[ifac][icoo]
                    * distance_cel_cel);          /*-n*I'J'/2*/
    }
  }

  ple_locator_exchange_point_var(cpl->localis_fbr,
                                 cpl->distant_of,
                                 cpl->local_of,
                                 nullptr,
                                 sizeof(cs_real_t),
                                 3,
                                 1);  /* reverse */

  BFT_FREE(local_xyzcen);
}

/*--------------------------------------------------------------------------- */
/*!
 * \brief Compute volumic coupling terms.
 */
/*--------------------------------------------------------------------------- */

static void
_sat_coupling_compute_data_at_cells
(
  cs_field_t       *f,       /*!<[in] pointer to cs_field_t */
  const int         reverse, /*!<[in] direction of coupling */
  const cs_lnum_t   n_elts,  /*!<[in] number of coupled cells  */
  const cs_lnum_t  *elt_ids, /*!<[in] list of coupled cells ids */
  const cs_coord_t *coords,  /*!<[in] coordinates of distant elements */
  cs_real_t        *warray1, /*!<[out] work array of explicit source terms */
  cs_real_t        *warray2  /*!<[out] work array of implicit source terms */
)
{
  /* For scalars */
  if (reverse == 0) {
    if (f->dim == 1) {
      cs_real_3_t *grad = nullptr;
      BFT_MALLOC(grad, cs_glob_mesh->n_cells_with_ghosts, cs_real_3_t);
      cs_field_gradient_scalar(f,
                               true, /* use previous */
                               1,    /* inc */
                               grad);

      /* Interpolation */
      cs_real_t *cell_cen = cs_glob_mesh_quantities->cell_cen;
      for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
        cs_lnum_t c_id = elt_ids[e_id];
        cs_real_t dxyz[3] = {0.};
        for (int i = 0; i < 3; i++)
          dxyz[i] = coords[3*e_id + i] - cell_cen[3*c_id + i];

        warray1[e_id] = f->val_pre[c_id]
                      + cs_math_3_dot_product(dxyz, grad[c_id]);

      }
      BFT_FREE(grad);
    }
    else if (f->dim == 3) {
      cs_real_33_t *grad = nullptr;
      BFT_MALLOC(grad, cs_glob_mesh->n_cells_with_ghosts, cs_real_33_t);
      cs_field_gradient_vector(f,
                               true, /* use previous */
                               1,    /* inc */
                               grad);

      /* Interpolation */
      cs_real_t *cell_cen = cs_glob_mesh_quantities->cell_cen;
      for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
        cs_lnum_t c_id = elt_ids[e_id];
        cs_real_t dxyz[3] = {0.};
        for (int i = 0; i < 3; i++)
          dxyz[i] = coords[3*e_id + i] - cell_cen[3*c_id + i];

        for (int j = 0; j < f->dim; j++)
          warray1[3 * e_id + j] = f->val_pre[3 * c_id + j]
                                + cs_math_3_dot_product(dxyz, grad[c_id][j]);
      }
      BFT_FREE(grad);
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("%s: Exchange is possible only for scalars and vectors.\n"),
                __func__);
  }
  else {
    cs_real_t *cvar_rho = CS_F_(rho)->val;
    cs_real_t *fluid_vol = cs_glob_mesh_quantities->cell_f_vol;
    if (f->dim == 1) {
      for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
        cs_lnum_t c_id = elt_ids[e_id];

        /* Mass weighted averaged value */
        warray1[e_id] = fluid_vol[c_id] * cvar_rho[c_id] * f->val_pre[c_id];

        /* Mass */
        warray2[e_id] = fluid_vol[c_id] * cvar_rho[c_id];
      }
    }
    else if (f->dim == 3) {
      int dim2 = f->dim * f->dim;
      for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
        cs_lnum_t c_id = elt_ids[e_id];

        for (int i = 0; i < f->dim; i++) {
          /* Mass weighted averaged value */
          warray1[3 * e_id + i] = fluid_vol[c_id] * cvar_rho[c_id]
                                * f->val_pre[3 * c_id + i];
          /* implicit part */
          for (int j = 0; j < f->dim; j++) {
            if (i == j)
              warray2[dim2 * e_id + 3 * i + j] =   fluid_vol[c_id]
                                                 * cvar_rho[c_id];
            else
              warray2[dim2 * e_id + 3 * i + j] = 0.;
          }
        }
      }
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("%s: Exchange is possible only for scalars and vectors.\n"),
                __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the maximum value of an integer variable associated
 *        to a coupling.
 *
 * It is assumed that the integer value is the same for each group of
 * processus (local and distant).
 */
/*----------------------------------------------------------------------------*/

static void
_sat_coupling_int_max
(
 cs_sat_coupling_t *coupl,
 const cs_lnum_t   *vardis,
 cs_lnum_t         *varmax
)
{
  bool  distant = false;

#if defined(HAVE_MPI)

  /* Initializations and verifications */

  if (coupl->comm != MPI_COMM_NULL) {

    distant = true;

    MPI_Allreduce(vardis, varmax, 1, CS_MPI_LNUM, MPI_MAX, coupl->comm);

  }

#endif /* defined(HAVE_MPI) */

  if (distant == false) {

    *varmax = *vardis;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Array of integers exchange, associated to a given coupling.
 *
 * It is assumed that the arrays have the same size and the same values on
 * each group of processus (local and distant).
 */
/*----------------------------------------------------------------------------*/

static void
_sat_coupling_check_turbulence_coherency
(
  cs_sat_coupling_t *coupl, /*!<[in] pointer to cs_sat_coupling_t */
  int                iturb  /*!<[in] local turbulence model */
)
{
  bool  distant = false;

#if defined(HAVE_MPI)

  MPI_Status  status;

  /* Initializations and verifications */

  if (coupl->comm != MPI_COMM_NULL) {

    distant = true;

    /* Exchange between the groups master node */

    if (cs_glob_rank_id < 1)
      MPI_Sendrecv(&iturb, 1, CS_MPI_LNUM, coupl->sat_root_rank, 0,
                   &(coupl->iturb), 1, CS_MPI_LNUM, coupl->sat_root_rank, 0,
                   coupl->comm, &status);

    /* Synchronization inside a group */

    if (cs_glob_n_ranks > 1)
      MPI_Bcast (&(coupl->iturb), 1, CS_MPI_LNUM, 0, cs_glob_mpi_comm);

  }

#endif /* defined(HAVE_MPI) */

  if (distant == false) {

    coupl->iturb = iturb;

  }

  /* Check coherency, since only RANS and/or laminar models are handled. */

  if (iturb == CS_TURB_V2F_PHI && coupl->iturb != CS_TURB_V2F_PHI) {
    bft_error(__FILE__, __LINE__, 0,
              _("V2F PHI FBAR can only be coupled with itself.\n"));
  }
  else if (iturb == CS_TURB_V2F_BL_V2K && coupl->iturb != CS_TURB_V2F_BL_V2K) {
    bft_error(__FILE__, __LINE__, 0,
              _("V2F BL-V2/K can only be coupled with itself.\n"));
  }
  else if (  (iturb >= CS_TURB_LES_SMAGO_CONST && iturb <= CS_TURB_LES_WALE)
           && !(   coupl->iturb >= CS_TURB_LES_SMAGO_CONST
                && coupl->iturb <= CS_TURB_LES_WALE)) {
    bft_error(__FILE__, __LINE__, 0,
              _("LES/RANS coupling is not yet handled.\n"));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Scalar gradient base projection.
 */
/*----------------------------------------------------------------------------*/

static void
_sat_coupling_interpolate_scalar_simple
(
  cs_real_t   *vals,
  const cs_real_3_t *grad,
  const cs_lnum_t    n_elts,
  const cs_lnum_t   *elt_ids,
  const cs_real_t   *distant_dist_fbr,
  cs_real_t   *interp_vals
)
{
  /* Sanity checks */
  assert(vals != nullptr);
  assert(grad != nullptr);

  /* Interpolate scalar using simple projection */
  for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
    cs_lnum_t c_id = elt_ids[e_id];
    const cs_real_t *xyzjjp = distant_dist_fbr + e_id * 3;

    interp_vals[e_id] = vals[c_id]
                      + cs_math_3_dot_product(xyzjjp, grad[c_id]);
  }
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Send varialbes data related to surface coupling.
 */
/*----------------------------------------------------------------------------*/

static void
_sat_coupling_send_bnd_data
(
  cs_sat_coupling_t  *coupl,            /*!<[in] pointer to coupling structure */
  const cs_coord_t   *coords,           /*!<[in] pointer to coordinates */
  const cs_lnum_t     n_b_faces_dist,   /*!<[in] number of coupled faces */
  const cs_lnum_t    *b_faces_dist_ids, /*!<[in] array of coupled faces ids */
  const cs_real_t    *distant_dist_fbr, /*!<[in] distance array */
  cs_real_t          *distant_pond_fbr, /*!<[in] ponderation array*/
  cs_real_t          *distant_of,       /*!<[in] OF array */
  cs_real_t         **rvdis             /*!<[in] data array */
)
{
  /* Pointers to field values */
  cs_real_t *cvar_vel = CS_F_(vel)->val;
  cs_real_t *cvar_p = CS_F_(p)->val;

  /* Global values used throughout this function */
  const int iprev = 0;
  const int inc   = 1;
  const int itytur_cpl = coupl->iturb / 10;

  cs_real_3_t *grad = nullptr;
  BFT_MALLOC(grad, cs_glob_mesh->n_cells_with_ghosts, cs_real_3_t);

  cs_real_33_t *gradv = nullptr;
  BFT_MALLOC(gradv, cs_glob_mesh->n_cells_with_ghosts, cs_real_33_t);

  cs_real_t *trav1 = nullptr;
  BFT_MALLOC(trav1, n_b_faces_dist, cs_real_t);

  cs_real_t *trav2 = nullptr;
  BFT_MALLOC(trav2, n_b_faces_dist, cs_real_t);

  cs_real_t *trav3 = nullptr;
  BFT_MALLOC(trav3, n_b_faces_dist, cs_real_t);

  cs_real_t *trav4 = nullptr;
  BFT_MALLOC(trav4, n_b_faces_dist, cs_real_t);

  cs_real_t *trav5 = nullptr;
  BFT_MALLOC(trav5, n_b_faces_dist, cs_real_t);

  cs_real_t *trav6 = nullptr;
  BFT_MALLOC(trav6, n_b_faces_dist, cs_real_t);

  cs_real_t *trav7 = nullptr;
  BFT_MALLOC(trav7, n_b_faces_dist, cs_real_t);

  cs_real_t *trav8 = nullptr;
  BFT_MALLOC(trav8, n_b_faces_dist, cs_real_t);

  int ipos = 0;

  /* Velocity
     -------- */

  cs_field_gradient_vector(CS_F_(vel), iprev, inc, gradv);

  if (cs_glob_sat_coupling_face_interpolation_type == 1) {

    /* For the velocity we want to impose a velocity Dirichlet which mimics
     * the behavior of an interior face. We can choose between
     * UPWIND, SOLU et CENTERED (commentet parts based on choice).
     * For now, only CENTERED behaves the same a a single domain with
     * respecte to diffusion.
     *
     * -- UPWIND
     *
     *        xjjp = djppts[ipt][0]
     *        yjjp = djppts[ipt][1]
     *        zjjp = djppts[ipt][2]
     *
     *        rvdis[ipos][ipt] = vel[c_id][isou]
     *
     * -- SOLU
     *
     *        xjf = coopts[ipt][0] - xyzcen[c_id][0]
     *        yjf = coopts[ipt][1] - xyzcen[c_id][1]
     *        zjf = coopts[ipt][2] - xyzcen[c_id][2]
     *
     *        rvdis[ipos][ipt] =   vel[c_id][isou]
     *                           + xjf*gradv[c_id][isou][0]
     *                           + yjf*gradv[c_id][isou][1]
     *                           + zjf*gradv[c_id][isou][2]
     *
     * -- CENTRE
     */

    for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
      cs_lnum_t c_id = b_faces_dist_ids[e_id];
      for (cs_lnum_t i = 0; i < 3; i++)
        rvdis[ipos+i][e_id]
          =   cvar_vel[3*c_id + i]
            + cs_math_3_dot_product(distant_dist_fbr + 3*e_id,
                                    gradv[c_id][i]);
    }
  }
  else {
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        cs_lnum_t c_id = b_faces_dist_ids[e_id];

        cs_real_t xyzjjp[3] = {0.};
        for (int j = 0; j < 3; j++)
          xyzjjp[j] = distant_of[3*e_id + j] + distant_dist_fbr[3*e_id + j];

        rvdis[ipos+i][e_id]
          =   cvar_vel[3*c_id + i]
            + cs_math_3_dot_product(xyzjjp, gradv[c_id][i]);
      }
    }
  }

  ipos += 3;

  /* -------- */
  /* Pressure */
  /* -------- */

  cs_field_gradient_scalar(CS_F_(p), iprev, inc, grad);
  if (cs_glob_sat_coupling_face_interpolation_type == 1) {
    cs_real_t *_rvdis = rvdis[ipos];

    /*
     * For the pressure we want to impose a Dirichlet such that the
     * pressure gradient is conserved between the 2 coupled domains.
     * For this we use a centered interpolation.
     */

    _sat_coupling_interpolate_scalar_simple(cvar_p, grad,
                                            n_b_faces_dist, b_faces_dist_ids,
                                            distant_dist_fbr, _rvdis);
    /* For a generic coupling no assumption can be made. */
  }
  else {
    cs_real_t *_rvdis = rvdis[ipos];

    cs_real_t *cell_cen = cs_glob_mesh_quantities->cell_cen;

    for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
      cs_lnum_t c_id = b_faces_dist_ids[e_id];

      cs_real_t xyzjpf[3] = {0.};
      for (int i = 0; i < 3; i++)
        xyzjpf[i] =   coords[3 * e_id + i]
                    - cell_cen[3 * c_id + i]
                    - distant_dist_fbr[3 * e_id + i];

      cs_real_t one_ov_jpf = 1./cs_math_3_norm(xyzjpf);

      if (distant_pond_fbr[e_id] >= 0. && distant_pond_fbr[e_id] <= 1.)
        one_ov_jpf *= -1.;

      _rvdis[e_id] = cs_math_3_dot_product(xyzjpf, grad[c_id]) * one_ov_jpf;
    }
  }

  ipos += 1;

  /* Turbulence
     ---------- */

  if (cs_glob_turb_model->itytur == 2) {

    /* k-epsilon
       --------- */

    /* Interpolate local values */

    /* k */
    cs_field_gradient_scalar(CS_F_(k), iprev, inc, grad);
    _sat_coupling_interpolate_scalar_simple(CS_F_(k)->val, grad,
                                            n_b_faces_dist, b_faces_dist_ids,
                                            distant_dist_fbr, trav1);

    /* epsilon */
    cs_field_gradient_scalar(CS_F_(eps), iprev, inc, grad);
    _sat_coupling_interpolate_scalar_simple(CS_F_(eps)->val, grad,
                                            n_b_faces_dist, b_faces_dist_ids,
                                            distant_dist_fbr, trav2);

    /* Transfer values depending on receiving model */

    if (itytur_cpl == 2) {
      /* Copy k and epsilon as is */
      cs_real_t *_k = rvdis[ipos];
      cs_real_t *_eps = rvdis[ipos+1];

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        _k[e_id] = trav1[e_id];
        _eps[e_id] = trav2[e_id];
      }

      ipos += 2;
    }
    else if (itytur_cpl == 3) {
      /* k-eps => Rij-eps */
      const cs_real_t two_ov_3 = 2. / 3.;

      /* R11, R22 and R33 */
      cs_real_t *_r11 = rvdis[ipos];
      cs_real_t *_r22 = rvdis[ipos + 1];
      cs_real_t *_r33 = rvdis[ipos + 2];

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        _r11[e_id] = two_ov_3 * trav1[e_id];
        _r22[e_id] = two_ov_3 * trav1[e_id];
        _r33[e_id] = two_ov_3 * trav1[e_id];
      }

      ipos += 3;

      /* R12, R13 and R23 using previously computed velocity gradient  */

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        cs_lnum_t c_id = b_faces_dist_ids[e_id];

        cs_real_3_t *_grad_uvw = gradv[c_id];

        trav3[e_id] = _grad_uvw[0][1];
        trav4[e_id] = _grad_uvw[0][2];

        trav5[e_id] = _grad_uvw[1][0];
        trav6[e_id] = _grad_uvw[1][2];

        trav7[e_id] = _grad_uvw[2][0];
        trav8[e_id] = _grad_uvw[2][1];
      }

      cs_real_t *_r12 = rvdis[ipos];
      cs_real_t *_r13 = rvdis[ipos + 1];
      cs_real_t *_r23 = rvdis[ipos + 2];
      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        cs_real_t r0 =   -1.0 * trav1[e_id] * trav1[e_id] * cs_turb_cmu
                       / cs_math_fmax(1.e-10, trav2[e_id]);

        _r12[e_id] = r0 * (trav3[e_id] + trav5[e_id]);
        _r13[e_id] = r0 * (trav4[e_id] + trav7[e_id]);
        _r23[e_id] = r0 * (trav6[e_id] + trav8[e_id]);
      }
      ipos += 3;

      /* epsilon */
      cs_real_t *_eps = rvdis[ipos];
      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        _eps[e_id] = trav2[e_id];
      }

      ipos += 1;
    }
    else if (coupl->iturb == CS_TURB_V2F_PHI) {
      /* Option is unavailable...*/
      cs_assert(0);
    }
    else if (coupl->iturb == CS_TURB_K_OMEGA) {
      /* k-eps => k-omega */
      cs_real_t *_k = rvdis[ipos];
      cs_real_t *_om = rvdis[ipos + 1];

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        _k[e_id] = trav1[e_id];
        _om[e_id] =   trav2[e_id]
                    / (cs_turb_cmu * cs_math_fmax(1.e-10, trav1[e_id]));
      }

      ipos += 2;
    }
  }

  /* Rij-epsilon
     ----------- */

  else if (cs_glob_turb_model->itytur == 3) {

    cs_real_63_t *grad_rij = nullptr;
    BFT_MALLOC(grad_rij, cs_glob_mesh->n_cells_with_ghosts, cs_real_63_t);

    cs_field_gradient_tensor(CS_F_(rij), iprev, inc, grad_rij);

    cs_field_gradient_scalar(CS_F_(eps), iprev, inc, grad);

    const cs_real_6_t *v_rij = (const cs_real_6_t *)CS_F_(rij)->val;
    const cs_real_t *v_eps = CS_F_(eps)->val;

    for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
      cs_lnum_t c_id = b_faces_dist_ids[e_id];

      const cs_real_3_t *_grad_rij = grad_rij[c_id];
      const cs_real_t *_rij = v_rij[c_id];
      const cs_real_t _eps = v_eps[c_id];

      const cs_real_t *xyzjjp = distant_dist_fbr + 3 * e_id;

      trav1[e_id] = _rij[0] + cs_math_3_dot_product(xyzjjp, _grad_rij[0]);
      trav2[e_id] = _rij[1] + cs_math_3_dot_product(xyzjjp, _grad_rij[1]);
      trav3[e_id] = _rij[2] + cs_math_3_dot_product(xyzjjp, _grad_rij[2]);
      trav4[e_id] = _rij[3] + cs_math_3_dot_product(xyzjjp, _grad_rij[3]);
      trav5[e_id] = _rij[4] + cs_math_3_dot_product(xyzjjp, _grad_rij[4]);
      trav6[e_id] = _rij[5] + cs_math_3_dot_product(xyzjjp, _grad_rij[5]);

      trav7[e_id] = _eps + cs_math_3_dot_product(xyzjjp, grad[3 * c_id]);
    }

    /* Free array */
    BFT_FREE(grad_rij);

    if (itytur_cpl == 3) {
      /* Coupled Reynolds tensor */
      cs_real_t *_r11 = rvdis[ipos];
      cs_real_t *_r22 = rvdis[ipos + 1];
      cs_real_t *_r33 = rvdis[ipos + 2];
      cs_real_t *_r12 = rvdis[ipos + 3];
      cs_real_t *_r13 = rvdis[ipos + 4];
      cs_real_t *_r23 = rvdis[ipos + 5];
      cs_real_t *_eps = rvdis[ipos + 6];

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        _r11[e_id] = trav1[e_id];
        _r22[e_id] = trav2[e_id];
        _r33[e_id] = trav3[e_id];
        _r12[e_id] = trav4[e_id];
        _r13[e_id] = trav5[e_id];
        _r23[e_id] = trav6[e_id];
        _eps[e_id] = trav7[e_id];
      }

      ipos += 7;
    }
    else if (itytur_cpl == 2) {
      cs_real_t *_k = rvdis[ipos];
      cs_real_t *_eps = rvdis[ipos + 1];

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        _k[e_id] = 0.5*(trav1[e_id] + trav2[e_id] + trav3[e_id]);
        _eps[e_id] = trav7[e_id];
      }
      ipos += 2;
    }
    else if (coupl->iturb == CS_TURB_K_OMEGA) {
      cs_real_t *_k = rvdis[ipos];
      cs_real_t *_om = rvdis[ipos + 1];

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        cs_real_t k = 0.5*(trav1[e_id] + trav2[e_id] + trav3[e_id]);
        _k[e_id] = k;
        _om[e_id] = trav7[e_id] / (cs_turb_cmu * cs_math_fmax(1.e-10, k));
      }

      ipos += 2;
    }
  }

  /* V2F-Phi
     ------- */

  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {

    /* Interpolate k */
    cs_field_gradient_scalar(CS_F_(k), iprev, inc, grad);
    _sat_coupling_interpolate_scalar_simple(CS_F_(k)->val, grad,
                                            n_b_faces_dist, b_faces_dist_ids,
                                            distant_dist_fbr, trav1);

    /* Interpolate epsilon */
    cs_field_gradient_scalar(CS_F_(eps), iprev, inc, grad);
    _sat_coupling_interpolate_scalar_simple(CS_F_(eps)->val, grad,
                                            n_b_faces_dist, b_faces_dist_ids,
                                            distant_dist_fbr, trav2);

    /* Interpolate Phi */
    cs_field_gradient_scalar(CS_F_(phi), iprev, inc, grad);
    _sat_coupling_interpolate_scalar_simple(CS_F_(phi)->val, grad,
                                            n_b_faces_dist, b_faces_dist_ids,
                                            distant_dist_fbr, trav3);
    /* Interpolate F-bar */
    cs_field_gradient_scalar(CS_F_(f_bar), iprev, inc, grad);
    _sat_coupling_interpolate_scalar_simple(CS_F_(f_bar)->val, grad,
                                            n_b_faces_dist, b_faces_dist_ids,
                                            distant_dist_fbr, trav4);


    /* Translation to coupled model */

    if (coupl->iturb == CS_TURB_V2F_PHI) {
      cs_real_t *_k    = rvdis[ipos];
      cs_real_t *_eps  = rvdis[ipos + 1];
      cs_real_t *_phi  = rvdis[ipos + 2];
      cs_real_t *_fbar = rvdis[ipos + 3];

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        _k[e_id]    = trav1[e_id];
        _eps[e_id]  = trav2[e_id];
        _phi[e_id]  = trav3[e_id];
        _fbar[e_id] = trav4[e_id];
      }

      ipos += 4;
    }
  }

  /* K-Omega
     ------- */

  else if (cs_glob_turb_model->iturb == CS_TURB_K_OMEGA) {

    /* Interpolate k */
    cs_field_gradient_scalar(CS_F_(k), iprev, inc, grad);
    _sat_coupling_interpolate_scalar_simple(CS_F_(k)->val, grad,
                                            n_b_faces_dist, b_faces_dist_ids,
                                            distant_dist_fbr, trav1);

    /* Interpolate omega */
    cs_field_gradient_scalar(CS_F_(omg), iprev, inc, grad);
    _sat_coupling_interpolate_scalar_simple(CS_F_(omg)->val, grad,
                                            n_b_faces_dist, b_faces_dist_ids,
                                            distant_dist_fbr, trav2);

    /* Translation to coupled model */

    if (coupl->iturb == CS_TURB_K_OMEGA) {
      cs_real_t *_k = rvdis[ipos];
      cs_real_t *_omg = rvdis[ipos + 1];

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        _k[e_id] = trav1[e_id];
        _omg[e_id] = trav2[e_id];
      }

      ipos += 2;
    }
    else if (itytur_cpl == 2) {
      /* k-omega => k-epsilon */
      cs_real_t *_k = rvdis[ipos];
      cs_real_t *_eps = rvdis[ipos + 1];

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        _k[e_id] = trav1[e_id];
        _eps[e_id] = trav2[e_id] * cs_turb_cmu * trav1[e_id];
      }

      ipos += 2;
    }
    else if (itytur_cpl == 3) {
      /* k-omega => Rij-epsilon */
      cs_real_t *_r11 = rvdis[ipos];
      cs_real_t *_r22 = rvdis[ipos + 1];
      cs_real_t *_r33 = rvdis[ipos + 2];

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        cs_real_t r0 = trav1[e_id] * 2. / 3.;
        _r11[e_id] = r0;
        _r22[e_id] = r0;
        _r33[e_id] = r0;
      }

      ipos += 3;

      /* R12, R13 and R23 using previously computed velocity gradient  */

      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        cs_lnum_t c_id = b_faces_dist_ids[e_id];

        cs_real_3_t *_grad_uvw = gradv[c_id];

        trav3[e_id] = _grad_uvw[0][1];
        trav4[e_id] = _grad_uvw[0][2];

        trav5[e_id] = _grad_uvw[1][0];
        trav6[e_id] = _grad_uvw[1][2];

        trav7[e_id] = _grad_uvw[2][0];
        trav8[e_id] = _grad_uvw[2][1];
      }

      cs_real_t *_r12 = rvdis[ipos];
      cs_real_t *_r13 = rvdis[ipos + 1];
      cs_real_t *_r23 = rvdis[ipos + 2];
      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++) {
        cs_real_t r0 = -1.0 * trav1[e_id] / cs_math_fmax(1.e-10, trav2[e_id]);

        _r12[e_id] = r0 * (trav3[e_id] + trav5[e_id]);
        _r13[e_id] = r0 * (trav4[e_id] + trav7[e_id]);
        _r23[e_id] = r0 * (trav6[e_id] + trav8[e_id]);
      }
      ipos += 3;

      cs_real_t *_eps = rvdis[ipos];
      for (cs_lnum_t e_id = 0; e_id < n_b_faces_dist; e_id++)
        _eps[e_id] = trav2[e_id] * cs_turb_cmu * trav1[e_id];

      ipos += 1;
    }
  }

  /* Scalars
     ------- */

  const int keysca = cs_field_key_id("scalar_id");
  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    int is_scalar = cs_field_get_key_int(f, keysca);
    if (is_scalar > -1) {
      cs_real_t *_scalar = rvdis[ipos];
      cs_field_gradient_scalar(f, iprev, inc, grad);
      _sat_coupling_interpolate_scalar_simple(f->val, grad, n_b_faces_dist,
                                              b_faces_dist_ids,
                                              distant_dist_fbr,
                                              _scalar);

      ipos += 1;
    }
  }

  /* Free memory */
  BFT_FREE(grad);
  BFT_FREE(gradv);
  BFT_FREE(trav1);
  BFT_FREE(trav2);
  BFT_FREE(trav3);
  BFT_FREE(trav4);
  BFT_FREE(trav5);
  BFT_FREE(trav6);
  BFT_FREE(trav7);
  BFT_FREE(trav8);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send varialbes data related to surface coupling.
 */
/*----------------------------------------------------------------------------*/

static void
_sat_interpolate_bc_from_b_face_data
(
  cs_field_t        *f,
  cs_real_t         *dt,
  int               *bc_type,
  const cs_lnum_t    n_b_faces_loc,
  const cs_lnum_t   *b_faces_loc_ids,
  const cs_lnum_t    n_b_faces_not_loc,
  const cs_lnum_t   *b_faces_not_loc_ids,
  const cs_real_t   *pond_fbr,
  const cs_real_t   *dofcpl,
  cs_real_t        **rvfbr
)
{
  /* Variable id key */
  const int kv = cs_field_key_id("variable_id");

  /* variable_id numbering starts at 1 because shared with fortran */
  int ivar = cs_field_get_key_int(f, kv) - 1;

  /* Mesh quantities */
  cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
  cs_real_t *diipb = cs_glob_mesh_quantities->diipb;
  cs_real_t *b_face_cog = cs_glob_mesh_quantities->b_face_cog;
  cs_real_t *cell_cen = cs_glob_mesh_quantities->cell_cen;

  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  const int cpl_bc_type
    = (cs_glob_sat_coupling_face_interpolation_type == 0) ?
       CS_COUPLED : CS_COUPLED_FD;

  /* Allocate work arrays */
  cs_real_3_t *grad = nullptr;
  BFT_MALLOC(grad, cs_glob_mesh->n_cells_with_ghosts, cs_real_3_t);

  cs_real_33_t *gradv = nullptr;
  BFT_MALLOC(gradv, cs_glob_mesh->n_cells_with_ghosts, cs_real_33_t);

  cs_real_63_t *gradt = nullptr;
  BFT_MALLOC(gradt, cs_glob_mesh->n_cells_with_ghosts, cs_real_63_t);

  /* Default parameters for gradient computation */

  switch (f->dim) {
  case 1:
    cs_field_gradient_scalar(f, false, 1, grad);
    break;
  case 3:
    cs_field_gradient_vector(f, false, 1, gradv);
    break;
  case 6:
    cs_field_gradient_tensor(f, false, 1, gradt);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Field '%s' has dimension %d, which is not compatible\n"
                "with code_saturne/code_saturne coupling.\n"),
              f->name, f->dim);
  }

  if (cs_glob_sat_coupling_face_interpolation_type == 1) {

    int fld_icodcl = (f->id == CS_F_(p)->id) ? -1 : 1;

    for (cs_lnum_t e_id = 0; e_id < n_b_faces_loc; e_id++) {
      cs_lnum_t f_id = b_faces_loc_ids[e_id];
      cs_lnum_t c_id = b_face_cells[f_id];

      cs_real_t *xyziip = diipb + 3 * f_id;

      /* Compute bnd face value using gradient projection */

      cs_real_t xip[6] = {0.};
      if (f->dim == 1) {
        /*
         * We want to prescribe a Dirichlet for pressure so as to conserve
         * the pressure gradient through the coupling and remain consistent
         * with the resolution of the pressure gradient on an orthogonal mesh.
         */
        xip[0] = f->val[c_id]
               + cs_math_3_dot_product(xyziip, grad[c_id]);
      }
      else if (f->dim == 3) {
        /*
         * For all other variables, we want to prescribe a Dirichlet matching
         * the convective fluxes at the center. We resrve a choice between
         * UPWIND, SOLU, and CENTERED. Only the centered case respects
         * the diffusion at the domain's interior faces.
         * For UPWIND and SOLU, the decentering
         * is done here and in bilsc2.f90 for coupled faces.
         *
         * UPWIND
         *
         *     xip =  cvar_var(iel)
         *
         * SOLU
         *     xip =   cvar_var(iel)
         *         + <gradv|xyzif>
         */
        /* Centered scheme */
        for (int i = 0; i < 3; i++)
          xip[i] = f->val[c_id * 3 + i]
                 + cs_math_3_dot_product(xyziip, gradv[c_id][i]);
      }
      else if (f->dim == 6) {
        for (int i = 0; i < 6; i++)
          xip[i] = f->val[c_id * 6 + i]
                 + cs_math_3_dot_product(xyziip, gradt[c_id][i]);
      }

      /* We need alpha_ij for centered interpolation
         and flumab for decentering */
      cs_real_t pondj = pond_fbr[e_id];

      cs_real_t xjp[6];
      for (int i = 0; i < f->dim; i++)
        xjp[i] = rvfbr[ivar + i][e_id];

      bc_type[f_id] = cpl_bc_type;

      f->bc_coeffs->icodcl[f_id] = fld_icodcl;
      for (int i = 0; i < f->dim; i++)
        f->bc_coeffs->rcodcl1[n_b_faces * i + f_id]
          = (1. - pondj) * xjp[i] + pondj * xip[i];
    }
  }
  else {

    int fld_icodcl = (f->id == CS_F_(p)->id) ? 3 : 1;

    for (cs_lnum_t e_id = 0; e_id < n_b_faces_loc; e_id++) {
      cs_lnum_t f_id = b_faces_loc_ids[e_id];
      cs_lnum_t c_id = b_face_cells[f_id];

      /* Test on atmospheric only applicable to iautom ? since bc_type is
       * replaced again.
       */
      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > CS_ATMO_OFF)
        cs_glob_bc_pm_info->iautom[f_id] = 1;

      cs_real_t *xyziip = diipb + 3 * f_id;

      cs_real_3_t xyzif = {b_face_cog[3 * f_id]     - cell_cen[3 * c_id],
                           b_face_cog[3 * f_id + 1] - cell_cen[3 * c_id + 1],
                           b_face_cog[3 * f_id + 2] - cell_cen[3 * c_id + 2]};

      cs_real_3_t xyzipf = {xyzif[0] - xyziip[0],
                            xyzif[1] - xyziip[1],
                            xyzif[2] - xyziip[2]};

      cs_real_t ipf = cs_math_3_norm(xyzipf);

      const cs_real_t *xyzopf = dofcpl + 3 * e_id;

      cs_real_t xip[6] = {0.};
      cs_real_t gradi[6] = {0.};

      switch (f->dim) {
      case 1:
        {
          /*FIXME Check why different formula for scalars ? */
          cs_real_t xyziipof[3] = {xyziip[0] + xyzopf[0],
                                   xyziip[1] + xyzopf[1],
                                   xyziip[2] + xyzopf[2]};
          xip[0] = f->val[c_id]
                 + cs_math_3_dot_product(xyziipof, grad[c_id]);

          gradi[0] = cs_math_3_dot_product(xyzipf, grad[c_id]) / ipf;
        }
        break;
      case 3:
        {
          cs_real_t *_vals = f->val + f->dim * c_id;
          for (int i = 0; i < f->dim; i++) {
            xip[i] = _vals[i]
                   + cs_math_3_dot_product(xyziip, gradv[c_id][i]);

            gradi[i] = cs_math_3_dot_product(xyzipf, gradv[c_id][i]) / ipf;
          }
        }
        break;
      case 6:
        {
          cs_real_t *_vals = f->val + f->dim * c_id;
          for (int i = 0; i < f->dim; i++) {
            xip[i] = _vals[i]
                   + cs_math_3_dot_product(xyziip, gradt[c_id][i]);

            gradi[i] = cs_math_3_dot_product(xyzipf, gradt[c_id][i]) / ipf;
          }
        }
        break;
      }

      cs_real_t xjp[6] = {0.};
      for (int i = 0; i < f->dim; i++)
        xjp[i] = rvfbr[ivar + i][e_id];

      bc_type[f_id] = cpl_bc_type;

      f->bc_coeffs->icodcl[f_id] = fld_icodcl;

      if (f->id != CS_F_(p)->id)
        for (int i = 0; i < f->dim; i++)
          f->bc_coeffs->rcodcl1[n_b_faces * i + f_id]
            = 0.5 * (xip[i] + xjp[i]);
      else {
        cs_real_t _coeff = -0.5 * dt[c_id];
        for (int i = 0; i < f->dim; i++)
          f->bc_coeffs->rcodcl3[n_b_faces * i + f_id]
            = _coeff * (gradi[i] + xjp[i]);
      }
    }
  }

  BFT_FREE(grad);
  BFT_FREE(gradv);
  BFT_FREE(gradt);

  /* Homogeneous neumann BC for non located faces */
  for (cs_lnum_t e_id = 0; e_id < n_b_faces_not_loc; e_id++) {
    cs_lnum_t f_id = b_faces_not_loc_ids[e_id];

    bc_type[f_id] = cpl_bc_type;

    f->bc_coeffs->icodcl[f_id] = 3;
    for (int i = 0; i < f->dim; i++)
      f->bc_coeffs->rcodcl3[n_b_faces * i + f_id] = 0.;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define new code_saturne coupling.
 *
 * The arguments to \ref cs_sat_coupling_define are:
 * \param[in] saturne_name          matching code_saturne application name
 * \param[in] boundary_cpl_criteria boundary face selection criteria for coupled
 *                                  faces, or nullptr
 * \param[in] volume_cpl_criteria   cell selection criteria for coupled cells, or
                                    nullptr
 * \param[in] boundary_loc_criteria boundary face selection criteria for location
 *                                  (not functional)
 * \param[in] volume_loc_criteria   cell selection criteria for location
 * \param[in] reverse               reverse mode if 1
 * \param[in] verbosity             verbosity level
 *
 * In the case of only 2 code_saturne instances, the 'saturne_name' argument
 * is ignored, as there is only one matching possibility.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * code_saturne instances based on the 'saturne_name' argument.
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_define(const char  *saturne_name,
                       const char  *boundary_cpl_criteria,
                       const char  *volume_cpl_criteria,
                       const char  *boundary_loc_criteria,
                       const char  *volume_loc_criteria,
                       int          reverse,
                       int          verbosity)
{
  _cs_sat_coupling_builder_t *scb = nullptr;

  /* Add corresponding coupling to temporary code_saturne couplings array */

  BFT_REALLOC(_sat_coupling_builder,
              _sat_coupling_builder_size + 1,
              _cs_sat_coupling_builder_t);

  scb = &(_sat_coupling_builder[_sat_coupling_builder_size]);

  scb->match_id = -1;

  scb->app_name = nullptr;
  if (saturne_name != nullptr) {
    BFT_MALLOC(scb->app_name, strlen(saturne_name) + 1, char);
    strcpy(scb->app_name, saturne_name);
  }

  scb->face_cpl_sel_c = nullptr;
  if (boundary_cpl_criteria != nullptr) {
    BFT_MALLOC(scb->face_cpl_sel_c, strlen(boundary_cpl_criteria) + 1, char);
    strcpy(scb->face_cpl_sel_c, boundary_cpl_criteria);
  }

  scb->cell_cpl_sel_c = nullptr;
  if (volume_cpl_criteria != nullptr) {
    BFT_MALLOC(scb->cell_cpl_sel_c, strlen(volume_cpl_criteria) + 1, char);
    strcpy(scb->cell_cpl_sel_c, volume_cpl_criteria);
  }

  scb->face_loc_sel_c = nullptr;
  if (boundary_loc_criteria != nullptr) {
    BFT_MALLOC(scb->face_loc_sel_c, strlen(boundary_loc_criteria) + 1, char);
    strcpy(scb->face_loc_sel_c, boundary_loc_criteria);
  }

  scb->cell_loc_sel_c = nullptr;
  if (volume_loc_criteria != nullptr) {
    BFT_MALLOC(scb->cell_loc_sel_c, strlen(volume_loc_criteria) + 1, char);
    strcpy(scb->cell_loc_sel_c, volume_loc_criteria);
  }

  scb->reverse = reverse;
  scb->verbosity = verbosity;

  _sat_coupling_builder_size += 1;
}

/*----------------------------------------------------------------------------
 * Get number of code_saturne couplings.
 *
 * returns:
 *   number of code_saturne couplings
 *----------------------------------------------------------------------------*/

int
cs_sat_coupling_n_couplings(void)
{
  return cs_glob_sat_n_couplings;
}

/*----------------------------------------------------------------------------
 * Get pointer to code_saturne coupling.
 *
 * parameters:
 *   coupling_id <-- Id (0 to n-1) of code_saturne coupling
 *
 * returns:
 *   pointer to code_saturne coupling structure
 *----------------------------------------------------------------------------*/

cs_sat_coupling_t *
cs_sat_coupling_by_id(int coupling_id)
{
  cs_sat_coupling_t  *retval = nullptr;

  if (   coupling_id > -1
      && coupling_id < cs_glob_sat_n_couplings)
    retval = cs_glob_sat_couplings[coupling_id];

  return retval;
}

/*----------------------------------------------------------------------------
 * Initialize code_saturne couplings.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_all_init(void)
{
  /* First, try using other MPI communicators */

#if defined(HAVE_MPI)

  if (_sat_coupling_builder_size > 0)
    _init_all_mpi_sat();

#endif

  /* Print unmatched instances */

  if (_sat_coupling_builder_size > 0) {

    bft_printf("Unmatched code_saturne couplings:\n"
               "---------------------------------\n\n");

    _print_all_unmatched_sat();

    bft_error(__FILE__, __LINE__, 0,
              _("At least 1 code_saturne coupling was defined for which\n"
                "no communication with a code_saturne instance is possible."));
  }
}

/*----------------------------------------------------------------------------*/
/*! \brief Initialization of main variables for code_saturne/code_saturne
 * coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_initialize
(
  void
)
{
  int nvar = 0;
  for (int field_id = 0; field_id < cs_field_n_fields(); field_id++) {
    cs_field_t *f = cs_field_by_id(field_id);
    /* Only couple field variables */
    if (!(f->type & CS_FIELD_VARIABLE))
      continue;
    if (f->type & CS_FIELD_CDO)
      continue;
    /* No coupling for mesh velocity */
    if (strcmp(f->name, "mesh_velocity") == 0)
      continue;
    nvar += f->dim;
  }

  for (int cpl_id = 0; cpl_id < cs_glob_sat_n_couplings; cpl_id++) {
    cs_sat_coupling_t *cpl = cs_glob_sat_couplings[cpl_id];

    cpl->icorio = -1;
    cpl->ale    = -1;
    cpl->imajcp = -1;
    cpl->nvarcp = -1;
    cpl->nvarto = -1;

    /* Face to face coupling should be defined for all couplings in the
     * same manner
     */
    int _face_interpolation_type_max = 0;
    _sat_coupling_int_max(cpl,
                          &cs_glob_sat_coupling_face_interpolation_type,
                          &_face_interpolation_type_max);
    cs_glob_sat_coupling_face_interpolation_type = _face_interpolation_type_max;

    /* Check if one of the instances is solved in a relative reference frame */
    _sat_coupling_int_max(cpl,
                          &(cs_glob_physical_constants->icorio),
                          &(cpl->icorio));

    /* Check coherency for reference frames of resolution. */
    if (cs_glob_physical_constants->icorio != cpl->icorio)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: Coupling is not available for different reference frames.\n"),
         __func__);

    /* Check for ALE/deformable mesh */
    _sat_coupling_int_max(cpl,
                          (const cs_lnum_t *)&(cs_glob_ale),
                          &(cpl->ale));

    if (   cpl->ale == CS_ALE_LEGACY
        || cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT) {
      cpl->imajcp = 1;
    }
    else
      cpl->imajcp = 0;

    /* Determine the number of coupled variables between instances of
     * coupling. All variables of a given instance are coupled except
     * for ALE, where the mesh velocity is not coupled.
     * Something should be done for specific physics.
     */
    cpl->nvarcp = nvar;
    if (cs_glob_ale != CS_ALE_NONE)
      cpl->nvarcp -= 3;

    /* Total number of sent variables */
    _sat_coupling_int_max(cpl, &(cpl->nvarcp), &(cpl->nvarto));

    /* Turbulence models */
    _sat_coupling_check_turbulence_coherency(cpl, cs_glob_turb_model->iturb);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the list of cells and boundary faces associated to a coupling
 * and a cloud of point.
 *
 * The local "support" cells and boundary faces are used to locate
 * the values in the distant "coupled" cells and faces.
 * Depending on the role of sender and/or receiver of the current process
 * in the coupling, some of these sets can be empty or not.
 *
 * The cell values are always located and interpolated on the distant
 * "cells" support. The face values are located and interpolated on
 * the distant "face" support if present, or on the distant "cell" support
 * if not.
 *
 * If the input arrays LCESUP and LFBSUP are not ordered, they will be
 * orderd in output.
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_locate_all
(
 void
)
{
  for (int cpl_id = 0; cpl_id < cs_glob_sat_n_couplings; cpl_id++) {

    /* If not first pass and mesh is not transient, nothing to do
     * for this coupling.
     */
    if (   _cs_sat_coupling_initialized == 1
        && cs_glob_sat_couplings[cpl_id]->imajcp == 0)
      continue;

    cs_lnum_t  ind;
    cs_lnum_t  nbr_fbr_cpl = 0, nbr_cel_cpl = 0;

    int  indic_glob[2] = {0, 0};
    int  indic_loc[2] = {0, 0};

    char coupled_mesh_name[64];
    cs_lnum_t *c_elt_list = nullptr;
    cs_lnum_t *f_elt_list = nullptr;
    cs_sat_coupling_t  *coupl = nullptr;
    fvm_nodal_t  *support_fbr = nullptr;
    cs_mesh_quantities_t  *mesh_quantities = cs_glob_mesh_quantities;

    int locator_options[PLE_LOCATOR_N_OPTIONS];
    locator_options[PLE_LOCATOR_NUMBERING] = 0;

    int *point_tag = nullptr;

    /* Initializations and verifications */

    coupl = cs_glob_sat_couplings[cpl_id];

    /* Removing the connectivity and localization informations in case of
       coupling update */

    if (coupl->cells_sup != nullptr) fvm_nodal_destroy(coupl->cells_sup);
    if (coupl->faces_sup != nullptr) fvm_nodal_destroy(coupl->faces_sup);

    /* Create the local lists */

    if (coupl->cell_loc_sel != nullptr) {

      BFT_MALLOC(c_elt_list, cs_glob_mesh->n_cells, cs_lnum_t);

      /* get the number  and a list of coupled local cells */
      cs_selector_get_cell_list(coupl->cell_loc_sel,
                                &(coupl->nbr_cel_sup),
                                c_elt_list);

    }

    if (coupl->face_loc_sel != nullptr) {

      BFT_MALLOC(f_elt_list, cs_glob_mesh->n_b_faces, cs_lnum_t);

      /* get the number and a list of coupled local faces */
      cs_selector_get_b_face_list(coupl->face_loc_sel,
                                  &(coupl->nbr_fbr_sup),
                                  f_elt_list);

    }

    if (coupl->nbr_cel_sup > 0) indic_loc[0] = 1; /* have coupled cells */
    if (coupl->nbr_fbr_sup > 0) indic_loc[1] = 1; /* have coupled faces */

    for (ind = 0; ind < 2; ind++)
      indic_glob[ind] = indic_loc[ind];

  #if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1)
      MPI_Allreduce (indic_loc, indic_glob, 2, MPI_INT, MPI_MAX,
                     cs_glob_mpi_comm);
  #endif

    /* One rank has coupled cells */
    if (indic_glob[0] > 0) {

      sprintf(coupled_mesh_name, _("coupled_cells_%d"), cpl_id+1);

      coupl->cells_sup = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                                        coupled_mesh_name,
                                                        false,
                                                        coupl->nbr_cel_sup,
                                                        c_elt_list);

    }

    /* One rank has coupled faces */
    if (indic_glob[1] > 0) {

      sprintf(coupled_mesh_name, _("coupled_faces_%d"), cpl_id + 1);

      coupl->faces_sup = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                        coupled_mesh_name,
                                                        false,
                                                        0,
                                                        coupl->nbr_fbr_sup,
                                                        nullptr,
                                                        f_elt_list);

    }

    if (coupl->cell_loc_sel != nullptr) BFT_FREE(c_elt_list);
    if (coupl->face_loc_sel != nullptr) BFT_FREE(f_elt_list);

    /* Build and initialize associated locator */

  #if defined(PLE_HAVE_MPI)

    if (coupl->localis_cel == nullptr)
      coupl->localis_cel = ple_locator_create(coupl->comm,
                                              coupl->n_sat_ranks,
                                              coupl->sat_root_rank);

    if (coupl->localis_fbr == nullptr)
      coupl->localis_fbr = ple_locator_create(coupl->comm,
                                              coupl->n_sat_ranks,
                                              coupl->sat_root_rank);

  #else

    if (coupl->localis_cel == nullptr)
      coupl->localis_cel = ple_locator_create();

    if (coupl->localis_fbr == nullptr)
      coupl->localis_fbr = ple_locator_create();

  #endif

    /* Initialization of the distant point localization */

    if (coupl->cell_cpl_sel != nullptr) {
      BFT_MALLOC(c_elt_list, cs_glob_mesh->n_cells, cs_lnum_t);
      cs_selector_get_cell_list(coupl->cell_cpl_sel,
                                &nbr_cel_cpl,
                                c_elt_list);
    }

    if (coupl->tag_func != nullptr) {
      BFT_MALLOC(point_tag, nbr_cel_cpl, int);
      coupl->tag_func(coupl->tag_context,
                      coupl->cells_sup,
                      nbr_cel_cpl,
                      0,
                      c_elt_list,
                      point_tag);
    }

    ple_locator_set_mesh(coupl->localis_cel,
                         coupl->cells_sup,
                         locator_options,
                         0.,
                         coupl->tolerance,
                         3,
                         nbr_cel_cpl,
                         c_elt_list,
                         point_tag,
                         mesh_quantities->cell_cen,
                         nullptr,
                         cs_coupling_mesh_extents,
                         cs_coupling_point_in_mesh_p);

    ple_locator_shift_locations(coupl->localis_cel, -1);

    BFT_FREE(point_tag);

    if (coupl->cell_cpl_sel != nullptr)
      BFT_FREE(c_elt_list);

    if (coupl->face_cpl_sel != nullptr) {
      BFT_MALLOC(f_elt_list, cs_glob_mesh->n_b_faces, cs_lnum_t);

      cs_selector_get_b_face_list(coupl->face_cpl_sel,
                                  &nbr_fbr_cpl,
                                  f_elt_list);
    }

    if (indic_glob[1] > 0)
      support_fbr = coupl->faces_sup;
    else
      support_fbr = coupl->cells_sup;

    if (coupl->tag_func != nullptr) {
      BFT_MALLOC(point_tag, nbr_fbr_cpl, int);
      coupl->tag_func(coupl->tag_context,
                      support_fbr,
                      nbr_fbr_cpl,
                      0,
                      f_elt_list,
                      point_tag);
    }

    ple_locator_set_mesh(coupl->localis_fbr,
                         support_fbr,
                         locator_options,
                         0.,
                         coupl->tolerance,
                         3,
                         nbr_fbr_cpl,
                         f_elt_list,
                         point_tag,
                         mesh_quantities->b_face_cog,
                         nullptr,
                         cs_coupling_mesh_extents,
                         cs_coupling_point_in_mesh_p);

    ple_locator_shift_locations(coupl->localis_fbr, -1);

    BFT_FREE(point_tag);

    if (coupl->face_cpl_sel != nullptr)
      BFT_FREE(f_elt_list);

    /* Computed some quantities needed for a centered-like interpolation */

    if (coupl->localis_fbr != nullptr)
      _sat_coupling_interpolate(coupl);

  #if 0
    /* TODO: associate the FVM meshes to the post-processing,
       with a function giving a pointer to the associated FVM structures,
       and another enabling its compacting or removing */
    {
      fvm_writer_t *w = fvm_writer_init("coupled_mesh",
                                        nullptr,
                                        "EnSight Gold",
                                        "binary",
                                        FVM_WRITER_FIXED_MESH);

      fvm_writer_export_nodal(w, coupl->cells_sup);
      fvm_writer_finalize(w);
    }
  #endif

    /* Compacting the interpolation support (could be removed) */

    if (coupl->cells_sup != nullptr)
      fvm_nodal_reduce(coupl->cells_sup, 1);
    if (coupl->faces_sup != nullptr)
      fvm_nodal_reduce(coupl->faces_sup, 1);

  #if 0 && defined(DEBUG) && !defined(NDEBUG)
    ple_locator_dump(coupl->localis_cel);
    ple_locator_dump(coupl->localis_fbr);
  #endif
    }

  _cs_sat_coupling_initialized = 1;
}

/*----------------------------------------------------------------------------
 * Create a sat_coupling_t structure.
 *
 * parameters:
 *   ref_axis           <-- reference axis
 *   face_sel_criterion <-- criterion for selection of boundary faces
 *   cell_sel_criterion <-- criterion for selection of cells
 *   sat_name           <-- code_saturne application name
 *   reverse            <-- reverse mode if 1
 *   verbosity          <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_add(const char  *face_cpl_sel_c,
                    const char  *cell_cpl_sel_c,
                    const char  *face_loc_sel_c,
                    const char  *cell_loc_sel_c,
                    const char  *sat_name,
                    int          reverse,
                    int          verbosity)
{
  cs_sat_coupling_t *sat_coupling = nullptr;

  /* Allocate _cs_sat_coupling_t structure */

  BFT_REALLOC(cs_glob_sat_couplings,
              cs_glob_sat_n_couplings + 1, cs_sat_coupling_t*);
  BFT_MALLOC(sat_coupling, 1, cs_sat_coupling_t);

  sat_coupling->sat_name = nullptr;
  sat_coupling->tag_func = nullptr;
  sat_coupling->tag_context = nullptr;
  sat_coupling->reverse = reverse;

  if (sat_name != nullptr) {
    BFT_MALLOC(sat_coupling->sat_name, strlen(sat_name) + 1, char);
    strcpy(sat_coupling->sat_name, sat_name);
  }

  /* Selection criteria  */

  sat_coupling->face_cpl_sel = nullptr;
  sat_coupling->cell_cpl_sel = nullptr;
  sat_coupling->face_loc_sel = nullptr;
  sat_coupling->cell_loc_sel = nullptr;

  if (face_cpl_sel_c != nullptr) {
    BFT_MALLOC(sat_coupling->face_cpl_sel, strlen(face_cpl_sel_c) + 1, char);
    strcpy(sat_coupling->face_cpl_sel, face_cpl_sel_c);
  }
  if (cell_cpl_sel_c != nullptr) {
    BFT_MALLOC(sat_coupling->cell_cpl_sel, strlen(cell_cpl_sel_c) + 1, char);
    strcpy(sat_coupling->cell_cpl_sel, cell_cpl_sel_c);
  }

  if (face_loc_sel_c != nullptr) {
    BFT_MALLOC(sat_coupling->face_loc_sel, strlen(face_loc_sel_c) + 1, char);
    strcpy(sat_coupling->face_loc_sel, face_loc_sel_c);
  }
  if (cell_loc_sel_c != nullptr) {
    BFT_MALLOC(sat_coupling->cell_loc_sel, strlen(cell_loc_sel_c) + 1, char);
    strcpy(sat_coupling->cell_loc_sel, cell_loc_sel_c);
  }

  sat_coupling->faces_sup = nullptr;
  sat_coupling->cells_sup = nullptr;

  sat_coupling->localis_fbr = nullptr;
  sat_coupling->localis_cel = nullptr;

  sat_coupling->nbr_fbr_sup = 0;
  sat_coupling->nbr_cel_sup = 0;

  sat_coupling->tolerance = 0.1;
  sat_coupling->verbosity = verbosity;

  /* Geometric quantities arrays for interpolation */

  sat_coupling->distant_dist_fbr = nullptr;
  sat_coupling->distant_of = nullptr;
  sat_coupling->local_of = nullptr;
  sat_coupling->distant_pond_fbr = nullptr;
  sat_coupling->local_pond_fbr = nullptr;

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
 * Create a new internal code_saturne coupling.
 *
 * arguments:
 *   tag_func          <-- pointer to tagging function
 *   tag_context       <-- pointer to tagging function context
 *   boundary_criteria <-- boundary face selection criteria, or nullptr
 *   volume_criteria   <-- volume cell selection criteria, or nullptr
 *   loc_tolerance     <-- location tolerance factor (0.1 recommended)
 *   reverse           <-- reverse mode if 1
 *   verbosity         <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_add_internal(cs_sat_coupling_tag_t  *tag_func,
                             void                   *tag_context,
                             const char             *boundary_cpl_criteria,
                             const char             *volume_cpl_criteria,
                             const char             *boundary_loc_criteria,
                             const char             *volume_loc_criteria,
                             float                   loc_tolerance,
                             int                     reverse,
                             int                     verbosity)
{
  CS_UNUSED(loc_tolerance);

  cs_sat_coupling_add(boundary_cpl_criteria,
                      volume_cpl_criteria,
                      boundary_loc_criteria,
                      volume_loc_criteria,
                      nullptr,
                      reverse,
                      verbosity);

  cs_sat_coupling_t *sat_coupling
    = cs_sat_coupling_by_id(cs_sat_coupling_n_couplings() - 1);

  sat_coupling->tag_func = tag_func;
  sat_coupling->tag_context = tag_context;

#if defined(HAVE_MPI)

  sat_coupling->comm = cs_glob_mpi_comm;
  sat_coupling->sat_root_rank = 0;
  sat_coupling->n_sat_ranks = cs_glob_n_ranks;

#endif
}

/*----------------------------------------------------------------------------
 * Array of reals exchange, associated to a given coupling.
 *
 * It is assumed that the arrays have the same size and the same values on
 * each group of processes (local and distant).
 *
 * int          cpl_id      : --> : coupling id (0-based)
 * int          nbrdis      : --> : number of values to send
 * int          nbrloc      : --> : number of values to receive
 * cs_real_t    vardis      : --> : distant values (to send)
 * cs_real_t    varloc      : <-- : local values (to receive)
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_array_exchange(int         cpl_id,
                               cs_lnum_t   nbrdis,
                               cs_lnum_t   nbrloc,
                               cs_real_t  *vardis,
                               cs_real_t  *varloc)
{
  cs_lnum_t  nbr = 0;
  bool  distant = false;

#if defined(HAVE_MPI)

  MPI_Status  status;

  cs_sat_coupling_t  *coupl = nullptr;

  /* Initializations and verifications */

  if (cpl_id < 0 || cpl_id >= cs_glob_sat_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling id %d; there are %d couplings"),
              cpl_id, cs_glob_sat_n_couplings);
  else
    coupl = cs_glob_sat_couplings[cpl_id];

  if (coupl->comm != MPI_COMM_NULL) {

    distant = true;

    /* Exchange between the groups master node */

    if (cs_glob_rank_id < 1)
      MPI_Sendrecv(vardis, nbrdis, CS_MPI_REAL, coupl->sat_root_rank, 0,
                   varloc, nbrloc, CS_MPI_REAL, coupl->sat_root_rank, 0,
                   coupl->comm, &status);

    /* Synchronization inside a group */

    if (cs_glob_n_ranks > 1)
      MPI_Bcast(varloc, nbrloc, CS_MPI_REAL, 0, cs_glob_mpi_comm);

  }

#endif /* defined(HAVE_MPI) */

  if (distant == false) {

    nbr = CS_MIN(nbrdis, nbrloc);

    for (cs_lnum_t ind = 0; ind < nbr; ind++)
      varloc[ind] = vardis[ind];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief code_saturne/code_saturne coupling using volumic source terms.
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_exchange_at_cells
(
  cs_field_t *f,    /*!<[in] pointer to cs_field_t */
  cs_real_t  *rhs,  /*!<[out] Explicit terms (RHS) */
  cs_real_t  *fimp  /*!<[out] Implicit source terms */
)
{
  /* This function should only be called for a variable field */
  if (!(f->type & CS_FIELD_VARIABLE))
    bft_error(__FILE__, __LINE__, 0,
              _("'%s' Cannot be called for non variable fields.\n"),
              __func__);

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  if (eqp->verbosity > 0)
    bft_printf("\n  Code-code coupling, add volume source terms for field %s\n\n", f->name);

  /* Loop on the different couplings */
  for (int cpl_id = 0; cpl_id < cs_glob_sat_n_couplings; cpl_id++) {
    cs_sat_coupling_t *cpl = cs_glob_sat_couplings[cpl_id];

    /* If the code_saturne/code_saturne coupling is a volumic one, handle
     * the exchange.
     */
    if (cpl->localis_cel != nullptr) {
      const int dim2 = f->dim * f->dim;

      ple_locator_t *cells_locator = cpl->localis_cel;

      const cs_lnum_t n_cells_loc
        = ple_locator_get_n_interior(cells_locator);
      const cs_lnum_t *cpl_cells_ids
        = ple_locator_get_interior_list(cells_locator);

      /* Prepare data for the send/recv operation */
      const cs_lnum_t n_cells_dist
        = ple_locator_get_n_dist_points(cells_locator);

      /* Count global number of elements, to ensure that a coupling is needed */
      cs_lnum_t _n_g_cells[2] = {n_cells_loc, n_cells_dist};
      cs_parall_sum(2, CS_LNUM_TYPE, _n_g_cells);

      cs_lnum_t n_g_cells_loc  = _n_g_cells[0];
      cs_lnum_t n_g_cells_dist = _n_g_cells[1];

      /* Get coordinates and indexes */
      const cs_lnum_t *dist_elt_ids
        = ple_locator_get_dist_locations(cells_locator);
      const cs_coord_t *dist_coords
        = ple_locator_get_dist_coords(cells_locator);

      /* Allocate temporary arrays */
      cs_real_t *cw1_dis = nullptr;
      BFT_MALLOC(cw1_dis, n_cells_dist * f->dim, cs_real_t);
      cs_array_real_fill_zero(n_cells_dist * f->dim, cw1_dis);

      cs_real_t *cw2_dis = nullptr;
      BFT_MALLOC(cw2_dis, n_cells_dist * dim2, cs_real_t);
      cs_array_real_fill_zero(n_cells_dist * dim2, cw2_dis);

      cs_real_t *cw1_loc = nullptr;
      BFT_MALLOC(cw1_loc, n_cells_loc * f->dim, cs_real_t);
      cs_array_real_fill_zero(n_cells_loc * f->dim, cw1_loc);

      cs_real_t *cw2_loc = nullptr;
      BFT_MALLOC(cw2_loc, n_cells_loc * dim2, cs_real_t);
      cs_array_real_fill_zero(n_cells_loc * dim2, cw2_loc);

      if (n_g_cells_dist > 0 && cpl->reverse == 0) {
        _sat_coupling_compute_data_at_cells(f,
                                            cpl->reverse,
                                            n_cells_dist,
                                            dist_elt_ids,
                                            dist_coords,
                                            cw1_dis,
                                            cw2_dis);
      }

      if (n_g_cells_loc > 0 && cpl->reverse == 1) {
        _sat_coupling_compute_data_at_cells(f,
                                            cpl->reverse,
                                            n_cells_loc,
                                            cpl_cells_ids,
                                            dist_coords,
                                            cw1_loc,
                                            cw2_loc);
      }

      /* Symmetric call with a test on both loc and dist sum */
      if (n_g_cells_loc > 0 || n_g_cells_dist > 0) {
        /* First call */
        cs_real_t *_cw1_dis = (n_cells_dist > 0) ? cw1_dis : nullptr;
        cs_real_t *_cw1_loc = (n_cells_loc > 0) ? cw1_loc : nullptr;
        ple_locator_exchange_point_var(cells_locator,
                                       _cw1_dis,
                                       _cw1_loc,
                                       nullptr,
                                       sizeof(cs_real_t),
                                       f->dim,
                                       cpl->reverse);

        /* Second call for the implicit part in reverse mode */
        if (cpl->reverse == 1) {
          cs_real_t *_cw2_dis = (n_cells_dist > 0) ? cw2_dis : nullptr;
          cs_real_t *_cw2_loc = (n_cells_loc > 0) ? cw2_loc : nullptr;
          ple_locator_exchange_point_var(cells_locator,
                                         _cw2_dis,
                                         _cw2_loc,
                                         nullptr,
                                         sizeof(cs_real_t),
                                         dim2,
                                         cpl->reverse);
        }
      }

      /* -------------------------------------------------- */
      /* Compute source terms based on the exchanged values */
      /* -------------------------------------------------- */

      cs_real_t *cvar_rho = CS_F_(rho)->val;
      const cs_real_t one_ov_xtau = 1./(100.0 * cs_glob_time_step->dt_ref);
      cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

      /*  Compute source term (n_cells_loc contribution) */
      if (cpl->reverse == 0 && n_g_cells_loc > 0) {
        if (f->dim == 1) {
          /* Scalars */
          for (cs_lnum_t e_id = 0; e_id < n_cells_loc; e_id++) {
            cs_lnum_t c_id = cpl_cells_ids[e_id];
            cs_real_t rovtau = cell_f_vol[c_id] * cvar_rho[c_id] *  one_ov_xtau;

            rhs[c_id] += rovtau * cw1_loc[e_id];
            fimp[c_id] -= rovtau;
          }
        }
        else {
          /* Vectors */
          for (cs_lnum_t e_id = 0; e_id < n_cells_loc; e_id++) {
            cs_lnum_t c_id = cpl_cells_ids[e_id];
            cs_real_t rovtau = cell_f_vol[c_id] * cvar_rho[c_id] * one_ov_xtau;
            for (int i = 0; i < f->dim; i++) {
              rhs[f->dim * c_id + i] += rovtau * cw1_loc[f->dim * e_id + i];
              fimp[dim2 * c_id + f->dim * i + i] -= rovtau;
            }
          }
        }
      }

      /* Compute source term in reverse mode (n_cells_dist contribution) */
      if (cpl->reverse == 1 && n_g_cells_dist > 0) {
        if (f->dim == 1) {
          /* Scalars */
          for (cs_lnum_t e_id = 0; e_id < n_cells_dist; e_id++) {
            cs_lnum_t c_id = dist_elt_ids[e_id];

            rhs[c_id]  += cw1_dis[e_id] * one_ov_xtau;
            fimp[c_id] -= cw2_dis[e_id] * one_ov_xtau;
          }
        }
        else {
          /* Vectors */
          for (cs_lnum_t e_id = 0; e_id < n_cells_dist; e_id++) {
            cs_lnum_t c_id = dist_elt_ids[e_id];
            for (int i = 0; i < f->dim; i++) {
              rhs[c_id * f->dim + i] += cw1_dis[e_id * f->dim + i] * one_ov_xtau;
              for (int j = 0; j < f->dim; j++) {
                fimp[c_id * dim2 + i * f->dim + j]
                  -= cw2_dis[e_id * dim2 + i * f->dim + j] * one_ov_xtau;
              }
            }
          }
        }
      }

      /* Free temporary work arrays */
      BFT_FREE(cw1_dis);
      BFT_FREE(cw2_dis);
      BFT_FREE(cw1_loc);
      BFT_FREE(cw2_loc);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief code_saturne/code_saturne boundary coupling initialization call
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_bnd_initialize
(
  int *bc_type /*!<[in] boundary face types */
)
{
  /* Global value */
  const int cpl_bc_type
    = (cs_glob_sat_coupling_face_interpolation_type == 0) ?
      CS_COUPLED : CS_COUPLED_FD;

  /* Loop on the different couplings */
  for (int cpl_id = 0; cpl_id < cs_glob_sat_n_couplings; cpl_id++) {
    cs_sat_coupling_t *cpl = cs_glob_sat_couplings[cpl_id];

    /* Check that this is a surface coupling */
    if (cpl->localis_fbr != nullptr) {
      ple_locator_t *locator = cpl->localis_fbr;

      /* Local coupled faces */
      cs_lnum_t n_b_faces_loc = ple_locator_get_n_interior(locator);
      const cs_lnum_t *b_faces_loc_ids = ple_locator_get_interior_list(locator);

      cs_lnum_t n_b_faces_not_loc = ple_locator_get_n_exterior(locator);
      const cs_lnum_t *b_faces_not_loc_ids
        = ple_locator_get_exterior_list(locator);

      if (n_b_faces_loc > 0) {
        cs_field_build_bc_codes_all();

        for (int field_id = 0; field_id < cs_field_n_fields(); field_id++) {
          cs_field_t *f = cs_field_by_id(field_id);

          /* Only couple field variables */
          if (!(f->type & CS_FIELD_VARIABLE))
            continue;

          /* No coupling for mesh velocity */
          if (strcmp(f->name, "mesh_velocity") == 0)
            continue;

          int fld_icodcl = 1;
          if (cs_glob_sat_coupling_face_interpolation_type != 1) {
            if (f->id == CS_F_(p)->id) {
              fld_icodcl = 3;
            }
          }

          for (cs_lnum_t e_id = 0; e_id < n_b_faces_loc; e_id++) {
            cs_lnum_t f_id = b_faces_loc_ids[e_id];
            bc_type[f_id] = cpl_bc_type;
            f->bc_coeffs->icodcl[f_id] = fld_icodcl;
          }

          /* Non-located boundary faces -> Homogeneous Neuman BC type */
          for (cs_lnum_t e_id = 0; e_id < n_b_faces_not_loc; e_id++) {
            cs_lnum_t f_id = b_faces_not_loc_ids[e_id];
            bc_type[f_id] = cpl_bc_type;
            f->bc_coeffs->icodcl[f_id] = 3;
          }
        }
      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief code_saturne/code_saturne coupling using boundary conditions
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_exchange_at_bnd_faces
(
  int       *bc_type, /*!<[in] boundary face types */
  cs_real_t *dt       /*!<[in] time step (per cell) */
)
{
  /* Loop on the different couplings */
  for (int cpl_id = 0; cpl_id < cs_glob_sat_n_couplings; cpl_id++) {
    cs_sat_coupling_t *cpl = cs_glob_sat_couplings[cpl_id];

    /* Check that this is a surface coupling */
    if (cpl->localis_fbr != nullptr) {
      ple_locator_t *locator = cpl->localis_fbr;

      /* Sanity check */
      if (cpl->nbr_fbr_sup > 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("%s: code_saturne/code_saturne coupling is not available "
                    "for boundary faces based values.\n"),
                  __func__);

      /* Local coupled faces */
      cs_lnum_t n_b_faces_loc = ple_locator_get_n_interior(locator);
      const cs_lnum_t *b_faces_loc_ids
        = ple_locator_get_interior_list(locator);

      /* Local non-coupled faces */
      cs_lnum_t n_b_faces_not_loc = ple_locator_get_n_exterior(locator);
      const cs_lnum_t *b_faces_not_loc_ids
        = ple_locator_get_exterior_list(locator);

      /* Distant coupled faces */
      cs_lnum_t n_b_faces_dist = ple_locator_get_n_dist_points(locator);
      const cs_lnum_t *b_faces_dist_ids
        = ple_locator_get_dist_locations(locator);

      /* Get coordinates */
      const cs_coord_t *coords = ple_locator_get_dist_coords(locator);

      /* Distance related pointers */
      cs_real_t *distant_dist_fbr = cpl->distant_dist_fbr; /* JJ' distance */
      cs_real_t *distant_of = cpl->distant_of;             /* OF distance */
      cs_real_t *distant_pond_fbr = cpl->distant_pond_fbr; /* Weighting
                                                              coefficient */

      /* Count global number of elements, to ensure that a coupling is needed */
      cs_gnum_t _n_g_b_faces[2] = {0, 0};
      _n_g_b_faces[0] = n_b_faces_loc;
      _n_g_b_faces[1] = n_b_faces_dist;
      cs_parall_sum(2, CS_GNUM_TYPE, _n_g_b_faces);

      cs_gnum_t n_g_b_faces_loc  = _n_g_b_faces[0];
      cs_gnum_t n_g_b_faces_dist = _n_g_b_faces[1];

      /* Allocate work arrays */
      cs_lnum_t _size_dis = (n_b_faces_dist > 0) ? n_b_faces_dist : 1;
      cs_real_t **rvdis = nullptr;
      BFT_MALLOC(rvdis, cpl->nvarto, cs_real_t *);

      cs_lnum_t _size_loc = (n_b_faces_loc > 0) ? n_b_faces_loc : 1;
      cs_real_t **rvfbr = nullptr;
      BFT_MALLOC(rvfbr, cpl->nvarto, cs_real_t *);

      for (int i = 0; i < cpl->nvarto; i++) {
        BFT_MALLOC(rvdis[i], _size_dis, cs_real_t);
        BFT_MALLOC(rvfbr[i], _size_loc, cs_real_t);
      }

      /* Exchange variables */
      if (n_g_b_faces_dist > 0)
        _sat_coupling_send_bnd_data(cpl, coords,
                                    n_b_faces_dist, b_faces_dist_ids,
                                    distant_dist_fbr, distant_pond_fbr,
                                    distant_of, rvdis);

      if (n_g_b_faces_dist > 0 || n_g_b_faces_loc > 0) {
        for (int i = 0; i < cpl->nvarto; i++) {
          cs_real_t *val_dist = (n_b_faces_dist > 0) ? rvdis[i] : nullptr;
          cs_real_t *val_loc = (n_b_faces_loc > 0) ? rvfbr[i] : nullptr;
          ple_locator_exchange_point_var(locator,
                                         val_dist,
                                         val_loc,
                                         nullptr,
                                         sizeof(cs_real_t),
                                         1, /* stride */
                                         cpl->reverse);
        }
      }

      for (int i = 0; i < cpl->nvarto; i++)
        BFT_FREE(rvdis[i]);
      BFT_FREE(rvdis);

      /* ----------------------------------------------- */
      /* Set boundary conditions based on exchanged data */
      /* ----------------------------------------------- */

      if (n_g_b_faces_loc > 0) {
        cs_real_t *local_pond_fbr = cpl->local_pond_fbr;
        cs_real_t *dofcpl = cpl->local_of;

        /* Ensure boundary conditions arrays are allocated */
        cs_field_build_bc_codes_all();

        for (int field_id = 0; field_id < cs_field_n_fields(); field_id++) {
          cs_field_t *f = cs_field_by_id(field_id);

          /* Only couple field variables */
          if (!(f->type & CS_FIELD_VARIABLE))
            continue;

          /* No coupling for mesh velocity */
          if (strcmp(f->name, "mesh_velocity") == 0)
            continue;

          _sat_interpolate_bc_from_b_face_data(f, dt, bc_type,
                                               n_b_faces_loc,
                                               b_faces_loc_ids,
                                               n_b_faces_not_loc,
                                               b_faces_not_loc_ids,
                                               local_pond_fbr,
                                               dofcpl,
                                               rvfbr);
        }
      }

      for (int i = 0; i < cpl->nvarto; i++)
        BFT_FREE(rvfbr[i]);
      BFT_FREE(rvfbr);
    }
  }
}

/*----------------------------------------------------------------------------
 * Destroy all couplings
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_all_finalize(void)
{
  for (int i = 0; i < cs_glob_sat_n_couplings; i++)
    _sat_coupling_destroy(cs_glob_sat_couplings[i]);

  BFT_FREE(cs_glob_sat_couplings);

  cs_glob_sat_n_couplings = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
