/*============================================================================
 * SYRTHES coupling
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>
#include <ple_coupling.h>
#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal.h"
#include "fvm_nodal_extract.h"
#include "fvm_nodal_project.h"

#if defined(HAVE_MPI)
#include "cs_coupling.h"
#endif

#include "cs_cf_thermo.h"
#include "cs_coupling.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_ht_convert.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_physical_model.h"
#include "cs_thermal_model.h"
#include "cs_timer_stats.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_syr_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

const int  cs_syr_coupling_tag = 'C'+'S'+'_'+'C'+'O'+'U'+'P'+'L'+'A'+'G'+'E';

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Structure associated with Syrthes coupling */

typedef struct {

  /* Mesh-related members */

  ple_locator_t     *locator;        /* Associated locator */

  int                elt_dim;        /* Element dimension */
  cs_lnum_t          n_elts;         /* Number of coupled elements */

  fvm_nodal_t       *elts;           /* Coupled elements */

  /* Saved arrays for post processing (float for reduced memory use) */

  int                post_mesh_id;   /* 0 if post-processing is not active,
                                        or post-processing mesh_id (< 0) */
  cs_real_t         *solid_temp;     /* Solid temperature received
                                        from SYRTHES */
  float             *flux;           /* Flux (calculated) */
  float             *tfluid_tmp;     /* Fluid temperature (points to flux in
                                        transient stage where solid_temp and
                                        fluid_temp are known, NULL once
                                        flux is calculated) */

  /* Saved array for volume coupling. Will be used to build source term. */

  double            *hvol;           /* Volumetric exchange coefficient. */

} cs_syr_coupling_ent_t ;

/* Structure associated with Syrthes coupling */

typedef struct {

  /* Mesh-related members */

  int                     dim;             /* Coupled mesh dimension */
  int                     ref_axis;        /* Axis for edge extraction */

  char                   *syr_name;        /* Application name */

  int                     n_b_locations;   /* Numbero of boundary locations */
  int                     n_v_locations;   /* Numbero of volume locations */
  int                    *b_location_ids;  /* Boundary location ids */
  int                    *v_location_ids;  /* Volume location ids */

  cs_syr_coupling_ent_t  *faces;           /* Wall coupling structure */
  cs_syr_coupling_ent_t  *cells;           /* Volume coupling structure */

  bool                    allow_nearest;   /* Allow nearest-neighbor
                                              mapping beyond basic matching
                                              tolerance */
  float                   tolerance;       /* Tolerance */
  int                     verbosity;       /* Verbosity level */
  int                     visualization;   /* Visualization output flag */

  /* Communication-related members */

#if defined(HAVE_MPI)

  MPI_Comm           comm;           /* Associated MPI communicator */

  int                n_syr_ranks;    /* Number of associated SYRTHES ranks */
  int                syr_root_rank;  /* First associated SYRTHES rank */

#endif

} cs_syr_coupling_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static int                  _syr_n_couplings = 0;
static cs_syr_coupling_t  **_syr_couplings = NULL;

/* Start and end (negative) numbers associated with
   dedicated post processing meshes */

static int  _syr_post_mesh_ext[2] = {0, 1};

static int  _syr_coupling_conservativity = 0; /* No forcing by default */
static int  _syr_coupling_implicit = 1;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize communicator for Syrthes coupling
 *
 * parameters:
 *   syr_coupling  <-> Syrthes coupling structure
 *   coupling_id   <-- id of this coupling (for log file message)
 *----------------------------------------------------------------------------*/

static void
_init_comm(cs_syr_coupling_t  *syr_coupling,
           int                 coupling_id)

{
#if defined(HAVE_MPI)

  int  mpi_flag = 0;
  int local_range[2] = {-1, -1};
  int distant_range[2] = {-1, -1};

  MPI_Initialized(&mpi_flag);

  if (mpi_flag == 0)
    return;

  bft_printf(_(" SYRTHES coupling %d: initializing MPI communication ... "),
             coupling_id);
  bft_printf_flush();

  ple_coupling_mpi_intracomm_create(MPI_COMM_WORLD,
                                    cs_glob_mpi_comm,
                                    syr_coupling->syr_root_rank,
                                    &(syr_coupling->comm),
                                    local_range,
                                    distant_range);

  bft_printf(_("[ok]\n"));
  bft_printf(_("  Local ranks = [%d..%d], distant ranks = [%d..%d].\n\n"),
             local_range[0], local_range[1] - 1,
             distant_range[0], distant_range[1] - 1);
  bft_printf_flush();

  syr_coupling->n_syr_ranks = distant_range[1] - distant_range[0];
  syr_coupling->syr_root_rank = distant_range[0];

#endif
}

/*----------------------------------------------------------------------------
 * Free communicator for Syrthes coupling
 *
 * parameters:
 *   syr_coupling  <-> Syrthes coupling structure
 *----------------------------------------------------------------------------*/

static void
_finalize_comm(cs_syr_coupling_t  *syr_coupling)
{
#if defined(HAVE_MPI)

  if (syr_coupling == NULL)
    return;

  if (syr_coupling->comm != MPI_COMM_NULL) {
    MPI_Comm_free(&(syr_coupling->comm));
    syr_coupling->comm = MPI_COMM_NULL;
  }

#endif
}

/*----------------------------------------------------------------------------
 * Exchange synchronization messages between code_saturne and SYRTHES.
 *
 * parameters:
 *   syr_coupling  <--  Syrthes coupling structure
 *   op_name_send  <--  operation name to send, or NULL. Only the 32
 *                      first characters are sent if the nae is longer.
 *   op_name_recv  <--  operation name to receive, or NULL (size: 33)
 *----------------------------------------------------------------------------*/

static void
_exchange_sync(cs_syr_coupling_t  *syr_coupling,
               const char         *op_name_send,
               char               *op_name_recv)
{
#if defined(HAVE_MPI)

  if (cs_glob_rank_id < 1) {

    MPI_Status status;

    if (op_name_send != NULL) {

      char _op_name_send[33];
      strncpy(_op_name_send, op_name_send, 32);
      _op_name_send[32] = '\0';

    /* Exchange command messages */
      if (op_name_recv != NULL) {
        MPI_Sendrecv(_op_name_send, 32, MPI_CHAR,
                     syr_coupling->syr_root_rank, cs_syr_coupling_tag,
                     op_name_recv, 32, MPI_CHAR,
                     syr_coupling->syr_root_rank, cs_syr_coupling_tag,
                     syr_coupling->comm, &status);
      }

      else
        MPI_Send(_op_name_send, 32, MPI_CHAR,
                 syr_coupling->syr_root_rank, cs_syr_coupling_tag,
                 syr_coupling->comm);

    }
    else if (op_name_recv != NULL) {
      MPI_Recv(op_name_recv, 32, MPI_CHAR,
               syr_coupling->syr_root_rank, cs_syr_coupling_tag,
               syr_coupling->comm, &status);
    }

  }

  if (op_name_recv != NULL && cs_glob_rank_id > -1) {
    MPI_Bcast(op_name_recv, 32, MPI_CHAR, 0, cs_glob_mpi_comm);
    op_name_recv[32] = '\0';
  }

#endif
}

/*----------------------------------------------------------------------------
 * Post process variables associated with Syrthes couplings
 *
 * parameters:
 *   coupling        <--  Void pointer to SYRTHES coupling structure
 *   ts              <--  time step status structure, or NULL
 *----------------------------------------------------------------------------*/

static void
_cs_syr_coupling_post_function(void                  *coupling,
                               const cs_time_step_t  *ts)
{
  int type_id;

  const cs_syr_coupling_t  *syr_coupling = coupling;
  cs_syr_coupling_ent_t *coupling_ent = NULL;

  for (type_id = 0; type_id < 2; type_id++) {

    if (type_id == 0)
      coupling_ent = syr_coupling->faces;
    else
      coupling_ent = syr_coupling->cells;

    if (coupling_ent != NULL) {

      if (coupling_ent->post_mesh_id != 0) {

        const cs_real_t *cell_temp = NULL;
        const cs_real_t *face_temp = NULL;

        if (type_id == 0)
          face_temp = coupling_ent->solid_temp;
        else
          cell_temp = coupling_ent->solid_temp;

        cs_post_write_var(coupling_ent->post_mesh_id,
                          CS_POST_WRITER_ALL_ASSOCIATED,
                          _("Solid T"),
                          1,
                          false,
                          false,
                          CS_POST_TYPE_cs_real_t,
                          cell_temp,
                          NULL,
                          face_temp,
                          ts);

        if (type_id == 1)
          cs_post_write_var(coupling_ent->post_mesh_id,
                            CS_POST_WRITER_ALL_ASSOCIATED,
                            _("Solid heat flux"),
                            1,
                            false,
                            false,
                            CS_POST_TYPE_float,
                            coupling_ent->flux,
                            NULL,
                            NULL,
                            ts);

      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Initialize post-processing of a SYRTHES coupling
 *
 * parameters:
 *   syr_coupling <-- partially initialized SYRTHES coupling structure
 *   coupling_ent <-- associated coupling mesh entity
 *----------------------------------------------------------------------------*/

static void
_post_init(cs_syr_coupling_t      *syr_coupling,
           cs_syr_coupling_ent_t  *coupling_ent)
{
  int dim_shift = 0;
  int coupling_id = -1;

  const int writer_id = -1;
  const int writer_ids[] = {writer_id};

  const cs_real_t t_ref = cs_glob_fluid_properties->t0;

  /* Determine coupling id */

  for (coupling_id = 0;
       (   coupling_id < _syr_n_couplings
        && _syr_couplings[coupling_id] != syr_coupling);
       coupling_id++);

  /* Exit silently if associated writer is not available */

  if (cs_post_writer_exists(writer_id) != true)
    return;

  int t_top_id
    = cs_timer_stats_switch(cs_timer_stats_id_by_name("postprocessing_stage"));

  /* Initialize post processing flag */

  coupling_ent->post_mesh_id = cs_post_get_free_mesh_id();

  /* Allocate arrays if not already present */

  if (coupling_ent->n_elts > 0) {
    if (coupling_ent->solid_temp == NULL) { /* surface coupling */
      BFT_MALLOC(coupling_ent->solid_temp, coupling_ent->n_elts, cs_real_t);
      for (cs_lnum_t i = 0; i < coupling_ent->n_elts; i++)
        coupling_ent->solid_temp[i] = t_ref;
    }
    if (coupling_ent->elt_dim == syr_coupling->dim) { /* volume coupling */
      if (coupling_ent->flux == NULL) {
        BFT_MALLOC(coupling_ent->flux, coupling_ent->n_elts, float);
        for (cs_lnum_t i = 0; i < coupling_ent->n_elts; i++)
          coupling_ent->flux[i] = 0;
      }
    }
  }
  coupling_ent->tfluid_tmp = NULL;

  /* Associate external mesh description with post processing subsystem */

  if (syr_coupling->dim == 2)
    dim_shift = 1;

  cs_post_define_existing_mesh(coupling_ent->post_mesh_id,
                               coupling_ent->elts,
                               dim_shift,
                               false,
                               false,
                               1,
                               writer_ids);

  /* Register post processing function */

  cs_post_add_time_dep_output(_cs_syr_coupling_post_function,
                              (void *)syr_coupling);

  /* Update start and end (negative) numbers associated with
     dedicated post processing meshes */

  if (_syr_post_mesh_ext[0] == 0)
    _syr_post_mesh_ext[0] = coupling_ent->post_mesh_id;

  _syr_post_mesh_ext[1] = coupling_ent->post_mesh_id;

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------
 * Check if coupling location is complete
 *
 * parameters:
 *   syr_coupling <-- partially initialized SYRTHES coupling structure
 *   coupling_ent <-- coupling entity
 *   n_ext        --> number of unlocated points
 *   ext_syr      --> 1 if SYRTHES has some unlocted elements, 0 otherwise
 *
 * returns:
 *   true if location is complete, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_is_location_complete(cs_syr_coupling_t      *syr_coupling,
                      cs_syr_coupling_ent_t  *coupling_ent,
                      cs_gnum_t              *n_ext,
                      bool                   *ext_syr)
{
  bool location_complete = true;

  /* Check that all points are effectively located */

  ple_lnum_t n_exterior = ple_locator_get_n_exterior(coupling_ent->locator);
  *n_ext = n_exterior;

  char  op_name_send[32 + 1];
  char  op_name_recv[32 + 1];

  cs_parall_counter(n_ext, 1);

  if (*n_ext > 0) {
    strcpy(op_name_send, "coupling:location:incomplete");
    location_complete = false;
  }
  else
    strcpy(op_name_send, "coupling:location:ok");

  _exchange_sync(syr_coupling, op_name_send, op_name_recv);
  if (!strcmp(op_name_recv, "coupling:location:incomplete")) {
    location_complete = false;
    *ext_syr = true;
  }
  else
    *ext_syr = false;

  return location_complete;
}

/*----------------------------------------------------------------------------
 * Define nodal mesh for Syrthes coupling from selection criteria.
 *
 * parameters:
 *   syr_coupling    <-- partially initialized SYRTHES coupling structure
 *   n_locations     <-- number of associated locations
 *   location_ids    <-- associated location ids
 *   elt_dim         <-- element dimension
 *
 * returns:
 *   pointer to created Syrthes coupling entity helper structure
 *----------------------------------------------------------------------------*/

static cs_syr_coupling_ent_t *
_create_coupled_ent(cs_syr_coupling_t  *syr_coupling,
                    int                 n_locations,
                    int                 location_ids[],
                    int                 elt_dim)
{
  char *coupled_mesh_name = NULL;
  bool      ext_syr = false;
  cs_lnum_t n_exterior = 0;
  cs_gnum_t n_ext = 0;
  cs_coord_t *elt_centers = NULL;
  fvm_nodal_t *location_elts = NULL;
  float *cs_to_syr_dist = NULL;
  float *syr_to_cs_dist = NULL;

  bool location_complete = false;
  cs_syr_coupling_ent_t *coupling_ent = NULL;

  int locator_options[PLE_LOCATOR_N_OPTIONS];
  locator_options[PLE_LOCATOR_NUMBERING] = 0;

  assert(syr_coupling != NULL);

  /* Initialization */

  BFT_MALLOC(coupling_ent, 1, cs_syr_coupling_ent_t);

  coupling_ent->locator = NULL;
  coupling_ent->elt_dim = elt_dim;

  coupling_ent->n_elts = 0;
  coupling_ent->elts = NULL;

  coupling_ent->post_mesh_id = 0;
  coupling_ent->solid_temp = NULL;
  coupling_ent->flux = NULL;
  coupling_ent->tfluid_tmp = NULL;

  coupling_ent->hvol = NULL;

  if (syr_coupling->verbosity > 0) {
    bft_printf(_("\nExtracting coupled mesh             ..."));
    bft_printf_flush();
  }

  /* Select elements */

  cs_lnum_t  n_elts = 0;
  cs_lnum_t *elt_list = NULL;

  for (int l_i = 0; l_i < n_locations; l_i++)
    n_elts += cs_mesh_location_get_n_elts(location_ids[l_i])[0];

  BFT_MALLOC(elt_list, n_elts, cs_lnum_t);

  n_elts = 0;
  for (int l_i = 0; l_i < n_locations; l_i++) {
    int loc_id = location_ids[l_i];
    const cs_lnum_t n = cs_mesh_location_get_n_elts(loc_id)[0];
    const cs_lnum_t *ids = cs_mesh_location_get_elt_list(loc_id);
    if (ids != NULL) {
      for (cs_lnum_t i = 0; i < n; i++)
        elt_list[n_elts++] = ids[i] + 1;
    }
    else {
      for (cs_lnum_t i = 0; i < n; i++)
        elt_list[n_elts++] = i + 1;
    }
  }

  /* Creation of a new nodal mesh from selected cells */

  if (elt_dim == syr_coupling->dim) {

    BFT_MALLOC(coupled_mesh_name,
                 strlen(_("SYRTHES %s cells"))
               + strlen(syr_coupling->syr_name) + 1, char);
    sprintf(coupled_mesh_name, _("SYRTHES %s cells"), syr_coupling->syr_name);

    coupling_ent->n_elts = n_elts;

    coupling_ent->elts
      = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                       coupled_mesh_name,
                                       false,
                                       coupling_ent->n_elts,
                                       elt_list);

    /* Allocate additional buffers */

    BFT_MALLOC(coupling_ent->hvol, coupling_ent->n_elts, double);
    BFT_MALLOC(coupling_ent->solid_temp, coupling_ent->n_elts, cs_real_t);

  }

  /* Creation of a new nodal mesh from selected border faces */

  else if (elt_dim == syr_coupling->dim - 1) {

    BFT_MALLOC(coupled_mesh_name,
               strlen("SYRTHES  faces") + strlen(syr_coupling->syr_name) + 1,
               char);
    sprintf(coupled_mesh_name, _("SYRTHES %s faces"), syr_coupling->syr_name);

    coupling_ent->n_elts = n_elts;

    coupling_ent->elts
      = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                       coupled_mesh_name,
                                       false,
                                       0,
                                       coupling_ent->n_elts,
                                       NULL,
                                       elt_list);

  }

  BFT_FREE(elt_list);

  BFT_FREE(coupled_mesh_name);

  if (syr_coupling->verbosity > 0) {
    bft_printf(" [ok]\n");
    bft_printf_flush();
  }

  if (fvm_nodal_get_n_g_vertices(coupling_ent->elts) == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" Selected mesh locations:\n"
                " lead to an empty mesh for SYRTHES coupling .\n"
                " \"%s\"\n"),
              syr_coupling->syr_name);

  /* In case of 2D coupling, project coupled elements to 2D */

  location_elts = coupling_ent->elts;

  if (syr_coupling->dim == 2) {

    double  a[6];
    cs_lnum_t  n_errors = 0;

    if (syr_coupling->verbosity > 0) {
      bft_printf(_("Projecting the extracted mesh to 2D ..."));
      bft_printf_flush();
    }

    fvm_nodal_project(coupling_ent->elts, syr_coupling->ref_axis, &n_errors);

    if (n_errors > 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error projecting the extracted mesh."));

    if (syr_coupling->verbosity > 0) {
      bft_printf(_(" [ok]\n"));
      bft_printf_flush();
    }

    location_elts = fvm_nodal_copy(coupling_ent->elts);

    if (syr_coupling->ref_axis == 0) {
      a[0] = 0.; a[1] = 1.; a[2] = 0.; a[3] = 0.; a[4] = 0.; a[5] = 1.;
    }
    else if (syr_coupling->ref_axis == 1) {
      a[0] = 1.; a[1] = 0.; a[2] = 0.; a[3] = 0.; a[4] = 0.; a[5] = 1.;
    }
    else if (syr_coupling->ref_axis == 2) {
      a[0] = 1.; a[1] = 0.; a[2] = 0.; a[3] = 0.; a[4] = 1.; a[5] = 0.;
    }

    fvm_nodal_project_coords(location_elts, a);
  }

  /* Element information */

  if (syr_coupling->verbosity > 0) {
    cs_gnum_t n_g_elts = coupling_ent->n_elts;
    cs_parall_counter(&n_g_elts, 1);
    bft_printf(_("\nExtracted mesh built of %llu elements.\n"),
               (unsigned long long)n_g_elts);
    bft_printf_flush();
  }

  /* Initialize post-processing */

  /* Precaution: deactivate visualization for time-dependent meshes,
     as this would require regenerating visualization at each time step */

  if (cs_post_get_writer_time_dep(-1) != FVM_WRITER_FIXED_MESH)
    syr_coupling->visualization = 0;

  if (syr_coupling->visualization != 0)
    _post_init(syr_coupling, coupling_ent);

  /* Build and initialize associated locator */

  if (syr_coupling->verbosity > 0) {
    bft_printf(_("\nLocator structure and mesh creation ..."));
    bft_printf_flush();
  }

  /* Retrieve coordinates using FVM functions rather than previous list and
     coordinates, in case the extracted nodal mesh contains elements in a
     different order (multiple element types) or coordinates are projected
     in 2D. */

  if (coupling_ent->n_elts > 0) {

    if (syr_coupling->visualization != 0)
      BFT_MALLOC(cs_to_syr_dist, coupling_ent->n_elts, float);

    BFT_MALLOC(elt_centers,
               coupling_ent->n_elts*syr_coupling->dim,
               cs_coord_t);
    fvm_nodal_get_element_centers(location_elts,
                                  CS_INTERLACE,
                                  coupling_ent->elt_dim,
                                  elt_centers);
  }

  /* Locate entities */

#if defined(PLE_HAVE_MPI)
  coupling_ent->locator = ple_locator_create(syr_coupling->comm,
                                             syr_coupling->n_syr_ranks,
                                             syr_coupling->syr_root_rank);
#else
  coupling_ent->locator = ple_locator_create();
#endif

  ple_locator_set_mesh(coupling_ent->locator,
                       location_elts,
                       locator_options,
                       0.,
                       syr_coupling->tolerance,
                       syr_coupling->dim,
                       coupling_ent->n_elts,
                       NULL,
                       NULL,
                       elt_centers,
                       cs_to_syr_dist,
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh);

  /* Check that all points are effectively located */

  location_complete = _is_location_complete(syr_coupling,
                                            coupling_ent,
                                            &n_ext,
                                            &ext_syr);

  if (syr_coupling->allow_nearest) {

    float tolerance = syr_coupling->tolerance;

    while (location_complete == false) {

      tolerance *= 4;

      if (syr_coupling->verbosity > 0) {
        bft_printf(_(" [failed]\n"));
        if (n_ext > 0)
          bft_printf(_(" %llu fluid mesh elements not located on solid mesh\n"),
                     (unsigned long long) n_ext);
        if (ext_syr)
          bft_printf(_(" Some solid mesh elements not located on fluid mesh\n"));
        bft_printf(_("\n   Extending search with tolerance factor %f..."),
                   tolerance);
        bft_printf_flush();
      }

      ple_locator_extend_search(coupling_ent->locator,
                                location_elts,
                                locator_options,
                                0.,
                                tolerance,
                                coupling_ent->n_elts,
                                NULL,
                                NULL,
                                elt_centers,
                                cs_to_syr_dist,
                                cs_coupling_mesh_extents,
                                cs_coupling_point_in_mesh);

      location_complete = _is_location_complete(syr_coupling,
                                                coupling_ent,
                                                &n_ext,
                                                &ext_syr);

    }

  }

  if (syr_coupling->verbosity > 0) {
    bft_printf(_(" [ok]\n"));
    bft_printf_flush();
  }

  /* Shift from 1-base to 0-based locations */

  ple_locator_shift_locations(coupling_ent->locator, -1);

  if (location_elts != coupling_ent->elts)
    fvm_nodal_destroy(location_elts);

  if (elt_centers != NULL)
    BFT_FREE(elt_centers);

  bool default_writer_is_active = false;

  if (syr_coupling->visualization != 0) {

    default_writer_is_active
      = cs_post_writer_is_active(CS_POST_WRITER_DEFAULT);

    cs_post_activate_writer(CS_POST_WRITER_DEFAULT, true);
    cs_post_write_meshes(cs_glob_time_step);

    const float *b_dist = NULL, *v_dist = NULL;

    if (coupling_ent->elt_dim == syr_coupling->dim - 1)
      b_dist = cs_to_syr_dist;
    else if (coupling_ent->elt_dim == syr_coupling->dim)
      v_dist = cs_to_syr_dist;

    cs_post_write_var(coupling_ent->post_mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      _("distance_to_solid"),
                      1,
                      false,
                      false, /* use_parent, */
                      CS_POST_TYPE_float,
                      v_dist,
                      NULL,
                      b_dist,
                      NULL);  /* time-independent variable */

    BFT_FREE(cs_to_syr_dist);

  }

  /* Post-process distances from SYRTHES points to code_saturne faces */

  if (elt_dim == syr_coupling->dim - 1) {

    cs_lnum_t n_dist_elts = ple_locator_get_n_dist_points(coupling_ent->locator);

    BFT_MALLOC(syr_to_cs_dist, n_dist_elts, float);

    ple_locator_exchange_point_var(coupling_ent->locator,
                                   syr_to_cs_dist,
                                   NULL,
                                   NULL,
                                   sizeof(float),
                                   1,
                                   1);

    if (   syr_coupling->visualization != 0
        && syr_coupling->allow_nearest == false) {

      cs_lnum_t i;
      int writer_ids[] = {CS_POST_WRITER_DEFAULT};
      int mesh_id = coupling_ent->post_mesh_id - 1;
      cs_lnum_t *p_vtx_num = NULL;
      fvm_io_num_t *vtx_io_num = NULL;
      fvm_nodal_t *syr_points = fvm_nodal_create("SYRTHES face centers",
                                                 syr_coupling->dim);

      BFT_MALLOC(p_vtx_num, n_dist_elts, cs_lnum_t);

      for (i = 0; i < (cs_lnum_t)n_dist_elts; i++)
        p_vtx_num[i] = i+1;

      fvm_nodal_define_vertex_list(syr_points, n_dist_elts, p_vtx_num);
      fvm_nodal_set_shared_vertices
        (syr_points,
         ple_locator_get_dist_coords(coupling_ent->locator));

      if (cs_glob_n_ranks > 1) {

        vtx_io_num = fvm_io_num_create_from_scan(n_dist_elts);

        fvm_nodal_init_io_num(syr_points,
                              fvm_io_num_get_global_num(vtx_io_num),
                              0);

      }

      cs_post_define_existing_mesh(mesh_id,
                                   syr_points,
                                   0,
                                   true,
                                   false,
                                   1,
                                   writer_ids);

      cs_post_write_meshes(cs_glob_time_step);

      cs_post_write_vertex_var(mesh_id,
                               CS_POST_WRITER_ALL_ASSOCIATED,
                               _("distance_to_fluid"),
                               1,
                               false,
                               false, /* use parent */
                               CS_POST_TYPE_float,
                               syr_to_cs_dist,
                               NULL); /* time-independent variable */

      cs_post_free_mesh(mesh_id);

      if (cs_glob_n_ranks > 1)
        fvm_io_num_destroy(vtx_io_num);

    } /* Do post-processing */

    BFT_FREE(syr_to_cs_dist);

  }

  if (n_ext) {

    int i;
    int writer_ids[] = {-1};
    int mesh_id = cs_post_get_free_mesh_id();
    cs_lnum_t *post_vtx_num = NULL;
    cs_coord_t *exterior_coords = NULL;
    cs_coord_t *el_list = NULL;
    fvm_io_num_t *vtx_io_num = NULL;
    fvm_nodal_t *ulck_points = fvm_nodal_create("unlocated elements (centers)",
                                                3);
    n_exterior = ple_locator_get_n_exterior(coupling_ent->locator);
    const ple_lnum_t *exterior_list
      = ple_locator_get_exterior_list(coupling_ent->locator);

    BFT_MALLOC(post_vtx_num, n_exterior, cs_lnum_t);
    BFT_MALLOC(exterior_coords, 3*n_exterior, cs_coord_t);
    BFT_MALLOC(el_list,
               coupling_ent->n_elts*3,
               cs_coord_t);

    fvm_nodal_get_element_centers(coupling_ent->elts,
                                  CS_INTERLACE,
                                  coupling_ent->elt_dim,
                                  el_list);

    for (i = 0; i < (cs_lnum_t)n_exterior; i++) {
      post_vtx_num[i] = i+1;
      if (exterior_list[i] >= coupling_ent->n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: invalid exterior elements selection."));
      exterior_coords[3*i   ] = el_list[3*exterior_list[i]   ];
      exterior_coords[3*i +1] = el_list[3*exterior_list[i] +1];
      exterior_coords[3*i +2] = el_list[3*exterior_list[i] +2];
    }

    fvm_nodal_define_vertex_list(ulck_points,
                                 (cs_lnum_t)n_exterior,
                                 post_vtx_num);
    fvm_nodal_set_shared_vertices
      (ulck_points,
       exterior_coords);

    if (cs_glob_n_ranks > 1) {
      vtx_io_num = fvm_io_num_create_from_scan(n_exterior);
      fvm_nodal_init_io_num(ulck_points,
                            fvm_io_num_get_global_num(vtx_io_num),
                            0);
    }

    cs_post_define_existing_mesh(mesh_id,
                                 ulck_points,
                                 0,
                                 true,
                                 false,
                                 1,
                                 writer_ids);

    cs_post_write_meshes(cs_glob_time_step);
    cs_post_free_mesh(mesh_id);

    if (cs_glob_n_ranks > 1)
      fvm_io_num_destroy(vtx_io_num);

    BFT_FREE(el_list);
    BFT_FREE(exterior_coords);

    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Coupling with SYRTHES impossible:\n"
                 "%llu element centers from mesh \"%s\"\n"
                 "not located on SYRTHES mesh."),
               (unsigned long long)n_ext,
               fvm_nodal_get_name(coupling_ent->elts));

  }

  /* Restore writer active status */

  if (syr_coupling->visualization != 0)
    cs_post_activate_writer(CS_POST_WRITER_DEFAULT, default_writer_is_active);

  /* Ensure clean stop */

  if (location_complete == false)
    cs_coupling_set_sync_flag(PLE_COUPLING_STOP);

  return coupling_ent;
}

/*----------------------------------------------------------------------------
 * Log timing info
 *----------------------------------------------------------------------------*/

static void
_all_comm_times(void)
{
  int coupl_id, ent_id;
  cs_syr_coupling_t *syr_coupling = NULL;

  if (_syr_n_couplings == 0)
    return;

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\nSYRTHES coupling overheads\n"));

  for (coupl_id = 0; coupl_id < _syr_n_couplings; coupl_id++) {

    syr_coupling = _syr_couplings[coupl_id];

    for (ent_id = 0; ent_id < 2; ent_id++) {

      cs_syr_coupling_ent_t
        *ce = (ent_id == 0) ? syr_coupling->faces : syr_coupling->cells;
      const char *ent_type[] = {N_("surface"), N_("volume")};

      if (ce != NULL) {

        double location_wtime, exchange_wtime;
        double location_comm_wtime, exchange_comm_wtime;

        if (syr_coupling->syr_name != NULL)
          cs_log_printf(CS_LOG_PERFORMANCE,
                        _("\n  %s (%s):\n\n"),
                        syr_coupling->syr_name, _(ent_type[ent_id]));
        else
          cs_log_printf(CS_LOG_PERFORMANCE,
                        _("\n  coupling %d (%s):\n\n"),
                        coupl_id, _(ent_type[ent_id]));

        ple_locator_get_times(ce->locator,
                              &location_wtime,
                              NULL,
                              &exchange_wtime,
                              NULL);

        ple_locator_get_comm_times(ce->locator,
                                   &location_comm_wtime,
                                   NULL,
                                   &exchange_comm_wtime,
                                   NULL);

        cs_log_printf(CS_LOG_PERFORMANCE,
                      _("    location time:                 %12.3f\n"
                        "      communication and wait:      %12.3f\n"
                        "    variable exchange time:        %12.3f\n"
                        "      communication and wait:      %12.3f\n"),
                      location_wtime, location_comm_wtime,
                      exchange_wtime, exchange_comm_wtime);

      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Destroy coupled entity helper structure.
 *
 * parameters:
 *   coupling ent <-> pointer to structure pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_coupled_ent(cs_syr_coupling_ent_t  **coupling_ent)
{
  cs_syr_coupling_ent_t *ce = *coupling_ent;

  if (ce == NULL)
    return;

  if (ce->locator != NULL)
    ce->locator = ple_locator_destroy(ce->locator);

  if (ce->solid_temp != NULL)
    BFT_FREE(ce->solid_temp);
  if (ce->flux != NULL)
    BFT_FREE(ce->flux);

  if (ce->hvol != NULL)
    BFT_FREE(ce->hvol);

  if (ce->elts != NULL)
    ce->elts = fvm_nodal_destroy(ce->elts);

  BFT_FREE(*coupling_ent);
}

/*----------------------------------------------------------------------------
 * Update post-processing variables of a Syrthes coupling
 *
 * Note that only the solid temperature is output for surface coupling,
 * while the heat flux is also output for volume coupling.
 *
 * parameters:
 *   coupling_ent <--  Syrthes coupling structure
 *   step         <--  0: var = wall temperature
 *                     1: var = fluid temperature
 *                     2: var = exchange coefficient
 *   var          <--  Pointer to variable values
 *----------------------------------------------------------------------------*/

static void
_post_var_update(cs_syr_coupling_ent_t  *coupling_ent,
                 int                     step,
                 const cs_real_t        *var)
{
  cs_lnum_t  n_elts, ii;

  assert(coupling_ent != NULL);

  if (coupling_ent->post_mesh_id == 0)
    return;

  assert(coupling_ent->solid_temp != NULL);
  assert(step == 0 || coupling_ent->flux != NULL);

  n_elts = coupling_ent->n_elts;

  /* Allocate arrays */

  switch(step) {

  case 0:
    for (ii = 0; ii < n_elts; ii++)
      coupling_ent->solid_temp[ii] = var[ii];
    break;

  case 1:
    coupling_ent->tfluid_tmp = coupling_ent->flux;
    for (ii = 0; ii < n_elts; ii++)
      coupling_ent->tfluid_tmp[ii] = var[ii];
    break;

  case 2:
    assert(coupling_ent->tfluid_tmp == coupling_ent->flux);
    for (ii = 0; ii < n_elts; ii++)
      coupling_ent->flux[ii] = var[ii] * (  coupling_ent->solid_temp[ii]
                                          - coupling_ent->flux[ii]);
    coupling_ent->tfluid_tmp = NULL;
    break;

  default:
    assert(0);
  }

}

/*----------------------------------------------------------------------------
 * Ensure conservativity thanks to a corrector coefficient computed by SYRTHES
 * SYRTHES computes a global flux for a given tfluid and hfluid field.
 * code_saturne sent before its computed global flux for this time step.
 *
 * parameters:
 *   syr_coupling   <-- SYRTHES coupling structure
 *   coupl_face_ids <-- ids of coupled boundary faces
 *----------------------------------------------------------------------------*/

static void
_ensure_conservativity(cs_syr_coupling_t   *syr_coupling,
                       const cs_lnum_t       coupl_face_ids[])
{
  cs_lnum_t ii, face_id;

  double g_flux = 0.0, _flux = 0.0, coef = 0.0;

  double  *surf = cs_glob_mesh_quantities->b_face_surf;
  cs_syr_coupling_ent_t  *coupling_ent = NULL;

  /* Sanity checks */

  assert(surf != NULL);
  assert(coupl_face_ids != NULL);
  assert(syr_coupling != NULL);
  coupling_ent = syr_coupling->faces;
  assert(coupling_ent != NULL);

  /* Compute code_saturne's global flux */

  for (ii = 0; ii < coupling_ent->n_elts; ii++) {
    face_id = coupl_face_ids[ii];
    _flux += coupling_ent->flux[ii] * surf[face_id];
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Reduce(&_flux, &g_flux, 1, MPI_DOUBLE, MPI_SUM, 0, cs_glob_mpi_comm);
#endif

  if (cs_glob_n_ranks == 1)
    g_flux = _flux;

  /* Send the computed global flux to SYRTHES and receive the corrector
     coefficient */

#if defined(HAVE_MPI)
  if (cs_glob_rank_id < 1) {

    MPI_Status  status;

    /* Send global flux */

    MPI_Send(&g_flux, 1, MPI_DOUBLE,
             syr_coupling->syr_root_rank,
             cs_syr_coupling_tag,
             syr_coupling->comm);

    if (syr_coupling->verbosity > 1)
      bft_printf(_(" Global heat flux exchanged with SYRTHES in W: %5.3e\n"),
                 g_flux);

    /* Receive corrector coefficient */

    MPI_Recv(&coef, 1, MPI_DOUBLE,
             syr_coupling->syr_root_rank,
             cs_syr_coupling_tag,
             syr_coupling->comm,
             &status);

  }
#endif

  /* Print message */

  if (syr_coupling->verbosity > 1)
    bft_printf(_(" Correction coefficient used to force conservativity during"
                 " coupling with SYRTHES: %5.3e\n"), coef);
}

/*----------------------------------------------------------------------------
 * Exchange location synchronization status
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *
 * returns:
 *   0 in case of success, 1 otherwise
 *----------------------------------------------------------------------------*/

static int
_sync_after_location(cs_syr_coupling_t  *syr_coupling)
{
  int retval = 1;

  char op_name_send[32 + 1];
  char op_name_recv[32 + 1];

  /* Communication with SYRTHES */
  /*----------------------------*/

  /* Ready to start time iterations */

  strcpy(op_name_send, "coupling:start");

  _exchange_sync(syr_coupling, op_name_send, op_name_recv);

  if (!strcmp(op_name_recv, "coupling:error:location")) {

    cs_coupling_set_sync_flag(PLE_COUPLING_STOP);

    cs_base_warn(__FILE__, __LINE__);

    bft_printf(_(" Message received from SYRTHES: \"%s\"\n"
                 " indicates meshes have not been matched correctly.\n\n"
                 " The calculation will not run.\n\n"),
               op_name_recv);

  }
  else if (strcmp(op_name_recv, "coupling:start"))
    bft_error(__FILE__, __LINE__, 0,
              _(" Message received from SYRTHES: \"%s\"\n"
                " indicates an error or is unexpected."),
              op_name_recv);

  else
    retval = 0;

  return retval;
}

/*----------------------------------------------------------------------------
 * Get pointer to SYRTHES coupling.
 *
 * parameters:
 *   coupling_id <-- Id (0 to n-1) of SYRTHES coupling
 *
 * returns:
 *   pointer to SYRTHES coupling structure
 *----------------------------------------------------------------------------*/

static cs_syr_coupling_t *
_syr_coupling_by_id(int  coupling_id)
{
  cs_syr_coupling_t  *retval = NULL;

  if (   coupling_id > -1
      && coupling_id < _syr_n_couplings)
    retval = _syr_couplings[coupling_id];

  return retval;
}

/*----------------------------------------------------------------------------
 * Create or redefine a syr_coupling_t structure.
 *
 * If a structure is redefined, associated locations are reset.
 *
 * parameters:
 *   dim                <-- spatial mesh dimension
 *   ref_axis           <-- reference axis
 *   syr_name           <-- SYRTHES application name
 *   allow_nonmatching  <-- nearest-neighbor search for non-matching faces flag
 *   tolerance          <-- addition to local extents of each element
 *                          extent = base_extent * (1 + tolerance)
 *   verbosity          <-- verbosity level
 *   visualization      <-- visualization output flag
 *----------------------------------------------------------------------------*/

static cs_syr_coupling_t *
_syr_coupling_define(int          dim,
                     int          ref_axis,
                     const char  *syr_name,
                     bool         allow_nonmatching,
                     float        tolerance,
                     int          verbosity,
                     int          visualization)
{
  cs_syr_coupling_t *syr_coupling = NULL;

  /* Search in existing couplings */

  for (int i = 0;
       i < _syr_n_couplings && syr_coupling == NULL;
       i++) {

    if (strcmp(_syr_couplings[i]->syr_name, syr_name) == 0) {
      syr_coupling = _syr_couplings[i];

      BFT_FREE(syr_coupling->syr_name);
      BFT_FREE(syr_coupling->b_location_ids);
      BFT_FREE(syr_coupling->v_location_ids);

      assert(syr_coupling->faces == NULL);  /* not built yet at this stage */
      assert(syr_coupling->cells == NULL);
    }
  }

  /* Allocate _cs_syr_coupling_t structure */

  if (syr_coupling == NULL) {
    BFT_REALLOC(_syr_couplings,
                _syr_n_couplings + 1, cs_syr_coupling_t *);
    BFT_MALLOC(syr_coupling, 1, cs_syr_coupling_t);

    _syr_couplings[_syr_n_couplings] = syr_coupling;

    _syr_n_couplings++;
  }

  syr_coupling->dim = dim;
  syr_coupling->ref_axis = ref_axis;

  syr_coupling->syr_name = NULL;

  if (syr_name != NULL) {
    BFT_MALLOC(syr_coupling->syr_name, strlen(syr_name) + 1, char);
    strcpy(syr_coupling->syr_name, syr_name);
  }
  else {
    BFT_MALLOC(syr_coupling->syr_name, 1, char);
    syr_coupling->syr_name[0] = '\0';
  }

  /* Selection criteria  */

  syr_coupling->n_b_locations = 0;
  syr_coupling->n_v_locations = 0;
  syr_coupling->b_location_ids = NULL;
  syr_coupling->v_location_ids = NULL;

  syr_coupling->faces = NULL;
  syr_coupling->cells = NULL;

  syr_coupling->allow_nearest = allow_nonmatching;
  syr_coupling->tolerance = tolerance;
  syr_coupling->verbosity = verbosity;
  syr_coupling->visualization = visualization;

  /* Initialize communicators */

#if defined(HAVE_MPI)

  syr_coupling->comm = MPI_COMM_NULL;
  syr_coupling->n_syr_ranks = 0;
  syr_coupling->syr_root_rank = -1;

#endif

  return  syr_coupling;
}

/*----------------------------------------------------------------------------
 * Add a mesh location to a syr_coupling_t structure.
 *
 * parameters:
 *   syr_coupling  <-- SYRTHES coupling structure
 *   location_id   <-- id of mesh location to add (boundary faces or cells)
 *----------------------------------------------------------------------------*/

static void
_add_mesh_location(cs_syr_coupling_t  *syr_coupling,
                   int                  location_id)
{
  cs_mesh_location_type_t l_type = cs_mesh_location_get_type(location_id);

  if (l_type == CS_MESH_LOCATION_BOUNDARY_FACES) {
    int i = syr_coupling->n_b_locations;
    syr_coupling->n_b_locations += 1;
    BFT_REALLOC(syr_coupling->b_location_ids, syr_coupling->n_b_locations, int);

    syr_coupling->b_location_ids[i] = location_id;
  }

  else if (l_type == CS_MESH_LOCATION_CELLS) {
    int i = syr_coupling->n_v_locations;
    syr_coupling->n_v_locations += 1;
    BFT_REALLOC(syr_coupling->v_location_ids, syr_coupling->n_v_locations, int);

    syr_coupling->v_location_ids[i] = location_id;
  }
}

/*----------------------------------------------------------------------------
 * Destroy cs_syr_coupling_t structures
 *----------------------------------------------------------------------------*/

static void
_syr_coupling_all_destroy(void)
{
  cs_lnum_t i_coupl;
  cs_syr_coupling_t *syr_coupling = NULL;

  if (_syr_n_couplings == 0)
    return;

  _all_comm_times();

  for (i_coupl = 0; i_coupl < _syr_n_couplings; i_coupl++) {

    syr_coupling = _syr_couplings[i_coupl];

    /* Free _cs_syr_coupling structure */

    BFT_FREE(syr_coupling->syr_name);
    BFT_FREE(syr_coupling->b_location_ids);
    BFT_FREE(syr_coupling->v_location_ids);

    if (syr_coupling->faces != NULL)
      _destroy_coupled_ent(&(syr_coupling->faces));
    if (syr_coupling->cells != NULL)
      _destroy_coupled_ent(&(syr_coupling->cells));

    /* Close communication */

    _finalize_comm(syr_coupling);

    BFT_FREE(syr_coupling);

  } /* End of loop on _syr_couplings */

  _syr_n_couplings = 0;
  BFT_FREE(_syr_couplings);

  bft_printf(_("\nStructures associated with SYRTHES coupling freed.\n"));
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Get name of SYRTHES coupling.
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *
 * returns:
 *   pointer to SYRTHES coupling name
 *----------------------------------------------------------------------------*/

static const char *
_syr_coupling_get_name(cs_syr_coupling_t  *syr_coupling)
{
  const char *retval = cs_empty_string;

  if (syr_coupling->syr_name != NULL)
    retval = syr_coupling->syr_name;

  return retval;
}

/*----------------------------------------------------------------------------
 * Initialize communicator for SYRTHES coupling
 *
 * parameters:
 *   syr_coupling  <-> SYRTHES coupling structure
 *   coupling_id   <-- id of this coupling (for log file message)
 *   syr_root_rank <-- SYRTHES root rank
 *   n_syr_ranks   <-- Number of ranks associated with SYRTHES
 *----------------------------------------------------------------------------*/

static void
_syr_coupling_init_comm(cs_syr_coupling_t  *syr_coupling,
                        int                 coupling_id,
                        int                 syr_root_rank,
                        int                 n_syr_ranks)
{
#if defined(HAVE_MPI)

  char  volume_flag = ' ';
  char  boundary_flag = ' ';
  char  conservativity_flag = '1';
  char  allow_nearest_flag = '1';
  char  op_name_send[32 + 1];
  char  op_name_recv[32 + 1];

  syr_coupling->n_syr_ranks = n_syr_ranks;
  syr_coupling->syr_root_rank = syr_root_rank;

  _init_comm(syr_coupling, coupling_id);

  /* Exchange coupling options */

  if (syr_coupling->n_b_locations > 0)
    boundary_flag = 'b';
  if (syr_coupling->n_v_locations > 0)
    volume_flag = 'v';
  if (_syr_coupling_conservativity == 0)
    conservativity_flag = '0';
  if (syr_coupling->allow_nearest == false)
    allow_nearest_flag = '0';

  snprintf(op_name_send, 32, "coupling:type:%c%c%c \2\2%c(%6.2g)",
           boundary_flag, volume_flag, conservativity_flag,
           allow_nearest_flag, (double)syr_coupling->tolerance);

  _exchange_sync(syr_coupling, op_name_send, op_name_recv);

  if (strncmp(op_name_recv, op_name_send, 16))
    bft_error
      (__FILE__, __LINE__, 0,
       _("========================================================\n"
         "   ** Incompatible SYRTHES coupling options:\n"
         "      ------------------------------------------------\n"
         "      code_saturne options: \"%s\"\n"
         "      SYRTHES options:      \"%s\"\n"
         "========================================================\n"),
       op_name_send, op_name_recv);

#endif
}

/*----------------------------------------------------------------------------
 * Define coupled mesh and send it to SYRTHES
 *
 * Optional post-processing output is also built at this stage.
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *----------------------------------------------------------------------------*/

static void
_syr_coupling_init_mesh(cs_syr_coupling_t  *syr_coupling)
{
  const cs_lnum_t verbosity = syr_coupling->verbosity;

  if (verbosity > 0)
    bft_printf(_("\n ** Processing the mesh for SYRTHES coupling "
                 "\"%s\"\n\n"),
                 syr_coupling->syr_name);

  /* Define coupled mesh */

  assert(syr_coupling->dim == 3 || syr_coupling->dim == 2);

  int match_flag = 0;

  if (syr_coupling->n_b_locations > 0) {
    syr_coupling->faces = _create_coupled_ent(syr_coupling,
                                              syr_coupling->n_b_locations,
                                              syr_coupling->b_location_ids,
                                              syr_coupling->dim - 1);
    match_flag += _sync_after_location(syr_coupling);
  }

  if (syr_coupling->n_v_locations > 0) {
    syr_coupling->cells = _create_coupled_ent(syr_coupling,
                                              syr_coupling->n_v_locations,
                                              syr_coupling->v_location_ids,
                                              syr_coupling->dim);
    match_flag += _sync_after_location(syr_coupling);
  }

  /* Communication with SYRTHES */
  /*----------------------------*/

  if (match_flag == 0 && verbosity > 0) {
    bft_printf(_("\n ** Mesh located for SYRTHES coupling \"%s\".\n\n"),
               syr_coupling->syr_name);
    bft_printf_flush();
  }
}

/*----------------------------------------------------------------------------
 * Return 1 if this coupling is a surface coupling else 0
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *
 * returns:
 *   1 or 0
 *----------------------------------------------------------------------------*/

static inline int
_syr_coupling_is_surf(const cs_syr_coupling_t  *syr_coupling)
{
  int retval = 0;

  if (syr_coupling->n_b_locations > 0)
    retval = 1;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return 1 if this coupling is a volume coupling else 0
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *
 * returns:
 *   1 or 0
 *----------------------------------------------------------------------------*/

static inline int
_syr_coupling_is_vol(const cs_syr_coupling_t  *syr_coupling)
{
  int retval = 0;

  if (syr_coupling->n_v_locations > 0)
    retval = 1;

  return retval;
}

/*----------------------------------------------------------------------------
 * Receive coupling variables from SYRTHES
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *   tsolid       --> solid temperature
 *   mode         <-- 0: surface coupling; 1: volume coupling
 *----------------------------------------------------------------------------*/

static void
_syr_coupling_recv_tsolid(cs_syr_coupling_t  *syr_coupling,
                          cs_real_t           tsolid[],
                          int                 mode)
{
  cs_syr_coupling_ent_t  *coupling_ent = NULL;

  assert(mode == 0 || mode == 1);

  if (mode == 0)
    coupling_ent = syr_coupling->faces;
  else
    coupling_ent = syr_coupling->cells;

  if (coupling_ent == NULL)
    return;

  /* Receive data */

  ple_locator_exchange_point_var(coupling_ent->locator,
                                 NULL,
                                 tsolid,
                                 NULL,
                                 sizeof(cs_real_t),
                                 1,
                                 0);

  if (coupling_ent->n_elts > 0) {
    if (mode == 1) { /* Save tsolid for a future used
                        in source term definition */
      cs_lnum_t i;
      assert(coupling_ent->solid_temp != NULL);
      for (i = 0; i < coupling_ent->n_elts; i++)
        coupling_ent->solid_temp[i] = tsolid[i];
    }
    else
      _post_var_update(coupling_ent, 0, tsolid);
  }
}

/*----------------------------------------------------------------------------
 * Send coupling variables to SYRTHES
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *   cpl_elt_ids  <-- ids of coupled elements
 *   tf           <-- fluid temperature
 *   hf           <-- fluid heat exchange coef. (numerical or user-defined)
 *   mode          <-- 0: surface coupling; 1: volume coupling
 *----------------------------------------------------------------------------*/

static void
_syr_coupling_send_tf_hf(cs_syr_coupling_t  *syr_coupling,
                         const cs_lnum_t     cpl_elt_ids[],
                         cs_real_t           tf[],
                         cs_real_t           hf[],
                         int                 mode)
{
  double *send_var = NULL;
  cs_syr_coupling_ent_t  *coupling_ent = NULL;

  assert(mode == 0 || mode == 1);

  if (mode == 0)
    coupling_ent = syr_coupling->faces;
  else
    coupling_ent = syr_coupling->cells;

  if (coupling_ent == NULL)
    return;

  const cs_lnum_t n_dist
    = ple_locator_get_n_dist_points(coupling_ent->locator);
  const cs_lnum_t *dist_loc
    = ple_locator_get_dist_locations(coupling_ent->locator);

  /* Prepare and send data */

  BFT_MALLOC(send_var, n_dist*2, double);

  for (cs_lnum_t ii = 0; ii < n_dist; ii++) {
    send_var[ii*2]     = tf[dist_loc[ii]];
    send_var[ii*2 + 1] = hf[dist_loc[ii]];
  }

  ple_locator_exchange_point_var(coupling_ent->locator,
                                 send_var,
                                 NULL,
                                 NULL,
                                 sizeof(double),
                                 2,
                                 0);

  BFT_FREE(send_var);

  if (mode == 1 && coupling_ent->n_elts > 0) {

    _post_var_update(coupling_ent, 1, tf);
    _post_var_update(coupling_ent, 2, hf);

    /* Saved hf for a future used in source term definition */

    assert(coupling_ent->hvol != NULL);
    for (cs_lnum_t ii = 0; ii < coupling_ent->n_elts; ii++)
      coupling_ent->hvol[ii] = hf[ii];

  }

  /* Exchange flux and corrector coefficient to ensure conservativity */

  if (_syr_coupling_conservativity > 0 && mode == 0)
    _ensure_conservativity(syr_coupling, cpl_elt_ids);
}

/*----------------------------------------------------------------------------
 * Compute the explicit/implicit contribution to source terms in case of
 * volume coupling with SYRTHES
 *
 * parameters:
 *   syr_coupling  <-- SYRTHES coupling structure
 *   tf            <-- fluid temperature
 *   ctbimp        <-> implicit contribution
 *   ctbexp        <-> explicit contribution
 *----------------------------------------------------------------------------*/

static void
_syr_coupling_ts_contrib(const cs_syr_coupling_t  *syr_coupling,
                         const cs_real_t           tf[],
                         cs_real_t                 ctbimp[],
                         cs_real_t                 ctbexp[])
{
  int  i;

  const double  *hvol = NULL;
  const cs_real_t  *solid_temp = NULL;
  cs_syr_coupling_ent_t  *ent = NULL;

  /* sanity checks */

  assert(syr_coupling != NULL);
  assert(syr_coupling->cells != NULL);

  ent = syr_coupling->cells;
  hvol = ent->hvol;
  solid_temp = ent->solid_temp;

  assert(hvol != NULL || ent->n_elts == 0);
  assert(solid_temp != NULL || ent->n_elts == 0);

  /* Compute contribution */

  if (_syr_coupling_implicit == 0) { /* Explicit treatment */

    for (i = 0; i < ent->n_elts; i++) {
      ctbexp[i] = -hvol[i] * (tf[i] - solid_temp[i]);
      ctbimp[i] = 0.0;
    }

  }
  else { /* Implicit treatment */

    for (i = 0; i < ent->n_elts; i++) {
      ctbexp[i] = hvol[i] * solid_temp[i];
      ctbimp[i] = hvol[i];
    }

  } /* Test if implicit */

}

/*----------------------------------------------------------------------------
 * Print information on yet unmatched SYRTHES couplings.
 *
 * parameters:
 *   n_unmatched    <--  number of unmatched couplings
 *   unmatched_ids  <--  array of unmatched couplings
 *----------------------------------------------------------------------------*/

static void
_print_all_unmatched_syr(int        n_unmatched,
                         const int  unmatched_ids[])
{
  /* Loop on defined SYRTHES instances */

  for (int i = 0; i < n_unmatched; i++) {

    cs_syr_coupling_t *syr_coupling
      = _syr_coupling_by_id(unmatched_ids[i]);
    const char *local_name = _syr_coupling_get_name(syr_coupling);

    bft_printf(_(" SYRTHES coupling:\n"
                 "   coupling id:              %d\n"
                 "   local name:               \"%s\"\n\n"),
               i, local_name);

  }

  bft_printf_flush();
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize MPI SYRTHES couplings using MPI.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *
 * parameters:
 *   n_unmatched    <->  pointer to number of unmatched couplings
 *   unmatched_ids  <->  pointer to array of unmatched couplings
 *----------------------------------------------------------------------------*/

static void
_init_all_mpi_syr(int  *n_unmatched,
                  int  **unmatched_ids)
{
  int _n_unmatched = *n_unmatched;
  int *_unmatched_ids = *unmatched_ids;

  const int n_couplings = _syr_n_couplings;

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == NULL)
    return;

  const int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);

  /* Loop on applications */

  for (int i = 0; i < n_apps; i++) {

    ple_coupling_mpi_set_info_t ai = ple_coupling_mpi_set_get_info(mpi_apps, i);

    if (strncmp(ai.app_type, "SYRTHES 4", 9) == 0) {

      int  match_queue_id = -1;
      int  coupling_id = -1;

      if (n_apps == 2 && n_couplings == 1 && _n_unmatched == 1) {
        match_queue_id = 0;
        coupling_id = 0;
      }
      else if (ai.app_name != NULL) {
        for (int j = 0; j < _n_unmatched; j++) {
          int k = _unmatched_ids[j];
          cs_syr_coupling_t *scpl = _syr_coupling_by_id(k);
          if (strcmp(ai.app_name, _syr_coupling_get_name(scpl)) == 0) {
            coupling_id = k;
            match_queue_id = j;
            break;
          }
        }
      }

      if (coupling_id > -1) {

        /* Remove from unmatched queue */
        _n_unmatched -= 1;
        for (int l = match_queue_id; l < _n_unmatched; l++)
          _unmatched_ids[l] = _unmatched_ids[l+1];
        if (_n_unmatched == 0)
          BFT_FREE(_unmatched_ids);

        /* Set communicator */
        _syr_coupling_init_comm(_syr_coupling_by_id(coupling_id),
                                coupling_id,
                                ai.root_rank,
                                ai.n_ranks);

        /* Print matching info */

        const char *syr_version = cs_empty_string;
        const char *local_name = cs_empty_string;
        const char *distant_name = cs_empty_string;

        if (ai.app_name != NULL)
          local_name = ai.app_name;
        if (ai.app_type != NULL)
          syr_version = ai.app_type;
        if (ai.app_name != NULL)
          distant_name = ai.app_name;

        bft_printf(_(" SYRTHES coupling:\n"
                     "   coupling id:              %d\n"
                     "   version:                  \"%s\"\n"
                     "   local name:               \"%s\"\n"
                     "   distant application name: \"%s\"\n"
                     "   MPI application id:       %d\n"
                     "   MPI root rank:            %d\n"
                     "   number of MPI ranks:      %d\n\n"),
                   coupling_id, syr_version, local_name, distant_name,
                   i, ai.root_rank, ai.n_ranks);
      }

      /* Note that a SYRTHES app may be present in the coupling set, but
         not coupled to the current code_saturne instance, so
         coupling_id < 0 here should not be reported as an error or
         complained about here. In case if missing matches, only the
         codes having defined and missing couplings should complain. */

    }

  } /* End of loop on applications */

  bft_printf_flush();

  /* Set return values */

  *n_unmatched = _n_unmatched;
  *unmatched_ids = _unmatched_ids;
}

/*----------------------------------------------------------------------------
 * Find name of single SYRTHES coupling using MPI.
 *
 * If no coupling or multiple couplings are present, the default cannot be
 * determined, so NULL is returned.
 *----------------------------------------------------------------------------*/

static const char *
_mpi_syr_default_name(void)
{
  const char *retval = NULL;

  int n_syr_apps = 0;

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == NULL)
    return NULL;

  int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);

  /* First pass to count available SYRTHES couplings */

  for (int i = 0; i < n_apps; i++) {
    const ple_coupling_mpi_set_info_t
      ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
    if (strncmp(ai.app_type, "SYRTHES 4", 9) == 0) {
      if (n_syr_apps == 0)
        retval = ai.app_name;
      else
        retval = NULL;
      n_syr_apps += 1;
    }
  }

  return retval;
}

#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define new SYRTHES coupling.
 *
 * \param[in] syrthes_name      matching SYRTHES application name
 * \param[in] boundary_criteria surface selection criteria, or NULL
 * \param[in] volume_criteria   volume selection criteria, or NULL
 * \param[in] projection_axis   x', 'y', or 'y' for 2D projection axis (case
 *                              independent), or ' ' for standard 3D coupling
 * \param[in] allow_nonmatching allow nearest-neighbor mapping where matching
 *                              within tolerance is not available (useful
 *                              when meshes have a different level of detail)
 * \param[in] tolerance         addition to local extents of each element
 *                              extent = base_extent * (1 + tolerance)
 * \param[in] verbosity         verbosity level
 * \param[in] visualization     visualization output level (0 or 1)
 *
 * In the case of a single code_saturne and single SYRTHES instance, the
 * 'syrthes_name' argument is ignored, as there is only one matching
 * possibility.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances based on the 'syrthes_name' argument.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_define(const char  *syrthes_name,
                       const char  *boundary_criteria,
                       const char  *volume_criteria,
                       char         projection_axis,
                       bool         allow_nonmatching,
                       float        tolerance,
                       int          verbosity,
                       int          visualization)
{
  int dim = 3;
  int ref_axis = -1;

  switch (projection_axis) {
  case 'x':
  case 'X':
    dim = 2;
    ref_axis = 0;
    break;
  case 'y':
  case 'Y':
    dim = 2;
    ref_axis = 1;
    break;
  case 'z':
  case 'Z':
    dim = 2;
    ref_axis = 2;
    break;
  default:
    break;
  }

  /* Ensure name is available */

#if defined(HAVE_MPI)
  if (syrthes_name == NULL)
    syrthes_name = _mpi_syr_default_name();
#endif

  if (syrthes_name == NULL)
    syrthes_name = cs_empty_string;

  /* Define additional coupling */

  cs_syr_coupling_t  *syr_coupling = _syr_coupling_define(dim,
                                                          ref_axis,
                                                          syrthes_name,
                                                          allow_nonmatching,
                                                          tolerance,
                                                          verbosity,
                                                          visualization);

  /* Add locations if done at that stage (deprecated) */

  int n_locations = cs_mesh_location_n_locations();

  const char *sel_criteria[2] = {boundary_criteria, volume_criteria};
  const char *type_name[2] = {"faces", "cells"};
  cs_mesh_location_type_t type_filter[2] = {CS_MESH_LOCATION_BOUNDARY_FACES,
                                            CS_MESH_LOCATION_CELLS};

  for (int i = 0; i < 2; i++) {

    if (sel_criteria[i] != NULL) {
      for (int j = 0; j < n_locations && sel_criteria[i] != NULL; j++) {
        cs_mesh_location_type_t l_type = cs_mesh_location_get_type(j);
        if (l_type & type_filter[i]) {
          const char *c = cs_mesh_location_get_selection_string(j);
          if (c != NULL) {
            if (strcmp(c, sel_criteria[i]) == 0) {
              _add_mesh_location(syr_coupling, j);
              sel_criteria[i] = NULL;
            }
          }
        }
      }
    }

    if (sel_criteria[i] != NULL) {

      char *name;
      size_t l = strlen(syrthes_name) + strlen(type_name[i]) + 2;
      BFT_MALLOC(name, l, char);
      snprintf(name, l, "%s_%s", syrthes_name, type_name[i]);

      int j = cs_mesh_location_add(name, type_filter[i], sel_criteria[i]);

      BFT_FREE(name);

      _add_mesh_location(syr_coupling, j);

    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associated a zone to a defined SYRTHES coupling.
 *
 * \param[in] syrthes_name  matching SYRTHES application name
 * \param[in] z             pointer to matching zone
 *                          (boundary or volume)
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_add_zone(const char       *syrthes_name,
                         const cs_zone_t  *z)
{
  /* Ensure name is available */

#if defined(HAVE_MPI)
  if (syrthes_name == NULL)
    syrthes_name = _mpi_syr_default_name();
#endif

  if (syrthes_name == NULL)
    syrthes_name = cs_empty_string;

  /* Search for matching name in existing couplings */

  int n_couplings = _syr_n_couplings;
  bool match = false;

  for (int i = 0; i < n_couplings; i++) {

    cs_syr_coupling_t  *syr_coupling = _syr_coupling_by_id(i);
    const char *cmp_name = _syr_coupling_get_name(syr_coupling);

    if (strcmp(syrthes_name, cmp_name) == 0) {
      _add_mesh_location(syr_coupling, z->location_id);
      match = true;
      break;
    }

  }

  if (match == false)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: no defined SYRTHES coupling named \"%s\"."),
              __func__, syrthes_name);
}

/*----------------------------------------------------------------------------
 * Initialize SYRTHES couplings.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_init(void)
{
  int n_unmatched = _syr_n_couplings;

  int *unmatched_ids;
  BFT_MALLOC(unmatched_ids, n_unmatched, int);

  for (int i = 0; i < n_unmatched; i++)
    unmatched_ids[i] = i;

  /* First try using MPI */

#if defined(HAVE_MPI)

  if (n_unmatched > 0)
    _init_all_mpi_syr(&n_unmatched, &unmatched_ids);

#endif

  if (n_unmatched > 0) {

    bft_printf("Unmatched SYRTHES couplings:\n"
               "----------------------------\n\n");

    _print_all_unmatched_syr(n_unmatched, unmatched_ids);

    BFT_FREE(unmatched_ids);

    bft_error(__FILE__, __LINE__, 0,
              _("At least 1 SYRTHES coupling was defined for which\n"
                "no communication with a SYRTHES instance is possible."));
  }
}

/*----------------------------------------------------------------------------
 * Finalize all SYRTHES couplings.
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_finalize(void)
{
  _syr_coupling_all_destroy();
}

/*----------------------------------------------------------------------------
 * Return number of SYRTHES couplings.
 *
 * return:
 *   number of SYRTHES couplings defined
 *----------------------------------------------------------------------------*/

int
cs_syr_coupling_n_couplings(void)
{
  return _syr_n_couplings;
}

/*----------------------------------------------------------------------------
 * Set conservativity forcing flag to True (1) or False (0) for all defined
 * SYRTHES couplings
 *
 * parameter:
 *   flag     <--  Conservativity forcing flag to set
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_set_conservativity(int  flag)
{
  assert(flag == 0 || flag == 1);
  _syr_coupling_conservativity = flag;
}

/*----------------------------------------------------------------------------
 * Set explicit treatment for the source terms in SYRTHES volume couplings
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_set_explicit_treatment(void)
{
  _syr_coupling_implicit = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log SYRTHES coupling setup information.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_log_setup(void)
{
  /* Get the number of SYRTHES couplings */
  int n_coupl = _syr_n_couplings;
  const int keysca = cs_field_key_id("scalar_id");
  const int kcpsyr = cs_field_key_id("syrthes_coupling");
  int icpsyr;

  if (n_coupl >= 1) {

    cs_log_printf
      (CS_LOG_SETUP,
       _("SYRTHES coupling\n"
         "----------------\n\n"
         "    number of couplings: %d\n"),
         n_coupl);

    int n_surf_coupl = 0, n_vol_coupl = 0, issurf, isvol;

    for (int coupl_id = 0; coupl_id < n_coupl; coupl_id++) {
      cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(coupl_id);

      /* Add a new surface coupling if detected */
      issurf = _syr_coupling_is_surf(syr_coupling);
      n_surf_coupl += issurf;

      /* Add a new volume coupling if detected */
      isvol = _syr_coupling_is_vol(syr_coupling);
      n_vol_coupl += isvol;
    }

    cs_log_printf
      (CS_LOG_SETUP,
       _("    with             %d surface coupling(s)\n"
         "    with             %d volume coupling(s)\n"),
         n_surf_coupl, n_vol_coupl);

    cs_log_printf
      (CS_LOG_SETUP,
       _("\n"
         "   Coupled scalars\n"
         "------------------------\n"
         " Scalar    Number icpsyr\n"
         "------------------------\n"));


    for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
      cs_field_t  *f = cs_field_by_id(f_id);
      if ((f->type & CS_FIELD_VARIABLE) || (f->type & CS_FIELD_USER)) {
        int ii = cs_field_get_key_int(f, keysca);
        if (ii > 0) {
          icpsyr = cs_field_get_key_int(f, kcpsyr);
          cs_log_printf
            (CS_LOG_SETUP,
             _(" %s %7d %7d\n"),cs_field_get_label(f),ii, icpsyr);
        }
      }
    }
    cs_log_printf
      (CS_LOG_SETUP,
       _("------------------------\n\n"
         "    icpsyr = 0 or 1         (1: scalar coupled to SYRTHES)\n"));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create coupled meshes and setup PLE locator for Syrthes couplings.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_init_meshes(void)
{
  int n_coupl = _syr_n_couplings;

  for (int coupl_id = 0; coupl_id < n_coupl; coupl_id++) {
    cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(coupl_id);
    _syr_coupling_init_mesh(syr_coupling);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the given SYRTHES coupling number is a surface couplings.
 *
 * \param[in] cpl_id   matching SYRTHES coupling id
 *
 * \return 1 if the coupling includes the surface, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_syr_coupling_is_surf(int  cpl_id)
{
  int retval = 0;  /* Default initialization */

  cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(cpl_id);

  if (syr_coupling == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling id %d impossible; "
                "there are %d couplings"),
              cpl_id, _syr_n_couplings);

  retval = _syr_coupling_is_surf(syr_coupling);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read boundary field/variable values relative to a SYRTHES coupling.
 *
 * \param[in]       nvar     number of variables
 * \param[in]       bc_type  boundary condition type
 * \param[in, out]  icodcl   boundary condition codes
 * \param[in, out]  rcodcl   boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_recv_boundary(int        nvar,
                              int        bc_type[],
                              int        icodcl[],
                              cs_real_t  rcodcl[])
{
  /* SYRTHES coupling: get wall temperature
     ====================================== */

  const int kcpsyr = cs_field_key_id("syrthes_coupling");

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t n_vars = nvar;  /* cast to cs_lnum_t because this
                                     is used in address computations here */

  /* Get number of coupling cases */

  int n_cpl = cs_syr_coupling_n_couplings();

  /* Loop on fields, handling only those coupled with Syrthes */

  int n_fields = cs_field_n_fields();
  for (int field_id = 0 ; field_id <  n_fields; field_id++) {

    cs_field_t  *f = cs_field_by_id(field_id);

    int icpsyr = 0;
    if (f->type & CS_FIELD_VARIABLE)
      icpsyr = cs_field_get_key_int(f, kcpsyr);

    if (icpsyr < 1)
      continue;

    /* Loop on couplings: get wall temperature array for each coupling
       and apply matching boundary condition. */

    for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {

      cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(cpl_id);
      cs_syr_coupling_ent_t *coupling_ent = syr_coupling->faces;

      if (coupling_ent == NULL)  /* ignore if volume-only */
        continue;

      cs_lnum_t n_cpl_faces = coupling_ent->n_elts;

      /* Get list of coupled faces */

      cs_lnum_t  *f_ids;
      BFT_MALLOC(f_ids, n_cpl_faces, cs_lnum_t);
      fvm_nodal_get_parent_id(coupling_ent->elts,
                              coupling_ent->elt_dim,
                              f_ids);

      /* Read wall temperature and interpolate if necessary */

      cs_real_t *t_solid;
      BFT_MALLOC(t_solid, n_cpl_faces, cs_real_t);
      _syr_coupling_recv_tsolid(syr_coupling, t_solid, 0);

      /*  For scalars coupled with SYRTHES, prescribe a Dirichlet
          condition at coupled faces.
          For the time being, pass here only once, as only one scalar is
          coupled with SYRTHES.
          For the compressible module, solve in energy, but save the
          temperature separately, for BC's to be clearer. */

      const int k_var_id = cs_field_key_id("variable_id");
      int var_id = cs_field_get_key_int(f, k_var_id) - 1;

      if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0) {
        if (f == CS_F_(e_tot)) {
          const cs_field_t *f_t_kelvin = CS_F_(t_kelvin);
          var_id = cs_field_get_key_int(f_t_kelvin, k_var_id);
        }
        else
          bft_error
            (__FILE__, __LINE__, 0,
             _("With the compressible module, only the \"total energy\"\n"
               "scalar field may be coupled with SYRTHES.\n"
               "Here, one tries to couple with the field \"%s\"."),
             f->name);
      }

      int  *_icodcl = icodcl + (var_id*n_b_faces);
      cs_real_t  *_rcodcl1 = rcodcl + (var_id*n_b_faces);
      cs_real_t  *_rcodcl2 = rcodcl + (n_b_faces*n_vars + var_id*n_b_faces);
      cs_real_t  *_rcodcl3 = rcodcl + (2*n_b_faces*n_vars + var_id*n_b_faces);

      for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {

        cs_lnum_t face_id = f_ids[i];

        if (   _icodcl[face_id] != CS_INDEF
            && _icodcl[face_id] != CS_SMOOTHWALL
            && _icodcl[face_id] != CS_ROUGHWALL) {
          if (bc_type[face_id] == CS_SMOOTHWALL)
            _icodcl[face_id] = CS_SMOOTHWALL;
          else if (bc_type[face_id] == CS_ROUGHWALL)
            _icodcl[face_id] = CS_ROUGHWALL;
        }

        _rcodcl1[face_id] = t_solid[i];
        _rcodcl2[face_id] = cs_math_infinite_r;
        _rcodcl3[face_id] = 0.;

      }

      /* Require temperature -> enthalpy conversion */

      if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

        if (f == cs_thermal_model_field()) {
          for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
            cs_lnum_t face_id = f_ids[i];
            _icodcl[face_id] *= -1;
          }
        }

      } /* End case for enthalpy */

      BFT_FREE(f_ids);
      BFT_FREE(t_solid);

    } /* End loop on couplings */

  } /* End loop on fields */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Send field/variable values relative to a SYRTHES coupling.
 *
 * \param[in]  h_wall   wall thermal exchange coefficient
 * \param[in]  v_fluid  near-wall fluid thermal variable
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_send_boundary(const cs_real_t  h_wall[],
                              cs_real_t        v_fluid[])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  /* Get number of coupling cases */

  int n_cpl = cs_syr_coupling_n_couplings();

  /* Check if we have a boundary coupling. */

  bool have_boundary_cpl = false;
  for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {
    cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(cpl_id);
    if (_syr_coupling_is_surf(syr_coupling)) {
      have_boundary_cpl = true;
      break;
    }
  }

  if (! have_boundary_cpl)
    return;

  /* Build arrays large enough for all cases */

  cs_lnum_t  *f_ids;
  cs_real_t  *t_fluid, *h_cpl;
  BFT_MALLOC(f_ids, n_b_faces, cs_lnum_t);
  BFT_MALLOC(t_fluid, n_b_faces, cs_real_t);
  BFT_MALLOC(h_cpl, n_b_faces, cs_real_t);

  /* Prepare conversion to temperature for enthalpy or energy
     (check for surface couplings to make sure it is needed,
     exit earlier otherwise) */

  cs_real_t  *wa = NULL;
  if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {
    BFT_MALLOC(wa, n_b_faces, cs_real_t);
    cs_ht_convert_h_to_t_faces(v_fluid, wa);
  }
  else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TOTAL_ENERGY) {
    /* Epsilon sup for perfect gas at cells */
    BFT_MALLOC(wa, n_cells, cs_real_t);
    cs_cf_thermo_eps_sup(CS_F_(rho)->val, wa, n_cells);
  }

  /* Loop on couplings */

  for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {

    cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(cpl_id);
    cs_syr_coupling_ent_t *coupling_ent = syr_coupling->faces;

    if (coupling_ent == NULL)  /* ignore if volume-only */
      continue;

    cs_lnum_t n_cpl_faces = coupling_ent->n_elts;

    /* Get list of coupled faces */

    fvm_nodal_get_parent_id(coupling_ent->elts,
                            coupling_ent->elt_dim,
                            f_ids);

    switch (cs_glob_thermal_model->itherm) {

    case CS_THERMAL_MODEL_TEMPERATURE:
      {
        for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
          cs_lnum_t face_id = f_ids[i];

          /* Saved fluid temperatures and exchange coefficients */
          t_fluid[i] = v_fluid[face_id];
          h_cpl[i] = h_wall[face_id];
        }
      }
      break;

    case CS_THERMAL_MODEL_ENTHALPY:
      {
        /* In enthalpy formulation, transform to temperatures for SYRTHES
         *  To conserve flux Phi = (lambda/d     ) Delta T
         *                 or Phi = (lambda/(d Cp)) Delta H
         * recall      hbord = lambda/d.
         *  Conservation is not guaranteed, so we add a warning. */

        for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
          cs_lnum_t face_id = f_ids[i];

          t_fluid[i] = wa[face_id];
          h_cpl[i] = h_wall[face_id];
        }
      }
      break;

    case CS_THERMAL_MODEL_TOTAL_ENERGY:
      {
        /* In energy formulation, transform to temperatures for SYRTHES
         *  To conserve flux Phi = (lambda/d     ) Delta T
         *                or Phi = (lambda/(d Cp)) Delta H
         *  Recall      hbord = lambda/ d
         *  Note that Ei = Cv Ti + 1/2 Ui*Ui + Epsilon_sup_i
         *  and  that Ep = Cv Tp + 1/2 Ui*Ui + Epsilon_sup_i
         *    (the difference is thus Cv Delta T)

         * Modify temperature and exchange coefficient

         * Compute e - CvT */

        const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

        cs_real_t   cv0 = cs_glob_fluid_properties->cv0;
        const cs_real_t  *cv = NULL;
        cs_lnum_t   cv_step = 0;

        const cs_real_3_t *cvar_vel = (const cs_real_3_t *)CS_F_(vel)->val;

        if (CS_F_(cv) != NULL) {
          cv = (const cs_real_t *)CS_F_(cv)->val;
          cv_step = 1;
        }
        else
          cv = &cv0;

        const cs_real_t *cvar_e_tot = (const cs_real_t *)CS_F_(e_tot)->val;

        for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
          cs_lnum_t face_id = f_ids[i];
          cs_lnum_t cell_id = b_face_cells[face_id];

          cs_real_t cvt =   cvar_e_tot[face_id]
                          - 0.5*cs_math_3_square_norm(cvar_vel[cell_id])
                          + wa[cell_id];

          t_fluid[i] = cvt / cv[cell_id * cv_step];
          h_cpl[i] = h_wall[face_id];
        }
      }
      break;

    default:
      break;
    }

    /* Fluxes are multiplied by porosity if present.
     * Here as the flux is expressed as h.(Tw-Tf), the exchange coefficient
     * is multipled by the porosity. */

    const cs_field_t *f_poro = CS_F_(poro);

    if (f_poro != NULL) {

      const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

      if (f_poro->dim == 1) {
        const cs_real_t * cpro_poro = (const cs_real_t *)f_poro->val;
        for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
          cs_lnum_t face_id = f_ids[i];
          cs_lnum_t cell_id = b_face_cells[face_id];
          h_cpl[i] *= cpro_poro[cell_id];
        }
      }
      else if (f_poro->dim == 6) {
        const cs_real_6_t * cpro_poro = (const cs_real_6_t *)f_poro->val;
        for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
          cs_lnum_t face_id = f_ids[i];
          cs_lnum_t cell_id = b_face_cells[face_id];
          /* TODO: using a product between the porosity
             and the boundary face normal would be more precise. */
          h_cpl[i] *= 1./3. * (  cpro_poro[cell_id][0]
                               + cpro_poro[cell_id][1]
                               + cpro_poro[cell_id][2]);
        }
      }

    }

    _syr_coupling_send_tf_hf(syr_coupling, f_ids, t_fluid, h_cpl, 0);

  } /* End loop on couplings */

  BFT_FREE(wa);

  BFT_FREE(f_ids);
  BFT_FREE(t_fluid);
  BFT_FREE(h_cpl);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Exchange volume values relative to a SYRTHES coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_exchange_volume(void)
{
  const int kcpsyr = cs_field_key_id("syrthes_coupling");

  /* Get number of coupling cases */

  int n_cpl = cs_syr_coupling_n_couplings();

  /* Loop on fields, handling only those coupled with Syrthes */

  int n_fields = cs_field_n_fields();
  for (int field_id = 0 ; field_id <  n_fields; field_id++) {

    cs_field_t  *f = cs_field_by_id(field_id);

    int icpsyr = 0;
    if (f->type & CS_FIELD_VARIABLE)
      icpsyr = cs_field_get_key_int(f, kcpsyr);

    if (icpsyr < 1)
      continue;

    /* Sanity check : only temperature is possible when doing a
       volume coupling with SYRTHES */
    if (f != cs_thermal_model_field())
      bft_error
        (__FILE__, __LINE__, 0,
         _("SYRTHES volume coupling possible only with temperature variable,\n"
           "not \"%s\"."),
         f->name);

    /* Loop on couplings: get wall temperature array for each coupling
       and apply matching boundary condition. */

    for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {

      cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(cpl_id);
      cs_syr_coupling_ent_t *coupling_ent = syr_coupling->cells;

      if (coupling_ent == NULL)  /* ignore if surface-only */
        continue;

      cs_lnum_t n_cpl_cells = coupling_ent->n_elts;

      /* Get list of coupled cells */

      cs_lnum_t  *c_ids;
      cs_real_t *t_fluid, *h_vol;
      BFT_MALLOC(c_ids, n_cpl_cells, cs_lnum_t);
      BFT_MALLOC(t_fluid, n_cpl_cells, cs_real_t);
      BFT_MALLOC(h_vol, n_cpl_cells, cs_real_t);

      fvm_nodal_get_parent_id(coupling_ent->elts,
                              coupling_ent->elt_dim,
                              c_ids);

      for (cs_lnum_t i = 0; i < n_cpl_cells; i++) {
        h_vol[i] = 0.;
      }

      /* Receive solid temperature.
       * This temperature is stored in a C structure for a future
       * use in source term definition. */

      _syr_coupling_recv_tsolid(syr_coupling, t_fluid, 1);

      const cs_real_t  *cvar_t = (const cs_real_t *)f->val;

      const char  *syrthes_name = _syr_coupling_get_name(syr_coupling);


      cs_user_syrthes_coupling_volume_h(cpl_id,
                                        syrthes_name,
                                        n_cpl_cells,
                                        c_ids,
                                        h_vol);

      for (cs_lnum_t i = 0; i < n_cpl_cells; i++)
        t_fluid[i] = cvar_t[c_ids[i]];

      _syr_coupling_send_tf_hf(syr_coupling, c_ids, t_fluid, h_vol, 1);

      BFT_FREE(c_ids);
      BFT_FREE(t_fluid);
      BFT_FREE(h_vol);

    } /* End loop on couplings */

  } /* End loop on fields */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the source term (implicit and/or explicit part) for a
 *         volume coupling with SYRTHES.
 *
 * \param[in]       field_id  field id
 * \param[in, out]  st_exp    explicit source term
 * \param[in, out]  st_imp    implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_volume_source_terms(int        field_id,
                                    cs_real_t  st_exp[],
                                    cs_real_t  st_imp[])
{
  cs_field_t  *f = cs_field_by_id(field_id);

  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  /* Get number of coupling cases */

  int n_cpl = cs_syr_coupling_n_couplings();

  /* Sanity check : only temperature is possible when doing a
     volume coupling with SYRTHES */
  if (f != cs_thermal_model_field())
    bft_error
      (__FILE__, __LINE__, 0,
       _("SYRTHES volume coupling possible only with temperature variable,\n"
         "not \"%s\"."),
       f->name);

  /* Loop on couplings: get wall temperature array for each coupling
     and apply matching boundary condition. */

  for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {

    cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(cpl_id);
    cs_syr_coupling_ent_t *coupling_ent = syr_coupling->cells;

    if (coupling_ent == NULL)  /* ignore if surface-only */
      continue;

    /* Get list of coupled cells */

    cs_lnum_t n_cpl_cells = coupling_ent->n_elts;

    cs_lnum_t  *c_ids;
    cs_real_t *t_fluid, *ctbimp, *ctbexp;
    BFT_MALLOC(c_ids, n_cpl_cells, cs_lnum_t);
    BFT_MALLOC(t_fluid, n_cpl_cells, cs_real_t);
    BFT_MALLOC(ctbimp, n_cpl_cells, cs_real_t);
    BFT_MALLOC(ctbexp, n_cpl_cells, cs_real_t);

    fvm_nodal_get_parent_id(coupling_ent->elts,
                            coupling_ent->elt_dim,
                            c_ids);

    /* Compute implicit and explicit contribution to source terms */

    const cs_real_t *cvara_vart = (const cs_real_t *)f->vals[1];

    for (cs_lnum_t i = 0; i < n_cpl_cells; i++) {
      t_fluid[i] = cvara_vart[c_ids[i]];
    }

    _syr_coupling_ts_contrib(syr_coupling, t_fluid, ctbimp, ctbexp);

    /* Loop on coupled cells to compute crvexp and crvimp */

    for (cs_lnum_t i = 0; i < n_cpl_cells; i++) {

      cs_lnum_t c_id = c_ids[i];

      cs_real_t tsexp = (ctbexp[i] - ctbimp[i]*t_fluid[i]) * cell_f_vol[c_id];
      cs_real_t tsimp =  ctbimp[i] * cell_f_vol[c_id];

      st_exp[c_id] += tsexp;
      st_imp[c_id] += tsimp;

    }

    BFT_FREE(c_ids);
    BFT_FREE(t_fluid);
    BFT_FREE(ctbimp);
    BFT_FREE(ctbexp);

  } /* End loop on couplings */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get number of coupled elements with SYRTHES.
 *
 * \param[in]   cpl_id  coupling id
 * \param[in]   mode    0 for boundary, 1 for volume
 *
 * \return  number of coupled elements for this coupling
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_syr_coupling_n_elts(int  cpl_id,
                       int  mode)
{
  cs_lnum_t retval = 0;

  cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(cpl_id);

  if (syr_coupling == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling id %d impossible; "
                "there are %d couplings"),
              cpl_id, _syr_n_couplings);

  else {

    cs_syr_coupling_ent_t  *coupling_ent = NULL;

    assert(mode == 0 || mode == 1);

    if (mode == 0)
      coupling_ent = syr_coupling->faces;
    else
      coupling_ent = syr_coupling->cells;

    if (coupling_ent != NULL)
      retval = coupling_ent->n_elts;

  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get local ids of elements coupled with SYRTHES
 *
 * \param[in]    cpl_id   coupling id
 * \param[in]    mode     0 for boundary, 1 for volume
 * \param[out]   elt_ids  ids of coupled elements (preallocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_elt_ids(int        cpl_id,
                        int        mode,
                        cs_lnum_t  elt_ids[])
{
  cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(cpl_id);

  if (syr_coupling == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling id %d impossible; "
                "there are %d couplings"),
              cpl_id, _syr_n_couplings);

  else {

    cs_syr_coupling_ent_t  *coupling_ent = NULL;

    assert(mode == 0 || mode == 1);

    if (mode == 0)
      coupling_ent = syr_coupling->faces;
    else
      coupling_ent = syr_coupling->cells;

    if (coupling_ent != NULL)
      fvm_nodal_get_parent_id(coupling_ent->elts,
                              coupling_ent->elt_dim,
                              elt_ids);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Receive coupling variables from SYRTHES.
 *
 * \param[in]    cpl_id   coupling id
 * \param[in]    mode     0 for boundary, 1 for volume
 * \param[out]   t_solid  solid temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_recv_tsolid(int        cpl_id,
                            int        mode,
                            cs_real_t  t_solid[])
{
  cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(cpl_id);

  if (syr_coupling == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling id %d impossible; "
                "there are %d couplings"),
              cpl_id, _syr_n_couplings);

  else
    _syr_coupling_recv_tsolid(syr_coupling, t_solid, mode);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Send coupling variables to SYRTHES.
 *
 * \param[in]    cpl_id   coupling id
 * \param[in]    mode     0 for boundary, 1 for volume
 * \param[in]    elt_ids  ids of coupled elements
 * \param[in]    t_fluid  fluid temperature
 * \param[in]    h_fluid  fluid exchange coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_send_tf_hf(int              cpl_id,
                           int              mode,
                           const cs_lnum_t  elt_ids[],
                           cs_real_t        t_fluid[],
                           cs_real_t        h_fluid[])
{
  cs_syr_coupling_t *syr_coupling = _syr_coupling_by_id(cpl_id);

  if (syr_coupling == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling id %d impossible; "
                "there are %d couplings"),
              cpl_id, _syr_n_couplings);

  else
    _syr_coupling_send_tf_hf(syr_coupling, elt_ids,
                             t_fluid, h_fluid, mode);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
