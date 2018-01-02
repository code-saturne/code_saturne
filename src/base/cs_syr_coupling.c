/*============================================================================
 * SYRTHES coupling
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

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#if defined(HAVE_MPI)
#include "cs_coupling.h"
#endif

#include "cs_log.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"
#include "cs_syr4_coupling.h"

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

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Structure associated with SYRTHES coupling */

typedef struct {

  int      match_id;       /* Id of matched application, -1 initially */
  int      dim;            /* Coupled mesh dimension */
  int      ref_axis;       /* Selected axis for edge extraction */
  char    *app_name;       /* Application name */
  char    *face_sel_c;     /* Face selection criteria */
  char    *cell_sel_c;     /* Cell selection criteria */
  bool     allow_nearest;  /* Allow nearest-neighbor beyond tolerance */
  float    tolerance;      /* Tolerance */
  int      verbosity;      /* Verbosity level */
  int      visualization;  /* Visualization level */
  int      conservativity; /* Conservativity forcing flag */

} _cs_syr_coupling_builder_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static int _cs_glob_n_syr_cp = -1;
static int _cs_glob_n_syr4_cp = -1;

static int                         _syr_coupling_builder_size = 0;
static _cs_syr_coupling_builder_t *_syr_coupling_builder = NULL;

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

  for (i = 0; i < _syr_coupling_builder_size; i++) {

    _cs_syr_coupling_builder_t *scb = _syr_coupling_builder + i;

    if (scb->match_id > -1) {
      if (scb->face_sel_c != NULL)
        BFT_FREE(scb->face_sel_c);
      if (scb->cell_sel_c != NULL)
        BFT_FREE(scb->cell_sel_c);
      if (scb->app_name != NULL)
        BFT_FREE(scb->app_name);
    }
  }

  /* Now, remove marked entries and resize */

  for (i = 0; i < _syr_coupling_builder_size; i++) {
    _cs_syr_coupling_builder_t *scb = _syr_coupling_builder + i;
    if (scb->match_id < 0) {
      *(_syr_coupling_builder + n_unmatched_entries) = *scb;
      n_unmatched_entries += 1;
    }
  }

  _syr_coupling_builder_size = n_unmatched_entries;

  BFT_REALLOC(_syr_coupling_builder,
              _syr_coupling_builder_size,
              _cs_syr_coupling_builder_t);
}

/*----------------------------------------------------------------------------
 * Print information on yet unmatched SYRTHES couplings.
 *----------------------------------------------------------------------------*/

static void
_print_all_unmatched_syr(void)
{
  int i;

  const char empty_string[] = "";

  /* Loop on defined SYRTHES instances */

  for (i = 0; i < _syr_coupling_builder_size; i++) {

    _cs_syr_coupling_builder_t *scb = _syr_coupling_builder + i;

    if (scb->match_id < 0) {

      const char *local_name = empty_string;

      if (scb->app_name != NULL)
        local_name = scb->app_name;

      bft_printf(_(" SYRTHES coupling:\n"
                   "   coupling id:              %d\n"
                   "   local name:               \"%s\"\n\n"),
                 i, local_name);
    }
  }

  bft_printf_flush();
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Add a SYRTHES 4 coupling using MPI.
 *
 * parameters:
 *   builder_id    <-- SYRTHES application id in coupling builder
 *   syr_root_rank <-- root rank associated with SYRTHES
 *   n_syr_ranks   <-- number of ranks associated with SYRTHES
 *----------------------------------------------------------------------------*/

static void
_syr4_add_mpi(int builder_id,
              int syr_root_rank,
              int n_syr_ranks)
{
  cs_syr4_coupling_t *syr_coupling = NULL;
  _cs_syr_coupling_builder_t *scb = _syr_coupling_builder + builder_id;

  /* Addition of coupling and association of communicator could
     be grouped together, but we prefer to keep the steps separate
     so that we could simply add SYRTHES 4 couplings directly (without
     resorting to a temporary builder), then match communications */

  cs_syr4_coupling_add(scb->dim,
                       scb->ref_axis,
                       scb->face_sel_c,
                       scb->cell_sel_c,
                       scb->app_name,
                       scb->allow_nearest,
                       scb->tolerance,
                       scb->verbosity,
                       scb->visualization);

  syr_coupling = cs_syr4_coupling_by_id(cs_syr4_coupling_n_couplings() - 1);

  cs_syr4_coupling_init_comm(syr_coupling,
                             builder_id,
                             syr_root_rank,
                             n_syr_ranks);
}

/*----------------------------------------------------------------------------
 * Print information on identified SYRTHES couplings using MPI.
 *
 * This function requires coupling_builder information, and must thus
 * be called before removing matched builder entries.
 *----------------------------------------------------------------------------*/

static void
_print_all_mpi_syr(void)
{
  int i;

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();
  const char empty_string[] = "";

  /* Loop on defined SYRTHES instances */

  for (i = 0; i < _syr_coupling_builder_size; i++) {

    _cs_syr_coupling_builder_t *scb = _syr_coupling_builder + i;

    if (scb->match_id > -1) {

      const char *syr_version = empty_string;
      const char *local_name = empty_string;
      const char *distant_name = empty_string;

      const ple_coupling_mpi_set_info_t
        ai = ple_coupling_mpi_set_get_info(mpi_apps, scb->match_id);

      if (scb->app_name != NULL)
        local_name = scb->app_name;
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
                 i, syr_version, local_name, distant_name,
                 scb->match_id, ai.root_rank, ai.n_ranks);
    }
  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Initialize MPI SYRTHES couplings using MPI.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

static void
_init_all_mpi_syr(void)
{
  int i;

  int n_apps = 0;
  int n_matched_apps = 0;
  int n_syr4_apps = 0;
  int syr_app_id = -1;

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == NULL)
    return;

  n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);

  /* First pass to count available SYRTHES couplings */

  for (i = 0; i < n_apps; i++) {
    const ple_coupling_mpi_set_info_t
      ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
    if (strncmp(ai.app_type, "SYRTHES 4", 9) == 0) {
      n_syr4_apps += 1;
      syr_app_id = i;
    }
  }

  /* In single-coupling mode, no identification necessary */

  if (n_syr4_apps == 1 && _syr_coupling_builder_size == 1) {

    ple_coupling_mpi_set_info_t ai
      = ple_coupling_mpi_set_get_info(mpi_apps, syr_app_id);

    _syr_coupling_builder->match_id = syr_app_id;

    BFT_REALLOC(_syr_coupling_builder->app_name, strlen(ai.app_name) + 1, char);
    strcpy(_syr_coupling_builder->app_name, ai.app_name);

    n_matched_apps += 1;

  }

  /* In multiple-coupling mode, identification is necessary */

  else {

    int j;
    ple_coupling_mpi_set_info_t ai;

    int n_syr_apps = 0;
    int *syr_appinfo = NULL;

    /* First, build an array of matched/unmatched SYRTHES applications, with
       2 entries per instance: matched indicator, app_id */

    BFT_MALLOC(syr_appinfo, n_syr4_apps*2, int);

    for (i = 0; i < n_apps; i++) {
      ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
      if (strncmp(ai.app_type, "SYRTHES 4", 9) == 0) {
        syr_appinfo[n_syr_apps*2] = 0;
        syr_appinfo[n_syr_apps*2 + 1] = i;
        n_syr_apps += 1;
      }
    }

    assert(n_syr_apps == n_syr4_apps);

    /* Loop on defined SYRTHES instances */

    for (i = 0; i < _syr_coupling_builder_size; i++) {

      _cs_syr_coupling_builder_t *scb = _syr_coupling_builder + i;

      /* First loop on available SYRTHES instances to match app_names */

      if (scb->app_name != NULL) {

        for (j = 0; j < n_syr_apps; j++) {

          if (syr_appinfo[j*2] != 0) /* Consider only unmatched applications */
            continue;

          ai = ple_coupling_mpi_set_get_info(mpi_apps, syr_appinfo[j*2 + 1]);
          if (ai.app_name == NULL)
            continue;

          if (strcmp(ai.app_name, scb->app_name) == 0) {
            scb->match_id = syr_appinfo[j*2 + 1];
            syr_appinfo[j*2] = i;
            n_matched_apps += 1;
            break;
          }

        }

      }

    } /* End of loop on defined SYRTHES instances */

    BFT_FREE(syr_appinfo);

  } /* End of test on single or multiple SYRTHES matching algorithm */

  /* Print matching info */

  _print_all_mpi_syr();

  /* Now initialize matched couplings */
  /*----------------------------------*/

  for (i = 0; i < _syr_coupling_builder_size; i++) {

    _cs_syr_coupling_builder_t *scb = _syr_coupling_builder + i;

    if (scb->match_id > -1) {
      const ple_coupling_mpi_set_info_t
        ai = ple_coupling_mpi_set_get_info(mpi_apps, scb->match_id);

      if (strncmp(ai.app_type, "SYRTHES 4", 9) == 0)
        _syr4_add_mpi(i, ai.root_rank, ai.n_ranks);
    }
  }

  /* Cleanup */

  _remove_matched_builder_entries();
}

#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get number of SYRTHES couplings.
 *
 * Fortran Interface:
 *
 * SUBROUTINE NBCSYR
 * *****************
 *
 * INTEGER          n_couplings     : <-- : number of SYRTHES couplings
 *----------------------------------------------------------------------------*/

void CS_PROCF(nbcsyr, NBCSYR)
(
 cs_lnum_t  *n_couplings
)
{
  *n_couplings = cs_syr_coupling_n_couplings();
}

/*----------------------------------------------------------------------------
 * Test if the given SYRTHES coupling number is a surface coupling
 * Return 1 if true else 0
 *
 * Fortran Interface:
 *
 * SUBROUTINE TSURSY
 * *****************
 *
 * INTEGER          cplnum     : <-- : number of the SYRTHES coupling
 * INTEGER          issurf     : --> : 1 if surface coupling else 0
 *----------------------------------------------------------------------------*/

void CS_PROCF(tsursy, TSURSY)
(
 cs_int_t  *const cplnum,
 cs_int_t  *issurf
)
{
  int n_couplings = 0;

  *issurf = 0; /* Default initialization */

  assert(_cs_glob_n_syr_cp > -1);

  if (_cs_glob_n_syr_cp == _cs_glob_n_syr4_cp) {

    n_couplings = _cs_glob_n_syr_cp;

    if (*cplnum < 1 || *cplnum > n_couplings)
      bft_error(__FILE__, __LINE__, 0,
                _("SYRTHES coupling number %d impossible; "
                  "there are %d couplings"), *cplnum, n_couplings);

    {
      cs_syr4_coupling_t *syr_coupling
        = cs_syr4_coupling_by_id(*cplnum - 1);

      *issurf = cs_syr4_coupling_is_surf(syr_coupling);
    }

  }
  else { /* Couplings are still defined in the builder structure */

    if (_syr_coupling_builder_size == _cs_glob_n_syr_cp) {

      _cs_syr_coupling_builder_t *scb = NULL;

      n_couplings = _cs_glob_n_syr_cp;

      if (*cplnum < 1 || *cplnum > n_couplings)
        bft_error(__FILE__, __LINE__, 0,
                  _("SYRTHES coupling number %d impossible; "
                    "there are %d couplings"), *cplnum, n_couplings);

      scb = _syr_coupling_builder + (*cplnum) - 1;

      if (scb->face_sel_c != NULL)
        *issurf = 1;
    }

  }

}

/*----------------------------------------------------------------------------
 * Test if the given SYRTHES coupling number is a volume coupling
 * Return 1 if true else 0
 *
 * Fortran Interface:
 *
 * SUBROUTINE TVOLSY
 * *****************
 *
 * INTEGER          cplnum     : <-- : number of the SYRTHES coupling
 * INTEGER          issurf     : --> : 1 if volume coupling else 0
 *----------------------------------------------------------------------------*/

void CS_PROCF(tvolsy, TVOLSY)
(
 cs_int_t  *const cplnum,
 cs_int_t  *isvol
)
{
  int n_couplings = 0;

  *isvol = 0; /* Default initialization */

  assert(_cs_glob_n_syr_cp > -1);

  if (_cs_glob_n_syr_cp == _cs_glob_n_syr4_cp) {

    n_couplings = _cs_glob_n_syr_cp;

    if (*cplnum < 1 || *cplnum > n_couplings)
      bft_error(__FILE__, __LINE__, 0,
                _("SYRTHES coupling number %d impossible; "
                  "there are %d couplings"), *cplnum, n_couplings);

    {
      cs_syr4_coupling_t *syr_coupling
        = cs_syr4_coupling_by_id(*cplnum - 1);

      *isvol = cs_syr4_coupling_is_vol(syr_coupling);
    }

  }
  else { /* Couplings are still defined in the builder structure */

    if (_syr_coupling_builder_size == _cs_glob_n_syr_cp) {

      _cs_syr_coupling_builder_t *scb = NULL;

      n_couplings = _cs_glob_n_syr_cp;

      if (*cplnum < 1 || *cplnum > n_couplings)
        bft_error(__FILE__, __LINE__, 0,
                  _("SYRTHES coupling number %d impossible; "
                    "there are %d couplings"), *cplnum, n_couplings);

      scb = _syr_coupling_builder + (*cplnum) - 1;

      if (scb->cell_sel_c != NULL)
        *isvol = 1;

    }

  }

}

/*----------------------------------------------------------------------------
 * Get number of coupled elements with SYRTHES.
 *
 * Fortran Interface:
 *
 * SUBROUTINE NBESYR
 * *****************
 *
 * INTEGER          coupl_num       : --> : coupling number
 * INTEGER          mode            : --> : 0 (surface); 1 (volume)
 * INTEGER          n_coupl_elts    : <-- : number of coupled elements
 *----------------------------------------------------------------------------*/

void CS_PROCF(nbesyr, NBESYR)
(
 const cs_int_t  *coupl_num,
 const cs_int_t  *mode,
       cs_int_t  *n_coupl_elts
)
{
  int n_couplings = _cs_glob_n_syr4_cp;

  if (*coupl_num < 1 || *coupl_num > n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *coupl_num, n_couplings);

  else {
    cs_syr4_coupling_t *syr_coupling
      = cs_syr4_coupling_by_id(*coupl_num - 1);
    *n_coupl_elts = cs_syr4_coupling_get_n_elts(syr_coupling, *mode);
  }
}

/*----------------------------------------------------------------------------
 * Get local numbering of coupled elements
 *
 * Fortran interface:
 *
 * SUBROUTINE LELTSY
 * *****************
 *
 * INTEGER      coupl_num       : --> : coupling number
 * INTEGER      mode            : --> : 0 (surface); 1 (volume)
 * INTEGER      coupl_elt_list  : <-- : list of coupled elements
 *----------------------------------------------------------------------------*/

void CS_PROCF(leltsy, LELTSY)
(
 const cs_int_t    *coupl_num,
 const cs_int_t    *mode,
       cs_lnum_t   *coupl_elt_list
)
{
  int n_couplings = _cs_glob_n_syr4_cp;

  if (*coupl_num < 1 || *coupl_num > n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *coupl_num, n_couplings);
  else {
    cs_syr4_coupling_t *syr_coupling
      = cs_syr4_coupling_by_id(*coupl_num - 1);
    cs_syr4_coupling_get_elt_list(syr_coupling, coupl_elt_list, *mode);
  }
}

/*----------------------------------------------------------------------------
 * Receive coupling variables from SYRTHES
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYI
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * INTEGER          MODE        : --> : 0 (surface); 1 (volume)
 * DOUBLE PRECISION TSOLID      : <-- : Solid temperature
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyi, VARSYI)
(
 cs_int_t   *numsyr,
 cs_int_t   *mode,
 cs_real_t  *tsolid
)
{
  int n_couplings = _cs_glob_n_syr4_cp;

  if (*numsyr < 1 || *numsyr > n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *numsyr, n_couplings);
  else {
    cs_syr4_coupling_t *syr_coupling
      = cs_syr4_coupling_by_id(*numsyr - 1);
    cs_syr4_coupling_recv_tsolid(syr_coupling, tsolid, *mode);
  }
}

/*----------------------------------------------------------------------------
 * Send coupling variables to SYRTHES
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYO
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * INTEGER          MODE        : --> : 0 (surface); 1 (volume)
 * INTEGER          LSTELT      : --> : List of coupled elements
 * DOUBLE PRECISION TFLUID      : --> : Fluid temperature
 * DOUBLE PRECISION HFLUID      : --> : Exchange coefficient
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyo, VARSYO)
(
 cs_int_t   *numsyr,
 cs_int_t   *mode,
 cs_int_t   *lstelt,
 cs_real_t  *tfluid,
 cs_real_t  *hfluid
)
{
  int n_couplings = _cs_glob_n_syr4_cp;

  if (*numsyr < 1 || *numsyr > n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *numsyr, n_couplings);
  else {
    cs_syr4_coupling_t *syr_coupling
      = cs_syr4_coupling_by_id(*numsyr - 1);
    cs_syr4_coupling_send_tf_hf(syr_coupling, lstelt, tfluid, hfluid, *mode);
  }
}

/*----------------------------------------------------------------------------
 * Compute the explicit/implicit contribution to source terms in case of
 * volume coupling with SYRTHES4
 *
 * Fortran Interface:
 *
 * SUBROUTINE CTBVSY
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * DOUBLE PRECISION TFLUID      : --> : Fluid temperature
 * DOUBLE PRECISION CTBIMP      : <-> : Implicit contribution
 * DOUBLE PRECISION CTBEXP      : <-> : Explicit contribution
 *----------------------------------------------------------------------------*/

void CS_PROCF (ctbvsy, CTBVSY)
(
 cs_int_t   *numsyr,
 cs_real_t  *tfluid,
 cs_real_t  *ctbimp,
 cs_real_t  *ctbexp
)
{
  int n_couplings = _cs_glob_n_syr4_cp;

  if (*numsyr < 1 || *numsyr > n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *numsyr, n_couplings);

  {
    cs_syr4_coupling_t *syr_coupling
      = cs_syr4_coupling_by_id(*numsyr - 1);

    cs_syr4_coupling_ts_contrib(syr_coupling, tfluid, ctbimp, ctbexp);
  }

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define new SYRTHES coupling.
 *
 * The arguments to \ref cs_syr_coupling_define are:
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
 * In the case of a single Code_Saturne and single SYRTHES instance, the
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
  int  conservativity = 1; /* Default value */

  _cs_syr_coupling_builder_t *scb = NULL;

  /* Add corresponding coupling to temporary SYRTHES couplings array */

  BFT_REALLOC(_syr_coupling_builder,
              _syr_coupling_builder_size + 1,
              _cs_syr_coupling_builder_t);

  scb = &(_syr_coupling_builder[_syr_coupling_builder_size]);

  scb->match_id = -1;

  switch (projection_axis) {
  case 'x':
  case 'X':
    scb->dim = 2;
    scb->ref_axis = 0;
    break;
  case 'y':
  case 'Y':
    scb->dim = 2;
    scb->ref_axis = 1;
    break;
  case 'z':
  case 'Z':
    scb->dim = 2;
    scb->ref_axis = 2;
    break;
  default:
    scb->dim = 3;
    scb->ref_axis = -1;
  }

  scb->app_name = NULL;
  if (syrthes_name != NULL) {
    BFT_MALLOC(scb->app_name, strlen(syrthes_name) + 1, char);
    strcpy(scb->app_name, syrthes_name);
  }

  scb->face_sel_c = NULL;
  if (boundary_criteria != NULL) {
    BFT_MALLOC(scb->face_sel_c, strlen(boundary_criteria) + 1, char);
    strcpy(scb->face_sel_c, boundary_criteria);
  }

  scb->cell_sel_c = NULL;
  if (volume_criteria != NULL) {
    BFT_MALLOC(scb->cell_sel_c, strlen(volume_criteria) + 1, char);
    strcpy(scb->cell_sel_c, volume_criteria);
  }

  scb->allow_nearest = allow_nonmatching;
  scb->tolerance = tolerance;
  scb->verbosity = verbosity;
  scb->visualization = visualization;
  scb->conservativity = conservativity;

  _syr_coupling_builder_size += 1;
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
  /* First try using MPI */

#if defined(HAVE_MPI)

  if (_syr_coupling_builder_size > 0)
    _init_all_mpi_syr();

#endif

  if (_syr_coupling_builder_size > 0) {

    bft_printf("Unmatched SYRTHES couplings:\n"
               "----------------------------\n\n");

    _print_all_unmatched_syr();

    bft_error(__FILE__, __LINE__, 0,
              _("At least 1 SYRTHES coupling was defined for which\n"
                "no communication with a SYRTHES instance is possible."));
  }

  _cs_glob_n_syr4_cp = cs_syr4_coupling_n_couplings();
}

/*----------------------------------------------------------------------------
 * Finalize all SYRTHES couplings.
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_finalize(void)
{
  cs_syr4_coupling_all_destroy();
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
  if (_cs_glob_n_syr_cp < 0) {
    if (_syr_coupling_builder_size > 0)
      _cs_glob_n_syr_cp = _syr_coupling_builder_size;
    else
      _cs_glob_n_syr_cp = cs_syr4_coupling_n_couplings();
  }

  return _cs_glob_n_syr_cp;
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
  cs_syr4_coupling_set_conservativity(flag);
}

/*----------------------------------------------------------------------------
 * Set explicit treatment for the source terms in SYRTHES volume couplings
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_set_explicit_treatment(void)
{
  cs_syr4_coupling_set_explicit_treatment();
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
  int n_coupl = cs_syr_coupling_n_couplings();
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

    for (int ii = 1 ; ii <= n_coupl ; ii++) {
      /* Add a new surface coupling if detected */
      issurf = 0;
      CS_PROCF(tsursy, TSURSY)(&ii, &issurf);
      n_surf_coupl += issurf;

      /* Add a new volume coupling if detected */
      isvol = 0;
      CS_PROCF(tvolsy, TVOLSY)(&ii, &isvol);
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
  for (int coupl_id = 0; coupl_id < _cs_glob_n_syr4_cp; coupl_id++) {
    cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(coupl_id);
    cs_syr4_coupling_init_mesh(syr_coupling);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
