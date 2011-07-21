/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * SYRTHES coupling
 *============================================================================*/

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
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
#include "cs_coupling.h"
#endif

#include "cs_parall.h"
#include "cs_prototypes.h"

#include "cs_syr3_coupling.h"
#include "cs_syr3_messages.h"
#include "cs_syr4_coupling.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_syr_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
  int      verbosity;      /* Verbosity level */
  int      visualization;  /* Visualization level */

} _cs_syr_coupling_builder_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static int _cs_glob_n_syr_cp = -1;
static int _cs_glob_n_syr3_cp = -1;
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
     so that when SYRTHES 3 support is removed, we will simply be
     able to add SYRTHES 4 couplings directly (without resorting
     to a temporary builder), then match communications */

  cs_syr4_coupling_add(scb->dim,
                       scb->ref_axis,
                       scb->face_sel_c,
                       scb->cell_sel_c,
                       scb->app_name,
                       scb->verbosity,
                       scb->visualization);

  syr_coupling = cs_syr4_coupling_by_id(cs_syr4_coupling_n_couplings() - 1);

  cs_syr4_coupling_init_comm(syr_coupling,
                             builder_id,
                             syr_root_rank,
                             n_syr_ranks);
}

/*----------------------------------------------------------------------------
 * Add a SYRTHES 3 coupling using MPI.
 *
 * parameters:
 *   builder_id <-- SYRTHES application id in coupling builder
 *   syr_rank   <-- rank associated with SYRTHES
 *----------------------------------------------------------------------------*/

static void
_syr3_add_mpi(int builder_id,
              int syr_rank)
{
  cs_syr3_coupling_t *syr_coupling = NULL;
  _cs_syr_coupling_builder_t *scb = _syr_coupling_builder + builder_id;

  cs_syr3_coupling_add(scb->dim,
                       scb->ref_axis,
                       scb->face_sel_c,
                       scb->app_name,
                       syr_rank,
                       CS_SYR3_COMM_TYPE_MPI,
                       scb->verbosity,
                       scb->visualization);

  syr_coupling = cs_syr3_coupling_by_id(cs_syr3_coupling_n_couplings() - 1);

  cs_syr3_coupling_init_comm(syr_coupling, builder_id);
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
  int n_syr3_apps = 0, n_syr4_apps = 0;
  int syr_app_id = -1, syr_app_type = 0;

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
      syr_app_type = 4;
    }
    else if (strncmp(ai.app_type, "SYRTHES 3", 9) == 0) {
      n_syr3_apps += 1;
      syr_app_id = i;
      syr_app_type = 3;
    }
  }

  /* In single-coupling mode, no identification necessary */

  if ((n_syr3_apps + n_syr4_apps == 1) && _syr_coupling_builder_size == 1) {
    _syr_coupling_builder->match_id = syr_app_id;
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

    BFT_MALLOC(syr_appinfo, (n_syr3_apps + n_syr4_apps)*2, int);

    for (i = 0; i < n_apps; i++) {
      ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
      if (   (strncmp(ai.app_type, "SYRTHES 4", 9) == 0)
          || (strncmp(ai.app_type, "SYRTHES 3", 9) == 0)) {
        syr_appinfo[n_syr_apps*2] = 0;
        syr_appinfo[n_syr_apps*2 + 1] = i;
        n_syr_apps += 1;
      }
    }

    assert(n_syr_apps == (n_syr3_apps + n_syr4_apps));

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

      else if (strncmp(ai.app_type, "SYRTHES 3", 9) == 0)
        _syr3_add_mpi(i, ai.root_rank);
    }
  }

  /* Cleanup */

  _remove_matched_builder_entries();
}

#endif /* defined(HAVE_MPI) */

#if defined(HAVE_SOCKET)

/*----------------------------------------------------------------------------
 * Add a SYRTHES 3 coupling using sockets.
 *
 * parameters:
 *   builder_id <-- SYRTHES application id in coupling builder
 *----------------------------------------------------------------------------*/

static void
_syr3_add_socket(int builder_id)
{
  cs_syr3_coupling_t *syr_coupling = NULL;
  _cs_syr_coupling_builder_t *scb = _syr_coupling_builder + builder_id;

  cs_syr3_coupling_add(scb->dim,
                       scb->ref_axis,
                       scb->face_sel_c,
                       scb->app_name,
                       -1,
                       CS_SYR3_COMM_TYPE_SOCKET,
                       scb->verbosity,
                       scb->visualization);

  syr_coupling = cs_syr3_coupling_by_id(cs_syr3_coupling_n_couplings() - 1);

  cs_syr3_coupling_init_comm(syr_coupling, builder_id);
}

/*----------------------------------------------------------------------------
 * Initialize SYRTHES couplings using sockets.
 *
 * Currently, only SYRTHES 3.4 may be coupled using sockets
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *
 * parameters:
 *   port_num <-- port number (only used for rank 0; automatic on others)
 *----------------------------------------------------------------------------*/

static void
_init_all_socket_syr(int port_num)
{
  int i;

  /* Initialize socket server */

  if (_syr_coupling_builder_size > 0)
    cs_syr3_comm_init_socket(port_num);

  bft_printf
    ("SYRTHES couplings for which the socket interface will be used:\n"
     "--------------------------------------------------------------\n\n");

  _print_all_unmatched_syr();

  /* Loop on unmatched SYRTHES instances */

  for (i = 0; i < _syr_coupling_builder_size; i++) {

    _cs_syr_coupling_builder_t *scb = _syr_coupling_builder + i;

    _syr3_add_socket(i);

    scb->match_id = 0;
  }

  /* Cleanup (if we have not deadlocked before) */

  _remove_matched_builder_entries();
}

#endif /* defined(HAVE_SOCKETS) */

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
 fvm_lnum_t  *n_couplings
)
{
  if (_cs_glob_n_syr_cp < 0) {
    if (_syr_coupling_builder_size > 0)
      _cs_glob_n_syr_cp = _syr_coupling_builder_size;
    else
      _cs_glob_n_syr_cp =   cs_syr3_coupling_n_couplings()
                          + cs_syr4_coupling_n_couplings();
  }

  *n_couplings = _cs_glob_n_syr_cp;
}

/*----------------------------------------------------------------------------
 * Create nodal coupled mesh.
 *
 * Send vertices's coordinates and connectivity of coupled mesh for SYRTHES 3,
 * setup PLE locator for SYRTHES 4.
 *
 * Fortran Interface:
 *
 * SUBROUTINE GEOSYR
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF(geosyr, GEOSYR)
(
 void
)
{
  int coupl_id;

  _cs_glob_n_syr3_cp = cs_syr3_coupling_n_couplings();
  _cs_glob_n_syr4_cp = cs_syr4_coupling_n_couplings();

  for (coupl_id = 0; coupl_id < _cs_glob_n_syr3_cp; coupl_id++) {
    cs_syr3_coupling_t *syr_coupling = cs_syr3_coupling_by_id(coupl_id);
    cs_syr3_coupling_init_mesh(syr_coupling);
  }

  for (coupl_id = 0; coupl_id < _cs_glob_n_syr4_cp; coupl_id++) {
    cs_syr4_coupling_t *syr_coupling = cs_syr4_coupling_by_id(coupl_id);
    cs_syr4_coupling_init_mesh(syr_coupling);
  }
}

/*----------------------------------------------------------------------------
 * Check if SYRTHES 3 couplings continue or if we must stop.
 *
 * For each SYRTHES 3 coupling, A message (stop or new iteration) is
 * received. No iteration start message is sent, as this is done
 * by ITDSYR.
 *
 * Fortran Interface:
 *
 * SUBROUTINE TSTSY3 (IMSFIN)
 * *****************
 *
 * INTEGER          NTMABS      : <-> : Maximum iteration number
 * INTEGER          NTCABS      : --> : Current iteration numbern
 *----------------------------------------------------------------------------*/

void CS_PROCF(tstsy3, TSTSY3)
(
 cs_int_t *ntmabs,
 cs_int_t *ntcabs
)
{
  int nt_max_abs = *ntmabs;

  if ((*ntcabs < *ntmabs) && cs_syr3_coupling_n_couplings() > 0)
    cs_syr3_messages_test_iter(*ntcabs, &nt_max_abs);

  *ntmabs = nt_max_abs;
}

/*----------------------------------------------------------------------------
 * Synchronize new time step message for SYRTHES 3 couplings.
 *
 * For SYRTHES 3, it is necessary to distinguish the last iteration from
 * other iterations (to allow for SYRTHES 3 to determine in advance that it
 * will need to output postprocessing/restart data), so using this separate
 * function allows it to be placed after MODPAR in the main time loop,
 * in case NTMABS is changed by that function.
 *
 * Fortran Interface:
 *
 * SUBROUTINE ITDSY3 (NTCABS, NTMABS)
 * *****************
 *
 * INTEGER          NTCABS      : --> : Current iteration number
 * INTEGER          NTMABS      : --> : Maximum iteration number
 *----------------------------------------------------------------------------*/

void CS_PROCF(itdsy3, ITDSY3)
(
 cs_int_t   *ntcabs,
 cs_int_t   *ntmabs
)
{
  if (cs_syr3_coupling_n_couplings() > 0)
    cs_syr3_messages_new_time_step(*ntcabs, *ntmabs);
}

/*----------------------------------------------------------------------------
 * Get number of boundary faces coupled with SYRTHES.
 *
 * Fortran Interface:
 *
 * SUBROUTINE NBFSYR
 * *****************
 *
 * INTEGER          coupl_num       : --> : coupling number
 * INTEGER          n_coupl_faces   : <-- : number of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF(nbfsyr, NBFSYR)
(
 const cs_int_t  *coupl_num,
       cs_int_t  *n_coupl_faces
)
{
  int n_couplings = _cs_glob_n_syr3_cp + _cs_glob_n_syr4_cp;

  if (*coupl_num < 1 || *coupl_num > n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *coupl_num, n_couplings);

  else {
    if (*coupl_num <= _cs_glob_n_syr3_cp) {
      cs_syr3_coupling_t *syr_coupling
        = cs_syr3_coupling_by_id(*coupl_num - 1);
      *n_coupl_faces = cs_syr3_coupling_get_n_faces(syr_coupling);
    }
    else {
      cs_syr4_coupling_t *syr_coupling
        = cs_syr4_coupling_by_id(*coupl_num - _cs_glob_n_syr3_cp - 1);
      *n_coupl_faces = cs_syr4_coupling_get_n_faces(syr_coupling);
    }
  }
}


/*----------------------------------------------------------------------------
 * Get local numbering of coupled faces
 *
 * Fortran interface:
 *
 * SUBROUTINE LFASYR
 * *****************
 *
 * INTEGER      coupl_num       : --> : coupling number
 * INTEGER      coupl_face_list : <-- : list of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF(lfasyr, LFASYR)
(
 const cs_int_t    *coupl_num,
       fvm_lnum_t  *coupl_face_list
)
{
  int n_couplings = _cs_glob_n_syr3_cp + _cs_glob_n_syr4_cp;

  if (*coupl_num < 1 || *coupl_num > n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *coupl_num, n_couplings);
  else {
    if (*coupl_num <= _cs_glob_n_syr3_cp) {
      cs_syr3_coupling_t *syr_coupling
        = cs_syr3_coupling_by_id(*coupl_num - 1);
      cs_syr3_coupling_get_face_list(syr_coupling, coupl_face_list);
    }
    else {
      cs_syr4_coupling_t *syr_coupling
        = cs_syr4_coupling_by_id(*coupl_num - _cs_glob_n_syr3_cp - 1);
      cs_syr4_coupling_get_face_list(syr_coupling, coupl_face_list);
    }
  }
}

/*----------------------------------------------------------------------------
 * User function wrapper for definition of SYRTHES couplings
 *
 * Fortran Interface:
 *
 * SUBROUTINE USSYRC
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (ussyrc, USSYRC)
(
 void
)
{
  cs_user_syrthes_coupling();
}

/*----------------------------------------------------------------------------
 * Receive coupling variables from SYRTHES
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYI (NUMSYR, TWALL)
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * DOUBLE PRECISION TWALL       : <-- : Wall temerature
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyi, VARSYI)
(
 cs_int_t   *numsyr,
 cs_real_t  *twall
)
{
  int n_couplings = _cs_glob_n_syr3_cp + _cs_glob_n_syr4_cp;

  if (*numsyr < 1 || *numsyr > n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *numsyr, n_couplings);
  else {
    if (*numsyr <= _cs_glob_n_syr3_cp) {
      cs_syr3_messages_recv_twall(*numsyr, twall);
    }
    else {
      cs_syr4_coupling_t *syr_coupling
        = cs_syr4_coupling_by_id(*numsyr - _cs_glob_n_syr3_cp -1);
      cs_syr4_coupling_recv_twall(syr_coupling, twall);
    }
  }
}

/*----------------------------------------------------------------------------
 * Send coupling variables to SYRTHES
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYO (NUMSYR, TFLUID, HWALL)
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * DOUBLE PRECISION TFLUID      : --> : Fluid temperature
 * DOUBLE PRECISION HWALL       : --> : Exchange coefficient
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyo, VARSYO)
(
 cs_int_t   *numsyr,
 cs_real_t  *tfluid,
 cs_real_t  *hwall
)
{
  int n_couplings = _cs_glob_n_syr3_cp + _cs_glob_n_syr4_cp;

  if (*numsyr < 1 || *numsyr > n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("SYRTHES coupling number %d impossible; "
                "there are %d couplings"),
              *numsyr, n_couplings);
  else {
    if (*numsyr <= _cs_glob_n_syr3_cp)
      cs_syr3_messages_send_tf_hwall(*numsyr, tfluid, hwall);
    else {
      cs_syr4_coupling_t *syr_coupling
        = cs_syr4_coupling_by_id(*numsyr - _cs_glob_n_syr3_cp - 1);
      cs_syr4_coupling_send_tf_hwall(syr_coupling, tfluid, hwall);
    }
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define new SYRTHES coupling.
 *
 * In the case of a single Code_Saturne and single SYRTHES instance, the
 * syrthes_name argument is ignored.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances based on the syrthes_name argument.
 *
 * arguments:
 *   syrthes_name      <-- name of SYRTHES instance
 *   boundary_criteria <-- boundary face selection criteria, or NULL
 *   volume_criteria   <-- volume cell selection criteria, or NULL
 *   projection_axis   <-- 'x', 'y', or 'y' for 2D projection axis (case
 *                         independent), or ' ' for standard 3D coupling
 *   verbosity         <-- verbosity level
 *   visualization     <-- visualization output level (0 or 1)
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_define(const char  *syrthes_name,
                       const char  *boundary_criteria,
                       const char  *volume_criteria,
                       char         projection_axis,
                       int          verbosity,
                       int          visualization)
{
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

  scb->verbosity = verbosity;
  scb->visualization = visualization;

  _syr_coupling_builder_size += 1;
}

/*----------------------------------------------------------------------------
 * Initialize SYRTHES couplings.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *
 * parameters:
 *   port_num <-- port number for rank 0 to enable sockets,
 *                < 0 to disable sockets
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_init(int  port_num)
{
  /* First try using MPI */

#if defined(HAVE_MPI)

  if (_syr_coupling_builder_size > 0)
    _init_all_mpi_syr();

#endif

  /* If not all SYRTHES instances have been found, try using sockets */

#if defined(HAVE_SOCKET)

  if (_syr_coupling_builder_size > 0 && port_num > -1)
    _init_all_socket_syr(port_num);

#endif

  if (_syr_coupling_builder_size > 0) {

    bft_printf("Unmatched SYRTHES couplings:\n"
               "----------------------------\n\n");

    _print_all_unmatched_syr();

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
  cs_syr3_coupling_all_destroy();
  cs_syr4_coupling_all_destroy();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
