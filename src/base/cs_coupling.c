/*============================================================================
 * \file Common functionnality for various coupling types.
 *============================================================================*/

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>
#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_nodal_extract.h"
#include "fvm_point_location.h"

#include "cs_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

#if defined(PLE_HAVE_MPI)

static ple_coupling_mpi_set_t *_cs_glob_coupling_mpi_app_world = NULL;

#endif

/* Syncronization flag used for external couplings */

static int _cs_coupling_sync_flag = 0; /* Synchronized */


/* Optional time step multiplier for external couplings;
   the time step for the current instance times this multiplier
   is the apparent time step for coupled codes */

static double _cs_coupling_ts_multiplier = 1.0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Synchronize with applications in the same PLE coupling group.
 *
 * This function should be called before starting a new time step. The
 * current time step id is that of the last finished time step, or 0 at
 * initialization.
 *
 * Fortran Interface:
 *
 * subroutine cplsyn (ntcmabs, ntcabs, dtref)
 * *****************
 *
 * integer          ntmabs      : <-> : maximum iteration number
 * integer          ntcabs      : <-- : current iteration number
 * double precision dtref       : <-> : reference time step value
 *----------------------------------------------------------------------------*/

void CS_PROCF(cplsyn, CPLSYN)
(
 cs_int_t         *ntmabs,
 const cs_int_t   *ntcabs,
 cs_real_t        *dtref
)
{
  int  current_ts_id = *ntcabs;
  int  max_ts_id = *ntmabs;
  double  ts = *dtref;

  cs_coupling_sync_apps(0,
                        current_ts_id,
                        &max_ts_id,
                        &ts);

  *ntmabs = max_ts_id;
  *dtref = ts;
}

/*----------------------------------------------------------------------------
 * Indicate if there are synchronized applications in the same
 * PLE coupling group.
 *
 * Fortran Interface:
 *
 * subroutine cplact (isync)
 * *****************
 *
 * integer          isync       : <-- : 1 if synchronized, 0 otherwise
 *----------------------------------------------------------------------------*/

void CS_PROCF(cplact, CPLACT)
(
 cs_int_t         *isync
)
{
  if (cs_coupling_is_sync_active())
    *isync = 1;
  else
    *isync = 0;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Discover other applications in the same MPI root communicator.
 *
 * \param[in]  app_name  name of this instance of Code_Saturne
 */
/*----------------------------------------------------------------------------*/

void
cs_coupling_discover_mpi_apps(const char  *app_name,
                              const char  *forced_app_type)
{
  int mpi_flag;
  int world_size;

  MPI_Initialized(&mpi_flag);

  if (!mpi_flag)
    return;

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if (cs_glob_n_ranks < world_size) {

    int i, n_apps, app_id;

    /* App_type contains a string such as
       "Code_Saturne 2.1.0" or "NEPTUNE_CFD 1.2.1" */

    const char app_type[] = CS_APP_NAME " " CS_APP_VERSION;

    const char *sync_name[2] = {N_("point-to-point or not synchronized"),
                                N_("group synchronized")};

    const char local_add[] = N_(" (this instance)");
    const char nolocal_add[] = "";

    if (cs_glob_rank_id < 1) {
      bft_printf(_("\n"
                   "Applications accessible through MPI:\n"
                   "------------------------------------\n\n"));
      bft_printf_flush();
    }

    if (forced_app_type == NULL)
      _cs_glob_coupling_mpi_app_world
        = ple_coupling_mpi_set_create(_cs_coupling_sync_flag,
                                      app_type,
                                      app_name,
                                      MPI_COMM_WORLD,
                                      cs_glob_mpi_comm);
    else
      _cs_glob_coupling_mpi_app_world
        = ple_coupling_mpi_set_create(_cs_coupling_sync_flag,
                                      forced_app_type,
                                      app_name,
                                      MPI_COMM_WORLD,
                                      cs_glob_mpi_comm);

    n_apps = ple_coupling_mpi_set_n_apps(_cs_glob_coupling_mpi_app_world);
    app_id = ple_coupling_mpi_set_get_app_id(_cs_glob_coupling_mpi_app_world);

    if (cs_glob_rank_id < 1) {

      for (i = 0; i < n_apps; i++) {
        const char *is_local = nolocal_add;
        ple_coupling_mpi_set_info_t
          ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world,
                                             i);
        int sync_type = (ai.status & PLE_COUPLING_NO_SYNC) ? 0 : 1;
        if (i == app_id)
          is_local = _(local_add);
        bft_printf(_("  %d; type:      \"%s\"%s\n"
                     "     case name: \"%s\"\n"
                     "     lead rank: %d; n_ranks: %d\n"
                     "     (%s"),
                   i+1, ai.app_type, is_local,
                   ai.app_name, ai.root_rank, ai.n_ranks,
                   _(sync_name[sync_type]));
        if (ai.status & PLE_COUPLING_TS_MIN)
          bft_printf(_(", time step min."));
        if (ai.status & PLE_COUPLING_TS_LEADER)
          bft_printf(_(", time step leader"));
        if (ai.status & PLE_COUPLING_UNSTEADY)
          bft_printf(_(", unsteady"));
        if (ai.status & PLE_COUPLING_STEADY)
          bft_printf(_(", steady"));

        bft_printf(_(")\n\n"));
      }

      bft_printf_flush();
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize MPI coupling helper structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_coupling_finalize(void)
{
  if (_cs_glob_coupling_mpi_app_world != NULL)
    ple_coupling_mpi_set_destroy(&_cs_glob_coupling_mpi_app_world);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return info on other applications in the same MPI root communicator.
 *
 * \return  info on other applications structure
 */
/*----------------------------------------------------------------------------*/

const ple_coupling_mpi_set_t *
cs_coupling_get_mpi_apps(void)
{
  return _cs_glob_coupling_mpi_app_world;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the optional synchronization flag for external couplings.
 *
 * See \ref cs_coupling_set_sync_flag for details.
 *
 * \return  synchronization flag to apply to couplings
 */
/*----------------------------------------------------------------------------*/

int
cs_coupling_get_sync_flag(void)
{
  return _cs_coupling_sync_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define an optional synchronization flag for external couplings.
 *
 * This flag is used by all couplings based on the PLE (Parallel Location
 * and Exchange) group synchronization mechanism, which include couplings
 * with SYRTHES 4, Code_Saturne, and NEPTUNE_CFD.
 *
 * It is defined by a mask, so for example flags f1, f2, and f3 may be
 * combined using the "f1 | f2 | f2" syntax.
 *
 * Note also that for Code_Saturne, in the case of a variable time step,
 * the reference time step is synchronized at the beginning of each
 * iteration, but the actual time step is recomputed later.
 *
 * Possible flags are:
 *   PLE_COUPLING_TS_MIN        Use smallest time step
 *   PLE_COUPLING_TS_LEADER     Prescribe time step for the group
 *                              (only one member may set this flag)
 *   PLE_COUPLING_UNSTEADY      Inform others that this instance is
 *                              using an unsteady solution approach
 *   PLE_COUPLING_STEADY        Inform others that this instance is
 *                              using a teady solution approach
 *   PLE_COUPLING_USER_1        User definable flag
 *   PLE_COUPLING_USER_2        User definable flag
 *   PLE_COUPLING_USER_3        User definable flag
 *   PLE_COUPLING_USER_4        User definable flag
 *
 * To force stopping, PLE_COUPLING_STOP may be set. In this case,
 * the calculation will stop at the first synchronization, even if
 * this function is called again with another flag.
 *
 * \param[in]  flag  synchronization flag to apply to couplings
 */
/*----------------------------------------------------------------------------*/

void
cs_coupling_set_sync_flag(int flag)
{
  int stop_mask = _cs_coupling_sync_flag & PLE_COUPLING_STOP;

  _cs_coupling_sync_flag = flag | stop_mask;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the time step multiplier for external couplings.
 *
 * See \ref cs_coupling_get_ts_multiplier for details.
 *
 * \return  time step multiplier for external couplings
 */
/*----------------------------------------------------------------------------*/

double
cs_coupling_get_ts_multiplier(void)
{
  return _cs_coupling_ts_multiplier;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a time step multiplier for external couplings.
 *
 * The apparent time step for the current instance times (as viewed by
 * coupled codes) is equal to the true time step times this multiplier.
 *
 * If the synchronization flag contains "time step min" (PLE_COUPLING_TS_MIN),
 * the apparent time step is used to determine which code has the smallest
 * time step.
 *
 * \param[in]  m  time step multipier to aply to couplings
 */
/*----------------------------------------------------------------------------*/

void
cs_coupling_set_ts_multiplier(double m)
{
  _cs_coupling_ts_multiplier = m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize with applications in the same PLE coupling group.
 *
 * This function should be called before starting a new time step. The
 * current time step id is that of the last finished time step, or 0 at
 * initialization.
 *
 * Default synchronization flags indicating a new iteration or end of
 * calculation are set automatically, but the user may set additional flags
 * to this function if necessary.
 *
 * \param[in]       flags          optional additional synchronization flags
 * \param[in]       current_ts_id  current time step id
 * \param[in, out]  max_ts_id      maximum time step id
 * \param[in, out]  ts             suggested time step value
 */
/*----------------------------------------------------------------------------*/

void
cs_coupling_sync_apps(int      flags,
                      int      current_ts_id,
                      int     *max_ts_id,
                      double  *ts)
{
#if defined(PLE_HAVE_MPI)

  if (_cs_glob_coupling_mpi_app_world != NULL) {

    int i;

    int sync_flags = 0;
    int stop_mask = _cs_coupling_sync_flag & PLE_COUPLING_STOP;
    int leader_id = -1;
    double ts_min = -1.;

    double _ts = *ts * _cs_coupling_ts_multiplier;

    int n_apps
      = ple_coupling_mpi_set_n_apps(_cs_glob_coupling_mpi_app_world);
    int app_id
      = ple_coupling_mpi_set_get_app_id(_cs_glob_coupling_mpi_app_world);

    int reset_flags[] = {PLE_COUPLING_NEW_ITERATION,
                         PLE_COUPLING_REDO_ITERATION};

    const int *app_status = NULL;
    const double *app_ts = NULL;

    ple_coupling_mpi_set_info_t ai;

    /* Set synchronization flag */

    app_status
      = ple_coupling_mpi_set_get_status(_cs_glob_coupling_mpi_app_world);

    sync_flags = app_status[app_id];
    for (i = 0; i < 2; i++) {
      if (sync_flags & reset_flags[i])
        sync_flags -= reset_flags[i];
    }
    sync_flags = sync_flags | flags | stop_mask;

    if (current_ts_id >= *max_ts_id)
      sync_flags = sync_flags | PLE_COUPLING_STOP;
    else {
      sync_flags = sync_flags | PLE_COUPLING_NEW_ITERATION;
      if (current_ts_id == *max_ts_id - 1)
        sync_flags = sync_flags | PLE_COUPLING_LAST;
    }

    if (flags & PLE_COUPLING_REDO_ITERATION) {
      if (sync_flags & PLE_COUPLING_NEW_ITERATION)
        sync_flags -= PLE_COUPLING_NEW_ITERATION;
      if (sync_flags & PLE_COUPLING_STOP)
        sync_flags -= PLE_COUPLING_STOP;
    }

    /* Synchronize applications */

    ple_coupling_mpi_set_synchronize(_cs_glob_coupling_mpi_app_world,
                                     sync_flags,
                                     _ts);

    app_status
      = ple_coupling_mpi_set_get_status(_cs_glob_coupling_mpi_app_world);
    app_ts
      = ple_coupling_mpi_set_get_timestep(_cs_glob_coupling_mpi_app_world);

    /* Check if we should use the smallest time step */

    if (app_status[app_id] & PLE_COUPLING_TS_MIN)
      ts_min = _ts;

    /* Loop on applications */

    for (i = 0; i < n_apps; i++) {

      if (app_status[i] & PLE_COUPLING_NO_SYNC)
        continue;

      /* Handle leader or minimum time step update */

      if (app_status[i] & PLE_COUPLING_TS_LEADER) {
        if (leader_id > -1) {
          ple_coupling_mpi_set_info_t ai_prev
            = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world,
                                            ts_min);
          ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
          bft_error
            (__FILE__, __LINE__, 0,
             _("\nApplication \"%s\" (%s) tried to set the group time step, but\n"
               "application \"%s\" (%s) has already done so."),
             ai.app_name, ai.app_type, ai_prev.app_name, ai_prev.app_type);
        }
        else {
          leader_id = i;
          *ts = app_ts[i] / _cs_coupling_ts_multiplier;
        }
      }
      else if (app_status[i] & PLE_COUPLING_TS_MIN) {
        if (ts_min > 0)
          ts_min = CS_MIN(ts_min, app_ts[i]);
      }

      /* Handle time stepping behavior */

      if (app_status[i] & PLE_COUPLING_STOP) {
        if (*max_ts_id > current_ts_id) {
          ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
          bft_printf
            (_("\nApplication \"%s\" (%s) requested calculation stop.\n"),
             ai.app_name, ai.app_type);
          *max_ts_id = current_ts_id;
        }
      }
      else if (app_status[i] & PLE_COUPLING_REDO_ITERATION) {
        ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
        bft_error
          (__FILE__, __LINE__, 0,
           _("\nApplication \"%s\" (%s) requested restarting iteration,\n"
             "but this is not currently handled."),
           ai.app_name, ai.app_type);
      }
      else if (! (app_status[i] & PLE_COUPLING_NEW_ITERATION)) {
        ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
        bft_error
          (__FILE__, __LINE__, 0,
           _("\nApplication \"%s\" (%s) synchronized with status flag %d,\n"
             "which does not specify a known behavior."),
           ai.app_name, ai.app_type, app_status[i]);
      }

      if (app_status[i] & PLE_COUPLING_LAST) {
        if (*max_ts_id > current_ts_id + 1) {
          ai = ple_coupling_mpi_set_get_info(_cs_glob_coupling_mpi_app_world, i);
          bft_printf
            (_("\nApplication \"%s\" (%s) requested last iteration.\n"),
             ai.app_name, ai.app_type);
          *max_ts_id = current_ts_id + 1;
        }
      }

    } /* end of loop on applications */

    if (ts_min > 0)
      *ts = ts_min / _cs_coupling_ts_multiplier;
  }

#else

  return;

#endif /* PLE_HAVE_MPI */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate is synchronization with applications in the same
 *        PLE group is active.
 *
 * \return  true if synchronization is required, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_coupling_is_sync_active(void)
{
  bool retval = false;

#if defined(PLE_HAVE_MPI)

  if (_cs_glob_coupling_mpi_app_world != NULL) {

    int i;
    int n_apps
      = ple_coupling_mpi_set_n_apps(_cs_glob_coupling_mpi_app_world);
    int app_id
      = ple_coupling_mpi_set_get_app_id(_cs_glob_coupling_mpi_app_world);

    const int *app_status = NULL;

    /* Set synchronization flag */

    app_status
      = ple_coupling_mpi_set_get_status(_cs_glob_coupling_mpi_app_world);

    /* Synchronize applications */

    app_status
      = ple_coupling_mpi_set_get_status(_cs_glob_coupling_mpi_app_world);

    /* Check if we should use the smallest time step */

    if (app_status[app_id] & PLE_COUPLING_NO_SYNC)
      return retval;

    /* Loop on applications */

    for (i = 0; i < n_apps; i++) {
      if (!(app_status[i] & PLE_COUPLING_NO_SYNC))
        retval = true;
    }
  }

#endif /* PLE_HAVE_MPI */

  return retval;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute extents of a mesh representation.
 *
 * \param[in]       mesh           pointer to mesh representation structure
 * \param[in]       n_max_extents  maximum number of sub-extents (such as
 *                                 element extents) to compute, or -1 to query
 * \param[in]       tolerance      addition to local extents of each element:
 *                                 extent = base_extent * (1 + tolerance)
 * \param[in, out]  extents        extents associated with mesh:
 *                                 x_min, y_min, ..., x_max, y_max, ...
 *                                 (size: 2*dim)
 *
 * \return  the number of extents computed
 */
/*----------------------------------------------------------------------------*/

ple_lnum_t
cs_coupling_mesh_extents(const void  *mesh,
                         ple_lnum_t   n_max_extents,
                         double       tolerance,
                         double       extents[])
{
  const fvm_nodal_t  *m = mesh;
  ple_lnum_t retval = 0;

  if (m == NULL)
    return 0;

  /* In query mode, return maximum extents available
     (currently limited to 1) */

  if (n_max_extents < 0)
    retval = 1;

  /* If n_max_extents > 0 return global mesh extents */

  else if (n_max_extents > 0) {
    fvm_nodal_extents(m, tolerance, extents);
    retval = 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find elements in a given mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * Location is relative to the id of a given element + 1 in
 * concatenated sections of same element dimension.
 *
 * \param[in]       mesh                pointer to mesh representation structure
 * \param[in]       tolerance_base      associated base tolerance (for bounding
 *                                      box check only, not for location test)
 * \param[in]       tolerance_fraction  associated fraction of element bounding
 *                                      boxes added to tolerance
 * \param[in]       n_points            number of points to locate
 * \param[in]       point_coords        point coordinates
 * \param[in]       point_tag           optional point tag
 * \param[in, out]  location            number of element containing or closest
 *                                      to each point (size: n_points)
 * \param[in, out]  distance            distance from point to element indicated
 *                                      by location[]: < 0 if unlocated,
 *                                      0 - 1 if inside, and > 1 if outside a
 *                                      volume element, or absolute distance to
 *                                      a surface element (size: n_points)
 */
/*----------------------------------------------------------------------------*/

void
cs_coupling_point_in_mesh(const void         *mesh,
                          float               tolerance_base,
                          float               tolerance_fraction,
                          ple_lnum_t          n_points,
                          const ple_coord_t   point_coords[],
                          const int           point_tag[],
                          ple_lnum_t          location[],
                          float               distance[])
{
  fvm_point_location_nodal((const fvm_nodal_t *)mesh,
                           tolerance_base,
                           tolerance_fraction,
                           0, /* Do not locate on parents */
                           n_points,
                           point_tag,
                           point_coords,
                           location,
                           distance);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find elements in a given mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * Location is relative to parent element numbers.
 *
 * \param[in]       mesh                pointer to mesh representation structure
 * \param[in]       tolerance_base      associated base tolerance (for bounding
 *                                      box check only, not for location test)
 * \param[in]       tolerance_fraction  associated fraction of element bounding
 *                                      boxes added to tolerance
 * \param[in]       n_points            number of points to locate
 * \param[in]       point_coords        point coordinates
 * \param[in]       point_tag           optional point tag
 * \param[in, out]  location            number of element containing or closest
 *                                      to each point (size: n_points)
 * \param[in, out]  distance            distance from point to element indicated
 *                                      by location[]: < 0 if unlocated,
 *                                      0 - 1 if inside, and > 1 if outside a
 *                                      volume element, or absolute distance to
 *                                      a surface element (size: n_points)
 */
/*----------------------------------------------------------------------------*/

void
cs_coupling_point_in_mesh_p(const void         *mesh,
                            float               tolerance_base,
                            float               tolerance_fraction,
                            ple_lnum_t          n_points,
                            const ple_coord_t   point_coords[],
                            const int           point_tag[],
                            ple_lnum_t          location[],
                            float               distance[])
{
  fvm_point_location_nodal((const fvm_nodal_t *)mesh,
                           tolerance_base,
                           tolerance_fraction,
                           1, /* Locate on parents */
                           n_points,
                           point_tag,
                           point_coords,
                           location,
                           distance);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
