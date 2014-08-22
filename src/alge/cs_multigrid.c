/*============================================================================
 * Multigrid solver.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_grid.h"
#include "cs_halo.h"
#include "cs_log.h"
#include "cs_matrix.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_post.h"
#include "cs_sles.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_multigrid.h"

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

/* Minimum size for OpenMP loops (needs benchmarking to adjust) */
#define CS_THR_MIN 128

/* SIMD unit size to ensure SIMD alignement (2 to 4 required on most
 * current architectures, so 16 should be enough on most architectures
 * through at least 2012) */

#define CS_SIMD_SIZE(s) (((s-1)/16+1)*16)

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Local Structure Definitions
 *----------------------------------------------------------------------------*/

/* Basic per linear system options and logging */
/*---------------------------------------------*/

typedef struct _cs_multigrid_info_t {

  char                *name;                /* System name */
  cs_sles_type_t       type[3];             /* Descent/ascent smoother
                                               and solver type */

  unsigned             n_calls[2];          /* Number of times grids built
                                               (0) or solved (1) */

  unsigned long long   n_levels_tot;        /* Total accumulated number of
                                               grid levels built */
  unsigned             n_levels[3];         /* Number of grid levels:
                                               [last, min, max] */

  unsigned             n_cycles[3];         /* Number of cycles for
                                               system resolution:
                                               [min, max, total] */

  cs_timer_counter_t   t_tot[2];            /* Total time used:
                                               [build, solve] */

} cs_multigrid_info_t;

/* Per level info and logging */
/*----------------------------*/

typedef struct _cs_multigrid_level_info_t {

  unsigned long long   n_ranks[4];          /* Number of ranks for this level:
                                               [last, min, max, total] */
  unsigned long long   n_g_cells[4];        /* Global number of cells
                                               (last, min, max, total) */
  unsigned long long   n_elts[3][4];        /* Mean number of cells,
                                               cells + ghosts, and faces
                                               across ranks (last, min, max,
                                               total) */
  double               unbalance[3][4];     /* Unbalance for cells, cells
                                               + ghosts, and faces
                                               (last, min, max, total) */

  unsigned long long   n_it_solve[4];       /* Number of iterations for
                                               solving [last, min, max, total] */
  unsigned long long   n_it_ds_smoothe[4];  /* Number of iterations for
                                               descent smoothing:
                                                 [last, min, max, total] */
  unsigned long long   n_it_as_smoothe[4];  /* Number of iterations for
                                               ascent smoothing:
                                                 [last, min, max, total] */

  unsigned             n_calls[6];          /* Total number of calls:
                                               build, solve, descent smoothe,
                                               ascent smoothe, restrict from
                                               finer, prolong to finer */

  cs_timer_counter_t   t_tot[6];            /* Total timers count:
                                               [build, solve, descent smoothe,
                                               ascent smoothe, restrict from
                                               finer, prolong to finer] */

} cs_multigrid_level_info_t;

/* Grid hierarchy */
/*----------------*/

typedef struct _cs_multigrid_t {

  int         n_levels;        /* Current number of grid levels */
  int         n_levels_max;    /* Maximum number of grid levels */

  int         n_levels_post;   /* Current number of postprocessed levels */
  int         post_cell_max;   /* If > 0, activates postprocessing of
                                  coarsening, projecting coarse cell
                                  numbers (modulo post_cell_max)
                                  on the base grid */

  cs_grid_t **grid_hierarchy;  /* Array of grid pointers */

  int       **post_cell_num;   /* If post_cell_max > 0, array of
                                  (n_levels - 1) arrays of projected
                                  coarse cell numbers on the base grid */

  int       **post_cell_rank;  /* If post_cell_max > 0 and grid merging
                                  is active, array of (n_levels - 1) arrays
                                  of projected coarse cell ranks on the
                                  base grid */

  cs_multigrid_level_info_t  *lv_info;  /* Info for each level */
  cs_multigrid_info_t         info;     /* Base multigrid info */


} cs_multigrid_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static int cs_glob_multigrid_n_systems = 0;     /* Current number of systems */
static int cs_glob_multigrid_n_max_systems = 0; /* Max. number of sytems for
                                                   cs_glob_mgrid_systems. */

/* System info array */

static cs_multigrid_t **cs_glob_multigrid_systems = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize multigrid info structure.
 *
 * parameters:
 *   name <-- system name
 *   info <-- pointer to multigrid info structure
 *----------------------------------------------------------------------------*/

static void
_multigrid_info_init(cs_multigrid_info_t *info,
                     const char          *name)
{
  int i;

  BFT_MALLOC(info->name, strlen(name) + 1, char);

  strcpy(info->name, name);

  for (i = 0; i < 3; i++)
    info->type[i] = CS_SLES_N_TYPES;

  for (i = 0; i < 2; i++)
    info->n_calls[i] = 0;

  info->n_levels_tot = 0;

  for (i = 0; i < 3; i++) {
    info->n_levels[i] = 0;
    info->n_cycles[i] = 0;
  }

  for (i = 0; i < 2; i++)
    CS_TIMER_COUNTER_INIT(info->t_tot[i]);
}

/*----------------------------------------------------------------------------
 * Destroy multigrid info structure.
 *
 * parameters:
 *   this_info <-> pointer to linear system info structure pointer
 *----------------------------------------------------------------------------*/

static void
_multigrid_info_unset(cs_multigrid_info_t  *this_info)
{
  assert(this_info != NULL);

  BFT_FREE(this_info->name);
}

/*----------------------------------------------------------------------------
 * Initialize multigrid level info structure.
 *
 * parameters:
 *   name <-- system name
 *   info <-- pointer to multigrid info structure
 *----------------------------------------------------------------------------*/

static void
_multigrid_level_info_init(cs_multigrid_level_info_t *info)
{
  int i;

  memset(info, 0, sizeof(cs_multigrid_level_info_t));

  for (i = 0; i < 3; i++) {
    info->unbalance[i][0] = HUGE_VALF;
    info->unbalance[i][1] = 0.;
  }

  for (i = 0; i < 6; i++)
    CS_TIMER_COUNTER_INIT(info->t_tot[i]);
}

/*----------------------------------------------------------------------------
 * Output information regarding multigrid resolution.
 *
 * parameters:
 *   mg <-> pointer to multigrid structure
 *----------------------------------------------------------------------------*/

static void
_multigrid_info_dump(const cs_multigrid_t *mg)
{
  unsigned i;

  unsigned long long n_builds_denom = CS_MAX(mg->info.n_calls[0], 1);
  unsigned long long n_solves_denom = CS_MAX(mg->info.n_calls[1], 1);
  int n_lv_min = mg->info.n_levels[1];
  int n_lv_max = mg->info.n_levels[2];
  int n_lv_mean = (int)(mg->info.n_levels_tot / n_builds_denom);
  int n_cy_mean = (int)(mg->info.n_cycles[2] / n_solves_denom);

  char tmp_s[6][64] =  {"", "", "", "", "", ""};
  const char *stage_name[2] = {N_("Construction:"), N_("Resolution:")};
  const char *lv_stage_name[6] = {N_("build:"), N_("solve:"),
                                  N_("descent smoothe:"), N_("ascent smoothe:"),
                                  N_("restrict:"), N_("prolong:")};

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Summary of multigrid for \"%s\":\n\n"),
                mg->info.name);

  if (mg->info.type[0] != CS_SLES_N_TYPES) {

    const char *descent_smoother_name = cs_sles_type_name[mg->info.type[0]];
    const char *ascent_smoother_name = cs_sles_type_name[mg->info.type[1]];

    if (mg->info.type[0] == mg->info.type[1])
      cs_log_printf(CS_LOG_PERFORMANCE,
                    _("  Smoother: %s\n"),
                    _(descent_smoother_name));
    else
      cs_log_printf(CS_LOG_PERFORMANCE,
                    _("  Descent smoother:     %s\n"
                      "  Ascent smoother:      %s\n"),
                    _(descent_smoother_name), _(ascent_smoother_name));

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  Coarsest level solver:       %s\n"),
                  _(cs_sles_type_name[mg->info.type[2]]));

  }

  sprintf(tmp_s[0], "%-36s", "");
  cs_log_strpadl(tmp_s[1], _(" mean"), 12, 64);
  cs_log_strpadl(tmp_s[2], _("minimum"), 12, 64);
  cs_log_strpadl(tmp_s[3], _("maximum"), 12, 64);

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n  %s %s %s %s\n",
                tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

  cs_log_strpad(tmp_s[0], _("Number of coarse levels:"), 36, 64);
  cs_log_strpad(tmp_s[1], _("Number of cycles:"), 36, 64);

  cs_log_printf(CS_LOG_PERFORMANCE,
                "  %s %12d %12d %12d\n",
                tmp_s[0], n_lv_mean, n_lv_min, n_lv_max);
  cs_log_printf(CS_LOG_PERFORMANCE,
                "  %s %12d %12d %12d\n\n",
                tmp_s[1], n_cy_mean,
                (int)(mg->info.n_cycles[0]), (int)(mg->info.n_cycles[1]));

  cs_log_timer_array_header(CS_LOG_PERFORMANCE,
                            2,                  /* indent, */
                            "",                 /* header title */
                            true);              /* calls column */
  cs_log_timer_array(CS_LOG_PERFORMANCE,
                     2,                  /* indent, */
                     2,                  /* n_lines */
                     stage_name,
                     mg->info.n_calls,
                     mg->info.t_tot);

  sprintf(tmp_s[0], "%-36s", "");
  cs_log_strpadl(tmp_s[1], _(" mean"), 12, 64);
  cs_log_strpadl(tmp_s[2], _("minimum"), 12, 64);
  cs_log_strpadl(tmp_s[3], _("maximum"), 12, 64);

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n  %s %s %s %s\n",
                tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

  for (i = 0; i <= mg->info.n_levels[2]; i++) {

    const cs_multigrid_level_info_t *lv_info = mg->lv_info + i;
    unsigned long long n_lv_builds = lv_info->n_calls[0];

    if (n_lv_builds < 1)
      continue;

    cs_log_strpad(tmp_s[0], _("Number of cells:"), 34, 64);
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  Grid level %d:\n"
                    "    %s %12llu %12llu %12llu\n"),
                  i, tmp_s[0],
                  lv_info->n_g_cells[3] / n_lv_builds,
                  lv_info->n_g_cells[1], lv_info->n_g_cells[2]);

#if defined(HAVE_MPI)

    if (cs_glob_n_ranks == 1) {
      cs_log_strpad(tmp_s[1], _("Number of faces:"), 34, 64);
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "    %s %12llu %12llu %12llu\n",
                    tmp_s[1],
                    lv_info->n_elts[2][3] / n_lv_builds,
                    lv_info->n_elts[2][1], lv_info->n_elts[2][2]);
    }

#endif
    if (cs_glob_n_ranks > 1) {
      cs_log_strpad(tmp_s[0], _("Number of active ranks:"), 34, 64);
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "    %s %12llu %12llu %12llu\n",
                    tmp_s[0],
                    lv_info->n_ranks[3] / n_lv_builds,
                    lv_info->n_ranks[1], lv_info->n_ranks[2]);
      cs_log_strpad(tmp_s[0], _("Mean local cells:"), 34, 64);
      cs_log_strpad(tmp_s[1], _("Mean local cells + ghosts:"), 34, 64);
      cs_log_strpad(tmp_s[2], _("Mean local faces:"), 34, 64);
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "    %s %12llu %12llu %12llu\n"
                    "    %s %12llu %12llu %12llu\n"
                    "    %s %12llu %12llu %12llu\n",
                    tmp_s[0],
                    lv_info->n_elts[0][3] / n_lv_builds,
                    lv_info->n_elts[0][1], lv_info->n_elts[0][2],
                    tmp_s[1],
                    lv_info->n_elts[1][3] / n_lv_builds,
                    lv_info->n_elts[1][1], lv_info->n_elts[1][2],
                    tmp_s[2],
                    lv_info->n_elts[2][3] / n_lv_builds,
                    lv_info->n_elts[2][1], lv_info->n_elts[2][2]);
      cs_log_strpad(tmp_s[0], _("Cells unbalance:"), 34, 64);
      cs_log_strpad(tmp_s[1], _("Cells + ghosts unbalance:"), 34, 64);
      cs_log_strpad(tmp_s[2], _("Faces unbalance"), 34, 64);
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "    %-34s %12.3f %12.3f %12.3f\n"
                    "    %-34s %12.3f %12.3f %12.3f\n"
                    "    %-34s %12.3f %12.3f %12.3f\n",
                    tmp_s[0],
                    lv_info->unbalance[0][3] / n_lv_builds,
                    lv_info->unbalance[0][1], lv_info->unbalance[0][2],
                    tmp_s[1],
                    lv_info->unbalance[1][3] / n_lv_builds,
                    lv_info->unbalance[1][1], lv_info->unbalance[1][2],
                    tmp_s[2],
                    lv_info->unbalance[2][3] / n_lv_builds,
                    lv_info->unbalance[2][1], lv_info->unbalance[2][2]);
    }

    if (lv_info->n_calls[1] > 0) {
      cs_log_strpad(tmp_s[0], _("Iterations for solving:"), 34, 64);
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "    %s %12llu %12llu %12llu\n",
                    tmp_s[0],
                    lv_info->n_it_solve[3] / lv_info->n_calls[1],
                    lv_info->n_it_solve[1], lv_info->n_it_solve[2]);
    }

    if (lv_info->n_calls[2] > 0) {
      cs_log_strpad(tmp_s[1], _("Descent smoother iterations:"), 34, 64);
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "    %s %12llu %12llu %12llu\n",
                    tmp_s[1],
                    lv_info->n_it_ds_smoothe[3] / lv_info->n_calls[2],
                    lv_info->n_it_ds_smoothe[1], lv_info->n_it_ds_smoothe[2]);
    }

    if (lv_info->n_calls[3] > 0) {
      cs_log_strpad(tmp_s[2], _("Ascent smoother iterations:"), 34, 64);
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "    %s %12llu %12llu %12llu\n",
                    tmp_s[2],
                    lv_info->n_it_as_smoothe[3] / lv_info->n_calls[3],
                    lv_info->n_it_as_smoothe[1], lv_info->n_it_as_smoothe[2]);
    }
  }

  cs_log_timer_array_header(CS_LOG_PERFORMANCE,
                            2,                  /* indent, */
                            "",                 /* header title */
                            true);              /* calls column */

  for (i = 0; i <= mg->info.n_levels[2]; i++) {

    const cs_multigrid_level_info_t *lv_info = mg->lv_info + i;

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  Grid level %d:\n"), i);

    cs_log_timer_array(CS_LOG_PERFORMANCE,
                       4,                  /* indent, */
                       6,                  /* n_lines */
                       lv_stage_name,
                       lv_info->n_calls,
                       lv_info->t_tot);

  }

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);
}

/*----------------------------------------------------------------------------
 * Initialize multigrid info structure.
 *
 * parameters:
 *   name <-- system name
 *
 * returns:
 *   pointer to newly created multigrid structure
 *----------------------------------------------------------------------------*/

static cs_multigrid_t *
_multigrid_create(const char  *name)
{
  int ii;
  cs_multigrid_t *mg;

  BFT_MALLOC(mg, 1, cs_multigrid_t);

  _multigrid_info_init(&(mg->info), name);

  mg->n_levels = 0;
  mg->n_levels_max = 10;

  mg->n_levels_post = 0;
  mg->post_cell_max = 0;

  BFT_MALLOC(mg->grid_hierarchy, mg->n_levels_max, cs_grid_t *);
  BFT_MALLOC(mg->lv_info, mg->n_levels_max, cs_multigrid_level_info_t);

  for (ii = 0; ii < mg->n_levels_max; ii++) {
    mg->grid_hierarchy[ii] = NULL;
    _multigrid_level_info_init(mg->lv_info + ii);
  }

  mg->post_cell_num = NULL;
  mg->post_cell_rank = NULL;

  return mg;
}

/*----------------------------------------------------------------------------
 * Destroy multigrid structure.
 *
 * parameters:
 *   mg <-> pointer multigrid structure pointer
 *----------------------------------------------------------------------------*/

static void
_multigrid_destroy(cs_multigrid_t  **mg)
{
  int ii;
  cs_multigrid_t  *_mg = *mg;

  assert(*mg != NULL);

  BFT_FREE(_mg->lv_info);

  _multigrid_info_unset(&(_mg->info));

  for (ii = 0; ii < _mg->n_levels_max; ii++)
    cs_grid_destroy(_mg->grid_hierarchy + ii);

  if (_mg->post_cell_max > 0) {
    for (ii = 0; ii < _mg->n_levels_max - 1; ii++)
      if (_mg->post_cell_num[ii] != NULL)
        BFT_FREE(_mg->post_cell_num[ii]);
    BFT_FREE(_mg->post_cell_num);
  }

  if (_mg->post_cell_rank != NULL) {
    for (ii = 0; ii < _mg->n_levels_max - 1; ii++)
      if (_mg->post_cell_rank[ii] != NULL)
        BFT_FREE(_mg->post_cell_rank[ii]);
    BFT_FREE(_mg->post_cell_rank);
  }

  BFT_FREE(_mg->grid_hierarchy);

  BFT_FREE(*mg);
}

/*----------------------------------------------------------------------------
 * Add grid to multigrid structure hierarchy.
 *
 * parameters:
 *   mg   <-- multigrid structure
 *   grid <-- grid to add
 *----------------------------------------------------------------------------*/

static void
_multigrid_add_level(cs_multigrid_t  *mg,
                     cs_grid_t       *grid)
{
  int ii;

  /* Reallocate arrays if necessary */

  if (mg->n_levels == mg->n_levels_max) {

    if (mg->n_levels_max == 0)
      mg->n_levels_max = 10;
    mg->n_levels_max *= 2;

    BFT_REALLOC(mg->grid_hierarchy, mg->n_levels_max, cs_grid_t *);
    BFT_REALLOC(mg->lv_info, mg->n_levels_max, cs_multigrid_level_info_t);

    for (ii = mg->n_levels; ii < mg->n_levels_max; ii++) {
      mg->grid_hierarchy[ii] = NULL;
      _multigrid_level_info_init(mg->lv_info + ii);
    }

    if (mg->post_cell_num != NULL) {
      BFT_REALLOC(mg->post_cell_num, mg->n_levels_max, int *);
      for (ii = mg->n_levels; ii < mg->n_levels_max; ii++)
        mg->post_cell_num[ii] = NULL;
      if (mg->n_levels > 0)
        mg->post_cell_num[mg->n_levels - 1] = NULL;
    }

    if (mg->post_cell_rank != NULL) {
      BFT_REALLOC(mg->post_cell_rank, mg->n_levels_max, int *);
      for (ii = mg->n_levels; ii < mg->n_levels_max; ii++)
        mg->post_cell_rank[ii] = NULL;
      if (mg->n_levels > 0)
        mg->post_cell_rank[mg->n_levels - 1] = NULL;
    }
  }

  mg->grid_hierarchy[mg->n_levels] = grid;

  /* Update associated info */

  {
    int  n_ranks;
    cs_lnum_t  n_cells, n_cells_with_ghosts, n_faces;
    cs_gnum_t  n_g_cells;
    cs_multigrid_level_info_t  *lv_info = mg->lv_info + mg->n_levels;

    cs_grid_get_info(grid,
                     NULL,
                     NULL,
                     NULL,
                     NULL,
                     &n_ranks,
                     &n_cells,
                     &n_cells_with_ghosts,
                     &n_faces,
                     &n_g_cells);

    mg->info.n_levels[0] = mg->n_levels + 1;

    lv_info->n_ranks[0] = n_ranks;
    if (lv_info->n_ranks[1] > (unsigned)n_ranks || lv_info->n_ranks[1] == 0)
      lv_info->n_ranks[1] = n_ranks;
    else if (lv_info->n_ranks[2] < (unsigned)n_ranks)
      lv_info->n_ranks[2] = n_ranks;
    lv_info->n_ranks[3] += n_ranks;

    lv_info->n_g_cells[0] = n_g_cells;
    if (lv_info->n_g_cells[1] > n_g_cells || lv_info->n_calls[0] == 0)
      lv_info->n_g_cells[1] = n_g_cells;
    else if (lv_info->n_g_cells[2] < n_g_cells)
      lv_info->n_g_cells[2] = n_g_cells;
    lv_info->n_g_cells[3] += n_g_cells;

    lv_info->n_elts[0][0] = n_cells;
    lv_info->n_elts[1][0] = n_cells_with_ghosts;
    lv_info->n_elts[2][0] = n_faces;

    for (ii = 0; ii < 3; ii++) {
      if (   lv_info->n_elts[ii][1] > lv_info->n_elts[ii][0]
          || lv_info->n_calls[0] == 0)
        lv_info->n_elts[ii][1] = lv_info->n_elts[ii][0];
      else if (lv_info->n_elts[ii][2] < lv_info->n_elts[ii][0])
        lv_info->n_elts[ii][2] = lv_info->n_elts[ii][0];
      lv_info->n_elts[ii][3] += lv_info->n_elts[ii][0];
    }

#if defined(HAVE_MPI)

    if (cs_glob_n_ranks > 1) {
      cs_gnum_t tot_sizes[3], max_sizes[3];
      cs_gnum_t loc_sizes[3] = {n_cells, n_cells_with_ghosts, n_faces};
      MPI_Allreduce(loc_sizes, tot_sizes, 3, CS_MPI_GNUM, MPI_SUM,
                    cs_glob_mpi_comm);
      MPI_Allreduce(loc_sizes, max_sizes, 3, CS_MPI_GNUM, MPI_MAX,
                    cs_glob_mpi_comm);
      for (ii = 0; ii < 3; ii++) {
        lv_info->unbalance[ii][0] = (  max_sizes[ii]
                                     / (tot_sizes[ii]*1.0/n_ranks)) - 1.0;
        if (   lv_info->unbalance[ii][1] > lv_info->unbalance[ii][0]
            || lv_info->n_calls[0] == 0)
          lv_info->unbalance[ii][1] = lv_info->unbalance[ii][0];
        else if (lv_info->unbalance[ii][2] < lv_info->unbalance[ii][0])
          lv_info->unbalance[ii][2] = lv_info->unbalance[ii][0];
        lv_info->unbalance[ii][3] += lv_info->unbalance[ii][0];
      }
    }

#endif /* defined(HAVE_MPI) */

    lv_info->n_calls[0] += 1;
  }

  /* Ready for next level */

  mg->n_levels += 1;
}

/*----------------------------------------------------------------------------
 * Get multigrid structure's id in list of known systems
 *
 * parameters:
 *   mg <-- multigrid structure
 *
 * returns:
 *   id (0 to n-1) of structure in list of know systems, or -1 otherwise
 *----------------------------------------------------------------------------*/

static int
_multigrid_id(const cs_multigrid_t  *mg)
{
  int id;

  for (id = 0; id < cs_glob_multigrid_n_systems; id++) {
    if (mg == cs_glob_multigrid_systems[id])
      break;
  }

  if (id >= cs_glob_multigrid_n_systems)
    id = -1;

  return id;
}

/*----------------------------------------------------------------------------
 * Return pointer to linear system info.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems.
 *
 * parameters:
 *   name <-- system name
 *----------------------------------------------------------------------------*/

static cs_multigrid_t *
_find_or_add_system(const char  *name)
{
  int ii, start_id, end_id, mid_id;
  int cmp_ret = 1;

  /* Use binary search to find system */

  start_id = 0;
  end_id = cs_glob_multigrid_n_systems - 1;
  mid_id = start_id + ((end_id -start_id) / 2);

  while (start_id <= end_id) {
    cmp_ret = strcmp((cs_glob_multigrid_systems[mid_id])->info.name, name);
    if (cmp_ret < 0)
      start_id = mid_id + 1;
    else if (cmp_ret > 0)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }

  /* If found, return */

  if (cmp_ret == 0)
    return cs_glob_multigrid_systems[mid_id];

  /* Reallocate global array if necessary */

  if (cs_glob_multigrid_n_systems >= cs_glob_multigrid_n_max_systems) {

    if (cs_glob_multigrid_n_max_systems == 0)
      cs_glob_multigrid_n_max_systems = 10;
    else
      cs_glob_multigrid_n_max_systems *= 2;
    BFT_REALLOC(cs_glob_multigrid_systems,
                cs_glob_multigrid_n_max_systems,
                cs_multigrid_t *);

    for (ii = cs_glob_multigrid_n_systems;
         ii < cs_glob_multigrid_n_max_systems;
         ii++)
      cs_glob_multigrid_systems[ii] = NULL;

  }

  /* Insert in sorted list */

  for (ii = cs_glob_multigrid_n_systems; ii > mid_id; ii--)
    cs_glob_multigrid_systems[ii] = cs_glob_multigrid_systems[ii - 1];

  cs_glob_multigrid_systems[mid_id] = _multigrid_create(name);
  cs_glob_multigrid_n_systems += 1;

  return cs_glob_multigrid_systems[mid_id];
}

/*----------------------------------------------------------------------------
 * Add postprocessing info to multigrid hierarchy
 *
 * parameters:
 *   mg           <-> multigrid structure
 *   n_base_cells <-- number of cells in base grid
 *----------------------------------------------------------------------------*/

static void
_multigrid_add_post(cs_multigrid_t  *mg,
                    cs_lnum_t        n_base_cells)
{
  int ii;

  assert(mg != NULL);

  if (mg->post_cell_max < 1)
    return;

  mg->n_levels_post = mg->n_levels - 1;

  assert(mg->n_levels_post < mg->n_levels_max);

  /* Reallocate arrays if necessary */

  if (mg->post_cell_num == NULL) {
    BFT_MALLOC(mg->post_cell_num, mg->n_levels_max, int *);
    for (ii = 0; ii < mg->n_levels_max; ii++)
      mg->post_cell_num[ii] = NULL;
  }

  if (mg->post_cell_rank == NULL && cs_grid_get_merge_stride() > 1) {
    BFT_MALLOC(mg->post_cell_rank, mg->n_levels_max, int *);
    for (ii = 0; ii < mg->n_levels_max; ii++)
      mg->post_cell_rank[ii] = NULL;
  }

  for (ii = 0; ii < mg->n_levels_post; ii++) {
    BFT_REALLOC(mg->post_cell_num[ii], n_base_cells, int);
    cs_grid_project_cell_num(mg->grid_hierarchy[ii+1],
                             n_base_cells,
                             mg->post_cell_max,
                             mg->post_cell_num[ii]);
  }

  if (mg->post_cell_rank != NULL) {
    for (ii = 0; ii < mg->n_levels_post; ii++) {
      BFT_REALLOC(mg->post_cell_rank[ii], n_base_cells, int);
      cs_grid_project_cell_rank(mg->grid_hierarchy[ii+1],
                                n_base_cells,
                                mg->post_cell_rank[ii]);
    }
  }
}

/*----------------------------------------------------------------------------
 * Post process variables associated with Multigrid hierarchy
 *
 * parameters:
 *   mgh <-- multigrid hierarchy
 *   ts  <-- time step status structure
 *----------------------------------------------------------------------------*/

static void
_cs_multigrid_post_function(void                  *mgh,
                            const cs_time_step_t  *ts)
{
  int ii;
  size_t name_len;
  char *var_name = NULL;
  cs_multigrid_t *mg = mgh;
  const char *base_name = NULL;
  const int nt_cur = (ts != NULL) ? ts->nt_cur : -1;

  /* Return if necessary structures inconsistent or have been destroyed */

  if (mg == NULL)
    return;

  if (mg->post_cell_num == NULL || cs_post_mesh_exists(-1) != true)
    return;

  /* Allocate name buffer */

  base_name = mg->info.name;
  name_len = 3 + strlen(base_name) + 1 + 3 + 1 + 4 + 1;
  BFT_MALLOC(var_name, name_len, char);

  /* Loop on grid levels */

  for (ii = 0; ii < mg->n_levels_post; ii++) {

    sprintf(var_name, "mg %s %2d %3d",
            base_name, (ii+1), nt_cur);

    cs_post_write_var(-1,
                      var_name,
                      1,
                      false,
                      true,
                      CS_POST_TYPE_int,
                      mg->post_cell_num[ii],
                      NULL,
                      NULL,
                      NULL);

    BFT_FREE(mg->post_cell_num[ii]);

    if (mg->post_cell_rank != NULL) {

      sprintf(var_name, "rk %s %2d %3d",
              base_name, (ii+1), nt_cur);

      cs_post_write_var(-1,
                        var_name,
                        1,
                        false,
                        true,
                        CS_POST_TYPE_int,
                        mg->post_cell_rank[ii],
                        NULL,
                        NULL,
                        NULL);

      BFT_FREE(mg->post_cell_rank[ii]);

    }

  }
  mg->n_levels_post = 0;

  BFT_FREE(var_name);
}

/*----------------------------------------------------------------------------
 * Compute dot product, summing result over all ranks.
 *
 * parameters:
 *   n_elts <-- local number of elements
 *   x      <-- first vector in s = x.y
 *   y      <-- second vector in s = x.y
 *
 * returns:
 *   result of s = x.y
 *----------------------------------------------------------------------------*/

inline static double
_dot_product(cs_int_t          n_elts,
             const cs_real_t  *x,
             const cs_real_t  *y)
{
  double s = cs_dot(n_elts, x, y);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    double _sum;
    MPI_Allreduce(&s, &_sum, 1, MPI_DOUBLE, MPI_SUM, cs_glob_mpi_comm);
    s = _sum;
  }

#endif /* defined(HAVE_MPI) */

  return s;
}

/*----------------------------------------------------------------------------
 * Test if convergence is attained.
 *
 * parameters:
 *   var_name      <-- variable name
 *   n_f_cells     <-- number of cells on fine mesh
 *   n_max_cycles  <-- maximum number of cycles
 *   cycle_id      <-- number of current cycle
 *
 *   verbosity     <-- verbosity level
 *   n_iters       <-- number of iterations
 *   precision     <-- precision limit
 *   r_norm        <-- residue normalization
 *   residue       <-> residue
 *   rhs           <-- right-hand side
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   cycle number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

static int
_convergence_test(const char         *var_name,
                  cs_int_t            n_f_cells,
                  int                 n_max_cycles,
                  int                 cycle_id,
                  int                 verbosity,
                  int                 n_iters,
                  double              precision,
                  double              r_norm,
                  double             *residue,
                  const cs_real_t     rhs[])
{
  const char cycle_h_fmt[]
    = N_("  ---------------------------------------------------\n"
         "    n.     | Cumulative iterations | Norm. residual\n"
         "    cycles | on fine mesh          | on fine mesh\n"
         "  ---------------------------------------------------\n");
  const char cycle_t_fmt[]
    = N_("  ---------------------------------------------------\n");
  const char cycle_cv_fmt[]
    = N_("     %4d  |               %6d  |  %12.4e\n");

  const char cycle_fmt[]
    = N_("   N. cycles: %4d; Fine mesh cumulative iter: %5d; "
         "Norm. residual %12.4e\n");

  static double initial_residue = 0.;

  /* Compute residue */

  *residue = sqrt(_dot_product(n_f_cells, rhs, rhs));

  if (cycle_id == 1)
    initial_residue = *residue;

  if (*residue < precision*r_norm) {

    if (verbosity == 2)
      bft_printf(_(cycle_fmt), cycle_id, n_iters, *residue/r_norm);
    else if (verbosity > 2) {
      bft_printf(_(cycle_h_fmt));
      bft_printf(_(cycle_cv_fmt),
                 cycle_id, n_iters, *residue/r_norm);
      bft_printf(_(cycle_t_fmt));
    }
    return 1;
  }

  else if (cycle_id >= n_max_cycles) {

    if (verbosity > 0) {
      if (verbosity == 1)
        bft_printf(_(cycle_fmt), cycle_id, n_iters, *residue/r_norm);
      else if (verbosity > 1) {
        bft_printf(_(cycle_h_fmt));
        bft_printf(_(cycle_fmt),
                   cycle_id, n_iters, *residue/r_norm);
        bft_printf(_(cycle_t_fmt));
      }
      bft_printf(_(" @@ Warning: algebraic multigrid for [%s]\n"
                   "    ********\n"
                   "    Maximum number of cycles (%d) reached.\n"),
                 var_name, n_max_cycles);

    }
    return -1;
  }

  else {

    if (verbosity > 2)
      bft_printf(_(cycle_fmt), cycle_id, n_iters, *residue/r_norm);

    if (*residue > initial_residue * 10000.0 && *residue > 100.)
      return -2;

#if (__STDC_VERSION__ >= 199901L)
    if (isnan(*residue) || isinf(*residue))
      return -2;
#endif
  }

  return 0;
}

/*----------------------------------------------------------------------------
 * Handle error output in case if divergence.
 *
 * Depending on the multigrid level at which divergence is detected,
 * output matrix hierarchy information, RHS, and variable info,
 * and finally abort on error.
 *
 * parameters:
 *   mg              <-- multigrid system
 *   level           <-- multigrid level at which divergence is detected
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   cycle_id        <-- id of currect cycle
 *   initial_residue <-- initial residue
 *   residue         <-- residue
 *   rhs             <-- right hand side
 *   vx              <-- system solution
 *   c_rhs           <-- right hand side for levels > 0
 *   c_vx            <-- system solution for levels > 0
 *----------------------------------------------------------------------------*/

static void
_abort_on_divergence(cs_multigrid_t      *mg,
                     int                  level,
                     cs_halo_rotation_t   rotation_mode,
                     int                  cycle_id,
                     double               initial_residue,
                     double               residue,
                     const cs_real_t      rhs[],
                     cs_real_t            vx[],
                     cs_real_t          **c_rhs,
                     cs_real_t          **c_vx)
{
  int mesh_id = cs_post_init_error_writer_cells();

  if (mesh_id != 0) {

    char var_name[32];

    int lv_id = 0;
    cs_real_t *var = NULL, *da = NULL;

    int i;
    int db_size[4] = {1, 1, 1, 1};
    int eb_size[4] = {1, 1, 1, 1};

    const cs_grid_t *g = mg->grid_hierarchy[0];
    const cs_lnum_t n_base_cells = cs_grid_get_n_cells(g);
    const cs_matrix_t  *_matrix = NULL;

    BFT_MALLOC(var, cs_grid_get_n_cells_ext(g), cs_real_t);
    BFT_MALLOC(da, cs_grid_get_n_cells_ext(g), cs_real_t);

    /* Output info on main level */

    _matrix = cs_grid_get_matrix(g);

    cs_sles_post_error_output_def(mg->info.name,
                                  mesh_id,
                                  rotation_mode,
                                  _matrix,
                                  rhs,
                                  vx);

    /* Output diagonal and diagonal dominance for all coarse levels */

    for (lv_id = 1; lv_id < mg->n_levels; lv_id++) {

      g = mg->grid_hierarchy[lv_id];

      cs_grid_get_info(g,
                       NULL,
                       NULL,
                       db_size,
                       eb_size,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL);


      _matrix = cs_grid_get_matrix(g);

      cs_matrix_copy_diagonal(_matrix, da);
      cs_grid_project_var(g, n_base_cells, da, var);
      sprintf(var_name, "Diag_%04d", lv_id);
      cs_sles_post_error_output_var(var_name, mesh_id, db_size[1], var);

      cs_grid_project_diag_dom(g, n_base_cells, var);
      sprintf(var_name, "Diag_Dom_%04d", lv_id);
      cs_sles_post_error_output_var(var_name, mesh_id, db_size[1], var);
    }

    /* Output info on current level if > 0 */

    if (level > 0) {

      cs_lnum_t ii;
      cs_lnum_t n_cells = 0;
      cs_lnum_t n_cells_ext = 0;

      cs_real_t *c_res = NULL;

      g = mg->grid_hierarchy[level];

      cs_grid_get_info(g,
                       NULL,
                       NULL,
                       db_size,
                       eb_size,
                       NULL,
                       &n_cells,
                       &n_cells_ext,
                       NULL,
                       NULL);

      cs_grid_project_var(g, n_base_cells, c_rhs[level], var);
      sprintf(var_name, "RHS_%04d", level);
      cs_sles_post_error_output_var(var_name, mesh_id, db_size[1], var);

      cs_grid_project_var(g, n_base_cells, c_vx[level], var);
      sprintf(var_name, "X_%04d", level);
      cs_sles_post_error_output_var(var_name, mesh_id, db_size[1], var);

      /* Compute residual */

      BFT_MALLOC(c_res, n_cells_ext*db_size[1], cs_real_t);

      _matrix = cs_grid_get_matrix(g);

      cs_matrix_vector_multiply(rotation_mode, _matrix, c_vx[level], c_res);

      for (ii = 0; ii < n_cells; ii++) {
        for (i = 0; i < db_size[0]; i++)
          c_res[ii*db_size[1] + i] = fabs( c_res[ii*db_size[1] + i]
                                         - c_rhs[level][ii*db_size[1] + i]);
      }

      cs_grid_project_var(g, n_base_cells, c_res, var);

      BFT_FREE(c_res);

      sprintf(var_name, "Residual_%04d", level);
      cs_sles_post_error_output_var(var_name, mesh_id, db_size[1], var);
    }

    cs_post_finalize();

    BFT_FREE(da);
    BFT_FREE(var);
  }

  /* Now abort */

  if (level == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("algebraic multigrid [%s]: divergence after %d cycles:\n"
                "  initial residual: %11.4e; current residual: %11.4e"),
              _(mg->info.name), cycle_id, initial_residue, residue);
  else
    bft_error(__FILE__, __LINE__, 0,
              _("algebraic multigrid [%s]: divergence after %d cycles\n"
                "  during resolution at level %d:\n"
                "  initial residual: %11.4e; current residual: %11.4e"),
              _(mg->info.name), cycle_id, level, initial_residue, residue);
}

/*----------------------------------------------------------------------------
 * Update level information iteration counts
 *
 * parameters:
 *   lv_info_it <-> logged number of iterations (last, min, max, total)
 *   n_iter     <-- current number of iterations
 *----------------------------------------------------------------------------*/

static inline void
_lv_info_update_stage_iter(unsigned long long  lv_info_it[],
                           unsigned            n_iter)
{
  lv_info_it[0] = n_iter;
  if (n_iter < lv_info_it[1] || lv_info_it[1] == 0)
    lv_info_it[1] = n_iter;
  else if (n_iter > lv_info_it[2])
    lv_info_it[2] = n_iter;
  lv_info_it[3] += n_iter;
}

/*----------------------------------------------------------------------------
 * Sparse linear system resolution using multigrid.
 *
 * parameters:
 *   mg                    <-- multigrid system
 *   descent_smoother_type <-- type of smoother for descent (PCG, Jacobi, ...)
 *   ascent_smoother_type  <-- type of smoother for ascent (PCG, Jacobi, ...)
 *   coarse_solver_type    <-- type of solver (PCG, Jacobi, ...)
 *   abort_on_divergence   <-- call errorhandler if devergence is detected
 *   poly_degree           <-- preconditioning polynomial degree (0: diagonal)
 *   rotation_mode         <-- halo update option for rotational periodicity
 *   verbosity             <-- verbosity level
 *   cycle_id              <-- id of currect cycle
 *   n_max_cycles          <-- maximum number of cycles
 *   n_max_iter            <-- maximum number of iterations per grid level
 *                             n_max_iter[level * 2]     for descent
 *                             n_max_iter[level * 2 + 1] for ascent
 *   n_equiv_iter          <-> equivalent number of iterations
 *   precision             <-- precision limit
 *   r_norm                <-- residue normalization
 *   residue               <-> residue
 *   rhs                   <-- right hand side
 *   vx                    --> system solution
 *   aux_size              <-- number of elements in aux_vectors
 *   aux_vectors           --- optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   cycle number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

static int
_multigrid_cycle(cs_multigrid_t      *mg,
                 cs_sles_type_t       descent_smoother_type,
                 cs_sles_type_t       ascent_smoother_type,
                 cs_sles_type_t       coarse_solver_type,
                 bool                 abort_on_divergence,
                 int                  poly_degree,
                 cs_halo_rotation_t   rotation_mode,
                 int                  verbosity,
                 int                  cycle_id,
                 int                  n_max_cycles,
                 const int            n_max_iter[],
                 int                 *n_equiv_iter,
                 double               precision,
                 double               r_norm,
                 double              *residue,
                 const cs_real_t     *rhs,
                 cs_real_t           *vx,
                 size_t               aux_size,
                 void                *aux_vectors)
{
  int level, coarsest_level;
  cs_lnum_t ii, jj;
  cs_timer_t t0, t1;

  int db_size[4] = {1, 1, 1, 1};
  int eb_size[4] = {1, 1, 1, 1};
  int cvg = 0, c_cvg = 0;
  int n_iter = 0;
  size_t alloc_size = 0, wr_size = 0;
  cs_real_t c_precision = precision;
  cs_real_t _residue = -1.;

  size_t _aux_size = aux_size;
  cs_lnum_t n_cells = 0, n_cells_ext = 0;
  cs_gnum_t n_g_cells = 0;
  cs_real_t r_norm_l = r_norm;

  char _var_lv_name[33];
  char *var_lv_name = _var_lv_name;

  double denom_n_g_cells_0 = 1.0;
  double _initial_residue = 0.;

  cs_real_t *_aux_vectors = aux_vectors;
  cs_real_t *wr = NULL;
  cs_real_t *_rhs_vx_val = NULL;
  cs_real_t **_rhs_vx = NULL, **_rhs = NULL, **_vx = NULL;
  cs_multigrid_info_t *mg_info = NULL;
  cs_multigrid_level_info_t  *lv_info = NULL;

  const char *var_name = NULL;
  const cs_real_t *_rhs_level = NULL;
  const cs_matrix_t  *_matrix = NULL;
  const cs_grid_t *f = NULL, *c= NULL;

  bool end_cycle = false;

  /* Initialization */

  mg_info = &(mg->info);
  var_name = mg_info->name;

  coarsest_level = mg->n_levels - 1;

  /* In theory, one should increase precision on coarsest mesh,
     but in practice, it is more efficient to have a lower precision */
  /* c_precision = precision * 0.01; */
  c_precision = precision;

  f = mg->grid_hierarchy[0];

  cs_grid_get_info(f,
                   NULL,
                   NULL,
                   db_size,
                   eb_size,
                   NULL,
                   &n_cells,
                   &n_cells_ext,
                   NULL,
                   &n_g_cells);

  denom_n_g_cells_0 = 1.0 / n_g_cells;

  if (strlen(var_name) + 5 > 32)
    BFT_MALLOC(var_lv_name, strlen(var_name) + 5 + 1, char);

  /* Allocate wr or use working area */

  for (level = 1, wr_size = n_cells_ext*db_size[1];
       level < mg->n_levels;
       level++) {
    cs_lnum_t n_cells_max
      = cs_grid_get_n_cells_max(mg->grid_hierarchy[level]);
    wr_size = CS_MAX(wr_size, (size_t)(n_cells_max*db_size[1]));
    wr_size = CS_SIMD_SIZE(wr_size);
  }

  if (aux_size >= wr_size) {
    wr = aux_vectors;
    _aux_vectors = wr + wr_size;
    _aux_size = aux_size - wr_size;
  }
  else
    BFT_MALLOC(wr, wr_size, cs_real_t);

  /* reserve memory for rhs and vx;
     for the finest level, simply point to input and output arrays */

  BFT_MALLOC(_rhs_vx, mg->n_levels*2*db_size[1], cs_real_t *);
  _rhs = _rhs_vx;
  _vx = _rhs_vx + mg->n_levels*db_size[1];

  _rhs[0] = NULL; /* Use _rhs_level when necessary to avoid const warning */
  _vx[0] = vx;

  /* Reserve memory for corrections and residues for coarse levels */

  if (mg->n_levels > 1) {

    alloc_size = 0;

    for (level = 1; level < mg->n_levels; level++)
      alloc_size
        += CS_SIMD_SIZE(cs_grid_get_n_cells_max(mg->grid_hierarchy[level]));

    BFT_MALLOC(_rhs_vx_val, alloc_size*2*db_size[1], cs_real_t);

    _rhs[1] = _rhs_vx_val;
    _vx[1] = _rhs_vx_val + alloc_size*db_size[1];

    for (level = 2; level < mg->n_levels; level++) {
      cs_lnum_t _n_cells_ext_prev
        = CS_SIMD_SIZE(cs_grid_get_n_cells_max(mg->grid_hierarchy[level-1]));
      _rhs[level] = _rhs[level - 1] + _n_cells_ext_prev*db_size[1];
      _vx[level] = _vx[level - 1] + _n_cells_ext_prev*db_size[1];
    }
  }

  /* Descent */
  /*---------*/

  if (verbosity > 2)
    bft_printf(_("  Multigrid cycle: descent\n"));

  for (level = 0; level < coarsest_level; level++) {

    lv_info = mg->lv_info + level;
    t0 = cs_timer_time();

    _rhs_level = (level == 0) ?  rhs : _rhs[level];

    sprintf(var_lv_name, "%s:%04d", var_name, level);

    c = mg->grid_hierarchy[level+1];

    /* Smoother pass */

    if (verbosity > 2)
      bft_printf(_("    level %3d: smoother\n"), level);

    _matrix = cs_grid_get_matrix(f);

    _initial_residue = _residue;

#if defined(HAVE_MPI)
    cs_sles_set_mpi_reduce_comm(cs_grid_get_comm(f));
#endif

    c_cvg = cs_sles_solve(var_lv_name,
                          descent_smoother_type,
                          false, /* Stats not updated here */
                          _matrix,
                          poly_degree,
                          rotation_mode,
                          verbosity - 2,
                          n_max_iter[level*2],
                          precision,
                          r_norm_l,
                          &n_iter,
                          &_residue,
                          _rhs_level,
                          _vx[level],
                          _aux_size,
                          _aux_vectors);

#if defined(HAVE_MPI)
    cs_sles_set_mpi_reduce_comm(cs_glob_mpi_comm);
#endif

    if (c_cvg == -2) {
      end_cycle = true;
      break;
    }

    /* Restrict residue
       TODO: get residue from cs_sles_solve(). This optimisation would
       require adding an argument and exercising caution to ensure the
       correct sign and meaning of the residue
       (regarding timing, this stage is part of the descent smoother) */

    cs_matrix_vector_multiply(rotation_mode,
                              _matrix,
                              _vx[level],
                              wr);

    if (db_size[0] == 1) {
#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (ii = 0; ii < n_cells; ii++)
        wr[ii] = _rhs_level[ii] - wr[ii];
    }
    else {
#     pragma omp parallel for private(jj) if(n_cells > CS_THR_MIN)
      for (ii = 0; ii < n_cells; ii++) {
        for (jj = 0; jj < db_size[0]; jj++)
        wr[ii*db_size[1] + jj] =   _rhs_level[ii*db_size[1] + jj]
                                 - wr[ii*db_size[1] + jj];
      }
    }

    /* Convergence test in beginning of cycle (fine mesh) */

    if (level == 0) {

      cvg = _convergence_test(var_name,
                              n_cells*db_size[1],
                              n_max_cycles,
                              cycle_id,
                              verbosity,
                              lv_info->n_it_ds_smoothe[0],
                              precision,
                              r_norm,
                              residue,
                              wr);

      /* If converged or cycle limit reached, break from descent loop */

      if (cvg != 0) {
        c_cvg = cvg;
        end_cycle = true;
        t1 = cs_timer_time();
        cs_timer_counter_add_diff(&(lv_info->t_tot[2]), &t0, &t1);
        lv_info->n_calls[2] += 1;
        _lv_info_update_stage_iter(lv_info->n_it_ds_smoothe, n_iter);

        *n_equiv_iter += n_iter * n_g_cells * denom_n_g_cells_0;
        break;
      }

    }

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[2]), &t0, &t1);
    lv_info->n_calls[2] += 1;
    _lv_info_update_stage_iter(lv_info->n_it_ds_smoothe, n_iter);

    *n_equiv_iter += n_iter * n_g_cells * denom_n_g_cells_0;

    /* Prepare for next level */

    cs_grid_restrict_cell_var(f, c, wr, _rhs[level+1]);

    cs_grid_get_info(c,
                     NULL,
                     NULL,
                     NULL,
                     NULL,
                     NULL,
                     &n_cells,
                     &n_cells_ext,
                     NULL,
                     &n_g_cells);

    f = c;

    /* Initialize correction */

    if (db_size[0] == 1) {
#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (ii = 0; ii < n_cells; ii++)
        _vx[level+1][ii] = 0.0;
    }
    else {
#     pragma omp parallel for private(jj) if(n_cells > CS_THR_MIN)
      for (ii = 0; ii < n_cells; ii++) {
        for (jj = 0; jj < db_size[0]; jj++)
          _vx[level+1][ii*db_size[1] + jj] = 0.0;
      }
    }

    t0 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[4]), &t1, &t0);
    lv_info->n_calls[4] += 1;

  } /* End of loop on levels (descent) */

  if (end_cycle == false) {

    /* Resolve coarsest level to convergence */
    /*---------------------------------------*/

    if (verbosity > 2)
      bft_printf(_("  Resolution on coarsest level\n"));

    assert(level == coarsest_level);
    assert(c == mg->grid_hierarchy[coarsest_level]);

    /* coarsest level == 0 should never happen, but we play it safe */
    _rhs_level = (level == 0) ?  rhs : _rhs[coarsest_level];

    sprintf(var_lv_name, "%s:%04d", var_name, coarsest_level);

    _matrix = cs_grid_get_matrix(c);

    _initial_residue = _residue;

    lv_info = mg->lv_info + level;
    t0 = cs_timer_time();

#if defined(HAVE_MPI)
    cs_sles_set_mpi_reduce_comm(cs_grid_get_comm(c));
#endif

    c_cvg = cs_sles_solve(var_lv_name,
                          coarse_solver_type,
                          false, /* Stats not updated here */
                          _matrix,
                          poly_degree,
                          rotation_mode,
                          verbosity - 2,
                          n_max_iter[level*2],
                          c_precision,
                          r_norm_l,
                          &n_iter,
                          &_residue,
                          _rhs_level,
                          _vx[level],
                          _aux_size,
                          _aux_vectors);

#if defined(HAVE_MPI)
    cs_sles_set_mpi_reduce_comm(cs_glob_mpi_comm);
#endif

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[1]), &t0, &t1);
    lv_info->n_calls[1] += 1;
    _lv_info_update_stage_iter(lv_info->n_it_solve, n_iter);

    *n_equiv_iter += n_iter * n_g_cells * denom_n_g_cells_0;

    if (c_cvg == -2)
      end_cycle = true;

  }

  if (end_cycle == false) {

    /* Ascent */
    /*--------*/

    if (verbosity > 2)
      bft_printf(_("  Multigrid cycle: ascent\n"));

    for (level = coarsest_level - 1; level > -1; level--) {

      cs_real_t *_f_vx = _vx[level];

      lv_info = mg->lv_info + level;

      c = mg->grid_hierarchy[level+1];
      f = mg->grid_hierarchy[level];

      cs_grid_get_info(f,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       &n_cells,
                       &n_cells_ext,
                       NULL,
                       &n_g_cells);

      /* Prolong correction */

      t0 = cs_timer_time();

      cs_grid_prolong_cell_var(c, f, _vx[level+1], wr);

      if (db_size[0] == 1) {
#       pragma omp parallel for if(n_cells > CS_THR_MIN)
        for (ii = 0; ii < n_cells; ii++)
          _f_vx[ii] += wr[ii];
      }
      else {
#       pragma omp parallel for private(jj) if(n_cells > CS_THR_MIN)
        for (ii = 0; ii < n_cells; ii++) {
          for (jj = 0; jj < db_size[0]; jj++)
            _f_vx[ii*db_size[1]+jj] += wr[ii*db_size[1]+jj];
        }
      }

      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&(lv_info->t_tot[5]), &t0, &t1);
      lv_info->n_calls[5] += 1;

      /* Smoother pass if level > 0
         (smoother not called for finest mesh, as it will be called in
         descent phase of the next cycle, before the convergence test). */

      if (level > 0) {

        if (verbosity > 2)
          bft_printf(_("    level %3d: smoother\n"), level);

        sprintf(var_lv_name, "%s:%04d", var_name, level);

        _matrix = cs_grid_get_matrix(f);

        _initial_residue = _residue;

#if defined(HAVE_MPI)
        cs_sles_set_mpi_reduce_comm(cs_grid_get_comm(f));
#endif

        c_cvg = cs_sles_solve(var_lv_name,
                              ascent_smoother_type,
                              false, /* Stats not updated here */
                              _matrix,
                              poly_degree,
                              rotation_mode,
                              verbosity - 2,
                              n_max_iter[level*2 + 1],
                              precision,
                              r_norm_l,
                              &n_iter,
                              &_residue,
                              _rhs[level],
                              _vx[level],
                              _aux_size,
                              _aux_vectors);

#if defined(HAVE_MPI)
        cs_sles_set_mpi_reduce_comm(cs_glob_mpi_comm);
#endif

        t0 = cs_timer_time();
        cs_timer_counter_add_diff(&(lv_info->t_tot[3]), &t1, &t0);
        lv_info->n_calls[3] += 1;
        _lv_info_update_stage_iter(lv_info->n_it_as_smoothe, n_iter);

        *n_equiv_iter += n_iter * n_g_cells * denom_n_g_cells_0;

        if (c_cvg == -2)
          break;
      }

    } /* End loop on levels (ascent) */

  } /* End of tests on end_cycle */

  if (c_cvg == -2) {
    cvg = -2;
    if (abort_on_divergence)
      _abort_on_divergence(mg, level,
                           rotation_mode, cycle_id,
                           _initial_residue, _residue,
                           rhs, vx, _rhs, _vx);
  }

  /* Free memory */

  if (var_lv_name != _var_lv_name)
    BFT_FREE(var_lv_name);

  if (wr != aux_vectors)
    BFT_FREE(wr);

  BFT_FREE(_rhs_vx);
  BFT_FREE(_rhs_vx_val);

  return cvg;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Build a hierarchy of meshes starting from a fine mesh, for an
 * ACM (Additive Corrective Multigrid) method.
 *----------------------------------------------------------------------------*/

void CS_PROCF(clmlga, CLMLGA)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname,     /* <-- variable name length */
 const cs_int_t   *isym,      /* <-- Symmetry indicator:
                                     1: symmetric; 2: not symmetric */
 const cs_int_t   *ibsize,    /* <-- Matrix block size */
 const cs_int_t   *iesize,    /* <-- Matrix extra diag block size */
 const cs_int_t   *nagmax,    /* <-- Agglomeration count limit */
 const cs_int_t   *ncpost,    /* <-- If > 0, postprocess coarsening, using
                                     coarse cell numbers modulo ncpost */
 const cs_int_t   *iwarnp,    /* <-- Verbosity level */
 const cs_int_t   *ngrmax,    /* <-- Maximum number of grid levels */
 const cs_int_t   *ncegrm,    /* <-- Maximum local number of cells on
                                     coarsest grid */
 const cs_real_t  *rlxp1,     /* <-- P0/P1 relaxation parameter */
 const cs_real_t  *dam,       /* <-- Matrix diagonal */
 const cs_real_t  *xam        /* <-- Matrix extra-diagonal terms */
)
{
  char *var_name;

  bool symmetric = (*isym == 1) ? true : false;
  int diag_block_size[4] = {*ibsize, *ibsize, *ibsize, (*ibsize)*(*ibsize)};
  int extra_diag_block_size[4] = {*iesize, *iesize, *iesize, (*iesize)*(*iesize)};

  var_name = cs_base_string_f_to_c_create(cname, *lname);

  cs_multigrid_build(var_name,
                     *iwarnp,
                     *ncpost,
                     *nagmax,
                     *ngrmax,
                     *ncegrm,
                     *rlxp1,
                     symmetric,
                     diag_block_size,
                     extra_diag_block_size,
                     dam,
                     xam);

  cs_base_string_f_to_c_free(&var_name);
}

/*----------------------------------------------------------------------------
 * Destroy a hierarchy of meshes starting from a fine mesh, keeping
 * the corresponding system and postprocessing information for future calls.
 *----------------------------------------------------------------------------*/

void CS_PROCF(dsmlga, DSMLGA)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname      /* <-- variable name length */
)
{
  char *var_name;

  var_name = cs_base_string_f_to_c_create(cname, *lname);

  cs_multigrid_destroy(var_name);

  cs_base_string_f_to_c_free(&var_name);
}

/*----------------------------------------------------------------------------
 * General sparse linear system resolution
 *----------------------------------------------------------------------------*/

void CS_PROCF(resmgr, RESMGR)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname,     /* <-- variable name length */
 const cs_int_t   *iresds,    /* <-- Descent smoother type:
                                     0: pcg; 1: Jacobi; 2: cg-stab,
                                     200: pcg_single reduction */
 const cs_int_t   *iresas,    /* <-- Ascent smoother type:
                                     0: pcg; 1: Jacobi; 2: cg-stab,
                                     200: pcg_single reduction */
 const cs_int_t   *ireslp,    /* <-- Coarse Resolution type:
                                     0: pcg; 1: Jacobi; 2: cg-stab,
                                     200: pcg_single reduction */
 const cs_int_t   *ipol,      /* <-- Preconditioning polynomial degree
                                     (0: diagonal, -1: none) */
 const cs_int_t   *ncymxp,    /* <-- Max number of cycles */
 const cs_int_t   *nitmds,    /* <-- Max number of iterations for descent */
 const cs_int_t   *nitmas,    /* <-- Max number of iterations for ascent */
 const cs_int_t   *nitmap,    /* <-- Max number of iterations for
                                     coarsest solution */
 const cs_int_t   *iinvpe,    /* <-- Indicator to cancel increments
                                     in rotational periodicity (2) or
                                     to exchange them as scalars (1) */
 const cs_int_t   *iwarnp,    /* <-- Verbosity level */
 cs_int_t         *ncyclf,    /* --> Number of cycles done */
 cs_int_t         *niterf,    /* --> Number of iterations done */
 const cs_real_t  *epsilp,    /* <-- Precision for iterative resolution */
 const cs_real_t  *rnorm,     /* <-- Residue normalization */
 cs_real_t        *residu,    /* --> Final non normalized residue */
 const cs_real_t  *rhs,       /* <-- System right-hand side */
 cs_real_t        *vx         /* <-> System solution */
)
{
  char *var_name;

  int res_type[3] = {*iresds, *iresas, *ireslp};

  cs_halo_rotation_t rotation_mode = CS_HALO_ROTATION_COPY;

  if (*iinvpe == 2)
    rotation_mode = CS_HALO_ROTATION_ZERO;
  else if (*iinvpe == 3)
    rotation_mode = CS_HALO_ROTATION_IGNORE;

  var_name = cs_base_string_f_to_c_create(cname, *lname);

  for (int i = 0; i < 3; i++) {
    switch(res_type[i]) {
    case 1:
      res_type[i] = CS_SLES_JACOBI;
      break;
    case 2:
      res_type[i] = CS_SLES_BICGSTAB;
      break;
    case 200:
      res_type[i] = CS_SLES_PCG_SR;
      break;
    default:
      res_type[i] = CS_SLES_PCG;
    }
  }

  cs_multigrid_solve(var_name,
                     res_type[0],
                     res_type[1],
                     res_type[2],
                     true,
                     *ipol,
                     rotation_mode,
                     *iwarnp,
                     *ncymxp,
                     *nitmds,
                     *nitmas,
                     *nitmap,
                     *epsilp,
                     *rnorm,
                     ncyclf,
                     niterf,
                     residu,
                     rhs,
                     vx,
                     0,
                     NULL);

  cs_base_string_f_to_c_free(&var_name);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize multigrid solver API.
 *
 * parameters:
 *   post_cell_max <-- if > 0, activates postprocessing of coarseninsg,
 *                     projecting coarse cell numbers (modulo post_cell_max)
 *                     on the base grid
 *----------------------------------------------------------------------------*/

void
cs_multigrid_initialize(void)
{
}

/*----------------------------------------------------------------------------
 * Finalize multigrid solver API.
 *----------------------------------------------------------------------------*/

void
cs_multigrid_finalize(void)
{
  int ii;

  /* Print system info */

  for (ii = 0; ii < cs_glob_multigrid_n_systems; ii++)
    _multigrid_info_dump(cs_glob_multigrid_systems[ii]);

  /* Free multigrid structures */

  for (ii = 0; ii < cs_glob_multigrid_n_systems; ii++)
    _multigrid_destroy(cs_glob_multigrid_systems + ii);

  BFT_FREE(cs_glob_multigrid_systems);

  cs_glob_multigrid_n_systems = 0;
  cs_glob_multigrid_n_max_systems = 0;

  cs_grid_finalize();
}

/*----------------------------------------------------------------------------
 * Build a hierarchy of meshes starting from a fine mesh, for an
 * ACM (Additive Corrective Multigrid) method.
 *
 * parameters:
 *   var_name               <-- variable name
 *   verbosity              <-- verbosity level
 *   postprocess_block_size <-- if > 0, postprocess coarsening, using
 *                              coarse cell numbers modulo ncpost
 *   aggregation_limit      <-- maximum allowed fine cells per coarse cell
 *   n_max_levels           <-- maximum number of grid levels
 *   n_g_cells_min          <-- global number of cells on coarsest grid
 *                              under which no merging occurs
 *   p0p1_relax             <-- p0/p1 relaxation_parameter
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *   da                     <-- diagonal values (NULL if zero)
 *   xa                     <-- extradiagonal values (NULL if zero)
 *----------------------------------------------------------------------------*/

void
cs_multigrid_build(const char       *var_name,
                   int               verbosity,
                   int               postprocess_block_size,
                   int               aggregation_limit,
                   int               n_max_levels,
                   cs_gnum_t         n_g_cells_min,
                   double            p0p1_relax,
                   bool              symmetric,
                   const int        *diag_block_size,
                   const int        *extra_diag_block_size,
                   const cs_real_t  *da,
                   const cs_real_t  *xa)
{
  cs_timer_t t0, t1, t2;

  int n_coarse_ranks = cs_glob_n_ranks;
  int n_coarse_ranks_prev = 0;
  cs_lnum_t n_cells = 0;
  cs_lnum_t n_cells_with_ghosts = 0;
  cs_lnum_t n_faces = 0;
  cs_gnum_t n_g_cells = 0;
  cs_gnum_t n_g_cells_prev = 0;

  cs_int_t grid_lv = 0;
  cs_multigrid_t *mg = NULL;
  cs_multigrid_level_info_t *mg_lv_info = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;

  cs_grid_t *g = NULL;

  /* Initialization */

  t0 = cs_timer_time();

  mg = _find_or_add_system(var_name);

  if (verbosity > 1)
    bft_printf(_("\n Construction of grid hierarchy for \"%s\"\n"),
               var_name);

  /* Destroy previous hierarchy if necessary */

  if (mg->n_levels > 0) {
    for (grid_lv = mg->n_levels - 1; grid_lv > -1; grid_lv--)
      cs_grid_destroy(mg->grid_hierarchy + grid_lv);
    mg->n_levels = 0;
  }

  /* Build coarse grids hierarchy */
  /*------------------------------*/

  g = cs_grid_create_from_shared(mesh->n_cells,
                                 mesh->n_cells_with_ghosts,
                                 mesh->n_i_faces,
                                 symmetric,
                                 diag_block_size,
                                 extra_diag_block_size,
                                 mesh->i_face_cells,
                                 mesh->halo,
                                 mesh->i_face_numbering,
                                 mq->cell_cen,
                                 mq->cell_vol,
                                 mq->i_face_normal,
                                 da,
                                 xa);

  _multigrid_add_level(mg, g); /* Assign to hierarchy */

  n_cells = mesh->n_cells;
  n_cells_with_ghosts = mesh->n_cells_with_ghosts;
  n_faces = mesh->n_i_faces;
  n_g_cells = mesh->n_g_cells;

  mg_lv_info = mg->lv_info;
  mg_lv_info->n_ranks[0] = cs_glob_n_ranks;
  mg_lv_info->n_elts[0][0] = n_cells;
  mg_lv_info->n_elts[1][0] = n_cells_with_ghosts;
  mg_lv_info->n_elts[2][0] = n_faces;

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t0, &t1);

  while (true) {

    n_g_cells_prev = n_g_cells;
    n_coarse_ranks_prev = n_coarse_ranks;

    /* Recursion test */

    if (grid_lv >= n_max_levels)
      break;

    /* Build coarser grid from previous grid */

    grid_lv += 1;

    if (verbosity > 2)
      bft_printf(_("\n   building level %2d grid\n"), grid_lv);

    g = cs_grid_coarsen(g, verbosity, aggregation_limit, p0p1_relax);

    cs_grid_get_info(g,
                     &grid_lv,
                     &symmetric,
                     NULL,
                     NULL,
                     &n_coarse_ranks,
                     &n_cells,
                     &n_cells_with_ghosts,
                     &n_faces,
                     &n_g_cells);

    _multigrid_add_level(mg, g); /* Assign to hierarchy */

    /* Print coarse mesh stats */

    if (verbosity > 2) {

#if defined(HAVE_MPI)

      if (cs_glob_n_ranks > 1) {

        int lcount[2], gcount[2];
        int n_c_min, n_c_max, n_f_min, n_f_max;

        lcount[0] = n_cells; lcount[1] = n_faces;
        MPI_Allreduce(lcount, gcount, 2, MPI_INT, MPI_MAX,
                      cs_glob_mpi_comm);
        n_c_max = gcount[0]; n_f_max = gcount[1];

        lcount[0] = n_cells; lcount[1] = n_faces;
        MPI_Allreduce(lcount, gcount, 2, MPI_INT, MPI_MIN,
                      cs_glob_mpi_comm);
        n_c_min = gcount[0]; n_f_min = gcount[1];

        bft_printf
          (_("                                  total       min        max\n"
             "     number of cells:     %12llu %10d %10d\n"
             "     number of faces:                  %10d %10d\n"),
           (unsigned long long)n_g_cells, n_c_min, n_c_max, n_f_min, n_f_max);
      }

#endif

      if (cs_glob_n_ranks == 1)
        bft_printf(_("     number of cells:     %10d\n"
                     "     number of faces:     %10d\n"),
                   (int)n_cells, (int)n_faces);

    }

    mg_lv_info = mg->lv_info + grid_lv;
    mg_lv_info->n_ranks[0] = n_coarse_ranks;
    mg_lv_info->n_elts[0][0] = n_cells;
    mg_lv_info->n_elts[1][0] = n_cells_with_ghosts;
    mg_lv_info->n_elts[2][0] = n_faces;

    t2 = cs_timer_time();
    cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t1, &t2);
    t1 = t2;

    /* If too few cells were grouped, we stop at this level */

    if (n_g_cells <= n_g_cells_min)
      break;
    else if (n_g_cells > (0.8 * n_g_cells_prev)
        && n_coarse_ranks == n_coarse_ranks_prev)
      break;
  }

  /* Print final info */

  if (verbosity > 1)
    bft_printf
      (_("   number of coarse grids:           %d\n"
         "   number of cells in coarsest grid: %llu\n\n"),
       grid_lv, (unsigned long long)n_g_cells);

  /* Prepare preprocessing info if necessary */

  if (postprocess_block_size > 0) {

    if (mg->post_cell_max == 0) {
      int mg_id = _multigrid_id(mg);
      if (mg_id > -1)
        cs_post_add_time_dep_output(_cs_multigrid_post_function,
                                    (void *)mg);
      mg->post_cell_max = postprocess_block_size;
    }

    _multigrid_add_post(mg, mesh->n_cells);

  }

  /* Update info */

#if defined(HAVE_MPI)

  /* In parallel, get global (average) values from local values */

  if (cs_glob_n_ranks > 1) {

    int i, j;
    cs_gnum_t *_n_elts_l = NULL, *_n_elts_s = NULL, *_n_elts_m = NULL;

    BFT_MALLOC(_n_elts_l, 3*grid_lv, cs_gnum_t);
    BFT_MALLOC(_n_elts_s, 3*grid_lv, cs_gnum_t);
    BFT_MALLOC(_n_elts_m, 3*grid_lv, cs_gnum_t);

    for (i = 0; i < grid_lv; i++) {
      cs_multigrid_level_info_t *mg_inf = mg->lv_info + i;
      for (j = 0; j < 3; j++)
        _n_elts_l[i*3 + j] = mg_inf->n_elts[j][0];
    }

    MPI_Allreduce(_n_elts_l, _n_elts_s, 3*grid_lv, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
    MPI_Allreduce(_n_elts_l, _n_elts_m, 3*grid_lv, CS_MPI_GNUM, MPI_MAX,
                  cs_glob_mpi_comm);

    for (i = 0; i < grid_lv; i++) {
      cs_multigrid_level_info_t *mg_inf = mg->lv_info + i;
      cs_gnum_t n_g_ranks = mg_inf->n_ranks[0];
      for (j = 0; j < 3; j++) {
        cs_gnum_t tmp_max = n_g_ranks * _n_elts_m[i*3+j];
        mg_inf->n_elts[j][0] = (_n_elts_s[i*3+j] + n_g_ranks/2) / n_g_ranks;
        mg_inf->unbalance[j][0] = (float)(tmp_max*1.0/_n_elts_s[i*3+j]);
      }
    }

    BFT_FREE(_n_elts_m);
    BFT_FREE(_n_elts_s);
    BFT_FREE(_n_elts_l);

  }

#endif

  mg->info.n_levels_tot += grid_lv;

  mg->info.n_levels[0] = grid_lv;

  if (mg->info.n_calls[0] > 0) {
    if (mg->info.n_levels[0] < mg->info.n_levels[1])
      mg->info.n_levels[1] = mg->info.n_levels[0];
    if (mg->info.n_levels[0] > mg->info.n_levels[2])
      mg->info.n_levels[2] = mg->info.n_levels[0];
  }
  else {
    mg->info.n_levels[1] = mg->info.n_levels[0];
    mg->info.n_levels[2] = mg->info.n_levels[0];
  }

  mg->info.n_calls[0] += 1;

  /* Update timers */

  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg->info.t_tot[0]), &t0, &t2);
}

/*----------------------------------------------------------------------------
 * Destroy a hierarchy of meshes starting from a fine mesh, keeping
 * the corresponding system and postprocessing information for future calls.
 *
 * parameters:
 *   var_name <-- variable name
 *----------------------------------------------------------------------------*/

void
cs_multigrid_destroy(const char  *var_name)
{
  int ii;
  cs_timer_t t0, t1;
  cs_multigrid_t *mg = NULL;

  /* Initialization */

  t0 = cs_timer_time();

  mg = _find_or_add_system(var_name);

  /* Destroy grid hierarchy */

  if (mg->n_levels > 0) {
    for (ii = mg->n_levels - 1; ii > -1; ii--)
      cs_grid_destroy(mg->grid_hierarchy + ii);
    mg->n_levels = 0;
  }

  /* Update timers */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg->info.t_tot[0]), &t0, &t1);
}

/*----------------------------------------------------------------------------
 * Sparse linear system resolution using multigrid.
 *
 * parameters:
 *   var_name              <-- variable name
 *   descent_smoother_type <-- type of smoother for descent (PCG, Jacobi, ...)
 *   ascent_smoother_type  <-- type of smoother for ascent (PCG, Jacobi, ...)
 *   coarse_solver_type    <-- type of solver (PCG, Jacobi, ...)
 *   abort_on_divergence   <-- call errorhandler if devergence is detected
 *   poly_degree           <-- preconditioning polynomial degree (0: diagonal)
 *   rotation_mode         <-- halo update option for rotational periodicity
 *   verbosity             <-- verbosity level
 *   n_max_cycles          <-- maximum number of cycles
 *   n_max_iter_descent    <-- maximum nb. of iterations for descent phases
 *   n_max_iter_ascent     <-- maximum nb. of iterations for ascent phases
 *   n_max_iter_coarse     <-- maximum nb. of iterations for coarsest solution
 *   precision             <-- precision limit
 *   r_norm                <-- residue normalization
 *   n_cycles              --> number of cycles
 *   n_equiv_iter          --> number of equivalent iterative solver iterations
 *   residue               <-> residue
 *   rhs                   <-- right hand side
 *   vx                    --> system solution
 *   aux_size              <-- number of elements in aux_vectors
 *   aux_vectors           --- optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   cycle number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

int
cs_multigrid_solve(const char          *var_name,
                   cs_sles_type_t       descent_smoother_type,
                   cs_sles_type_t       ascent_smoother_type,
                   cs_sles_type_t       coarse_solver_type,
                   bool                 abort_on_divergence,
                   int                  poly_degree,
                   cs_halo_rotation_t   rotation_mode,
                   int                  verbosity,
                   int                  n_max_cycles,
                   int                  n_max_iter_descent,
                   int                  n_max_iter_ascent,
                   int                  n_max_iter_coarse,
                   double               precision,
                   double               r_norm,
                   int                 *n_cycles,
                   int                 *n_equiv_iter,
                   double              *residue,
                   const cs_real_t     *rhs,
                   cs_real_t           *vx,
                   size_t               aux_size,
                   void                *aux_vectors)
{
  int ii;
  int db_size[4] = {1, 1, 1, 1};
  int eb_size[4] = {1, 1, 1, 1};

  int cvg = 0;
  cs_lnum_t n_cells = 0;

  cs_multigrid_t *mg = NULL;
  cs_multigrid_info_t *mg_info = NULL;
  cs_timer_t t0, t1;

  t0 = cs_timer_time();
  mg = _find_or_add_system(var_name);
  mg_info = &(mg->info);

  cs_grid_get_info(mg->grid_hierarchy[0],
                   NULL,
                   NULL,
                   db_size,
                   eb_size,
                   NULL,
                   &n_cells,
                   NULL,
                   NULL,
                   NULL);

  /* Initialize number of iterations and residue,
     check for immediate return,
     solve sparse linear system using multigrid algorithm. */

  *n_cycles = 0;
  *n_equiv_iter = 0;

  if (cs_sles_needs_solving(var_name,
                            _("Multigrid"),
                            n_cells*db_size[1],
                            verbosity,
                            r_norm,
                            residue,
                            rhs) != 0) {

    int cycle_id = 1;

    int *n_max_iter = NULL;
    size_t  _aux_size = n_cells * 6 * db_size[1];
    cs_real_t *_aux_vectors = aux_vectors;

    BFT_MALLOC(n_max_iter, mg->n_levels * 2, int);

    if (_aux_size <= aux_size)
      BFT_MALLOC(_aux_vectors, _aux_size, cs_real_t);
    else
      _aux_size = aux_size;

    for (ii = 0; ii < mg->n_levels; ii++) {
      n_max_iter[ii*2] = n_max_iter_descent;
      n_max_iter[ii*2 + 1] = n_max_iter_ascent;
    }
    n_max_iter[(mg->n_levels-1)*2]     = n_max_iter_coarse;
    n_max_iter[(mg->n_levels-1)*2 + 1] = n_max_iter_coarse;

    if (verbosity == 2) /* More detailed headers later if > 2 */
      bft_printf(_("Multigrid [%s]:\n"), var_name);

    /* Cycle to solution */

    while (cvg == 0) {

      if (verbosity > 2)
        bft_printf(_("Multigrid [%s]: cycle %4d\n"),
                   var_name, cycle_id);

      cvg = _multigrid_cycle(mg,
                             descent_smoother_type,
                             ascent_smoother_type,
                             coarse_solver_type,
                             abort_on_divergence,
                             poly_degree,
                             rotation_mode,
                             verbosity,
                             cycle_id,
                             n_max_cycles,
                             n_max_iter,
                             n_equiv_iter,
                             precision,
                             r_norm,
                             residue,
                             rhs,
                             vx,
                             aux_size,
                             _aux_vectors);

      cycle_id++;
      *n_cycles += 1;
    }

    if (_aux_vectors != aux_vectors)
      BFT_FREE(_aux_vectors);
    BFT_FREE(n_max_iter);
  }

  /* Update statistics */

  t1 = cs_timer_time();

  /* Update stats on number of iterations (last, min, max, total) */

  mg_info->type[0] = descent_smoother_type;
  mg_info->type[1] = ascent_smoother_type;
  mg_info->type[2] = coarse_solver_type;

  mg_info->n_cycles[2] += *n_cycles;

  if (mg_info->n_calls[1] > 0) {
    if (mg_info->n_cycles[0] > (unsigned)(*n_cycles))
      mg_info->n_cycles[0] = *n_cycles;
    if (mg_info->n_cycles[1] < (unsigned)(*n_cycles))
      mg_info->n_cycles[1] = *n_cycles;
  }
  else {
    mg_info->n_cycles[0] = *n_cycles;
    mg_info->n_cycles[1] = *n_cycles;
  }

  /* Update number of resolutions and timing data */

  mg_info->n_calls[1] += 1;
  cs_timer_counter_add_diff(&(mg->info.t_tot[1]), &t0, &t1);

  return cvg;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
