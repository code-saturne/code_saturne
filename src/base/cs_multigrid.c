/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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
 * Multigrid solver.
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(__STDC_VERSION__)      /* size_t */
#if (__STDC_VERSION__ == 199901L)
#    include <stddef.h>
#  else
#    include <stdlib.h>
#  endif
#else
#include <stdlib.h>
#endif

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_grid.h"
#include "cs_matrix.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_perio.h"
#include "cs_post.h"
#include "cs_sles.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_multigrid.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define EPZERO  1.E-12
#define RINFIN  1.E+30

#if !defined(HUGE_VAL)
#define HUGE_VAL  1.E+12
#endif


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

  unsigned             n_builds;            /* Number of times grids built */
  unsigned             n_solves;            /* Number of times system solved */

  unsigned long long   n_levels_tot;        /* Total accumulated number of
                                               grid levels built */
  unsigned             n_levels[3];         /* Number of grid levels:
                                               [last, min, max] */

  unsigned             n_iterations[3][4];  /* Number of iterations for
                                               system resolution:
                                                 [last, min, max]
                                                 [finest, coarsest, total,
                                                  fine grid equivalent] */
  unsigned long long   n_iterations_tot[4]; /* Total accumulated number of
                                               iterations:
                                                 [finest, coarsest, total,
                                                  fine grid equivalent] */

  double               wt_tot[2];           /* Total wall-clock time used:
                                                 [build, solve] */
  double               cpu_tot[2];          /* Total (local) CPU used:
                                                 [build, solve] */

} cs_multigrid_info_t;

/* Grid hierarchy */
/*----------------*/

typedef struct _cs_multigrid_t {

  cs_multigrid_info_t info;   /* Multigrid info */

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
  int i, j;

  BFT_MALLOC(info->name, strlen(name) + 1, char);

  strcpy(info->name, name);

  for (i = 0; i < 3; i++)
    info->type[i] = CS_SLES_N_TYPES;

  info->n_builds = 0;
  info->n_solves = 0;

  info->n_levels_tot = 0;

  for (i = 0; i < 3; i++) {
    info->n_levels[i] = 0;
    for (j = 0; j < 4; j++)
      info->n_iterations[i][j] = 0;
  }

  for (i = 0; i < 4; i++)
    info->n_iterations_tot[i] = 0;

  for (i = 0; i < 2; i++) {
    info->wt_tot[i] = 0.0;
    info->cpu_tot[i] = 0.0;
  }
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
 * Output information regarding multigrid resolution.
 *
 * parameters:
 *   this_info <-> pointer to multigrid info structure
 *----------------------------------------------------------------------------*/

static void
_multigrid_info_dump(const cs_multigrid_info_t *this_info)
{
  unsigned long long n_builds_denom = CS_MAX(this_info->n_builds, 1);
  unsigned long long n_solves_denom = CS_MAX(this_info->n_solves, 1);
  int n_builds = this_info->n_builds;
  int n_solves = this_info->n_solves;
  int n_lv_min = this_info->n_levels[1];
  int n_lv_max = this_info->n_levels[2];
  int n_lv_mean = (int)(this_info->n_levels_tot / n_builds_denom);
  int n_it_f_min = this_info->n_iterations[1][0];
  int n_it_f_max = this_info->n_iterations[2][0];
  int n_it_c_min = this_info->n_iterations[1][1];
  int n_it_c_max = this_info->n_iterations[2][1];
  int n_it_t_min = this_info->n_iterations[1][2];
  int n_it_t_max = this_info->n_iterations[2][2];
  int n_it_e_min = this_info->n_iterations[1][3];
  int n_it_e_max = this_info->n_iterations[2][3];
  int n_it_f_mean = (int)(this_info->n_iterations_tot[0] / n_solves_denom);
  int n_it_c_mean = (int)(this_info->n_iterations_tot[1] / n_solves_denom);
  int n_it_t_mean = (int)(this_info->n_iterations_tot[2] / n_solves_denom);
  int n_it_e_mean = (int)(this_info->n_iterations_tot[3] / n_solves_denom);

  bft_printf(_("\n"
               "Summary of multigrid for \"%s\":\n\n"),
               this_info->name);

  if (this_info->type[0] != CS_SLES_N_TYPES) {

    const char *descent_smoother_name = cs_sles_type_name[this_info->type[0]];
    const char *ascent_smoother_name = cs_sles_type_name[this_info->type[1]];

    if (this_info->type[0] == this_info->type[1])
      bft_printf(_("  Smoother: %s\n"), _(descent_smoother_name));
    else
      bft_printf(_("  Descent smoother:     %s\n"
                   "  Ascent smoother:      %s\n"),
                 _(descent_smoother_name), _(ascent_smoother_name));

    bft_printf(_("  Coarsest level solver:       %s\n"),
               _(cs_sles_type_name[this_info->type[2]]));

  }

  bft_printf(_("  Number of constructions:          %d\n"
               "  Number of resolutions:            %d\n"
               "  Number of levels:\n"
               "    minimum:                        %d\n"
               "    maximum:                        %d\n"
               "    mean:                           %d\n"
               "  Number of iterations:\n"
               "    on finest grid:\n"
               "      minimum:                      %d\n"
               "      maximum:                      %d\n"
               "      mean:                         %d\n"
               "    on coarsest grid:\n"
               "      minimum:                      %d\n"
               "      maximum:                      %d\n"
               "      mean:                         %d\n"
               "    total on grids:\n"
               "      minimum:                      %d\n"
               "      maximum:                      %d\n"
               "      mean:                         %d\n"
               "    equivalent (total weighted by number of cells) :\n"
               "      minimum:                      %d\n"
               "      maximum:                      %d\n"
               "      mean:                         %d\n"
               "  Associated times (construction, resolution)\n"
               "    total elapsed:                  %12.3f  %12.3f\n"),
             n_builds, n_solves, n_lv_min, n_lv_max, n_lv_mean,
             n_it_f_min, n_it_f_max, n_it_f_mean,
             n_it_c_min, n_it_c_max, n_it_c_mean,
             n_it_t_min, n_it_t_max, n_it_t_mean,
             n_it_e_min, n_it_e_max, n_it_e_mean,
             this_info->wt_tot[0], this_info->wt_tot[1]);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    double cpu_min[2], cpu_max[2], cpu_tot[2], cpu_loc[2];
    cpu_loc[0] = this_info->cpu_tot[0];
    cpu_loc[1] = this_info->cpu_tot[1];

    MPI_Allreduce(cpu_loc, cpu_min, 2, MPI_DOUBLE, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(cpu_loc, cpu_max, 2, MPI_DOUBLE, MPI_MAX,
                  cs_glob_mpi_comm);
    MPI_Allreduce(cpu_loc, cpu_tot, 2, MPI_DOUBLE, MPI_SUM,
                  cs_glob_mpi_comm);

    bft_printf(_("    Min local total CPU time:       %12.3f  %12.3f\n"
                 "    Max local total CPU time:       %12.3f  %12.3f\n"
                 "    Total CPU time:                 %12.3f  %12.3f\n"),
               cpu_min[0], cpu_min[1], cpu_max[0], cpu_max[1],
               cpu_tot[0], cpu_tot[1]);

  }

#endif

  if (cs_glob_n_ranks == 1)
    bft_printf(_("    Total CPU time:                 %12.3f  %12.3f\n"),
               this_info->cpu_tot[0], this_info->cpu_tot[1]);
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

  for (ii = 0; ii < mg->n_levels_max; ii++)
    mg->grid_hierarchy[ii] = NULL;

  mg->post_cell_num = NULL;

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

  _multigrid_info_unset(&(_mg->info));

  for (ii = 0; ii < _mg->n_levels_max; ii++)
    cs_grid_destroy(_mg->grid_hierarchy + ii);

  if (_mg->post_cell_max > 0) {
    for (ii = 0; ii < _mg->n_levels_max - 1; ii++)
      if (_mg->post_cell_num[ii] != NULL)
        BFT_FREE(_mg->post_cell_num[ii]);
    BFT_FREE(_mg->post_cell_num);
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
    for (ii = mg->n_levels; ii < mg->n_levels_max; ii++)
      mg->grid_hierarchy[ii] = NULL;

    if (mg->post_cell_num != NULL) {
      BFT_REALLOC(mg->post_cell_num, mg->n_levels_max, int *);
      for (ii = mg->n_levels; ii < mg->n_levels_max; ii++)
        mg->post_cell_num[ii] = NULL;
      if (mg->n_levels > 0)
        mg->post_cell_num[mg->n_levels - 1] = NULL;
    }
  }

  mg->grid_hierarchy[mg->n_levels] = grid;
  mg->n_levels += 1;
}

/*----------------------------------------------------------------------------
 * Get multigrid structure's id in list of know systems
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
                    fvm_lnum_t       n_base_cells)
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

  for (ii = 0; ii < mg->n_levels_post; ii++) {
    BFT_REALLOC(mg->post_cell_num[ii], n_base_cells, int);
    cs_grid_project_cell_num(mg->grid_hierarchy[ii+1],
                             n_base_cells,
                             mg->post_cell_max,
                             mg->post_cell_num[ii]);
  }
}

/*----------------------------------------------------------------------------
 * Post process variables associated with Syrthes couplings
 *
 * parameters:
 *   hierarchy_id        <--  Id of multigrid hierarchy
 *   nt_cur_abs          <--  Current time step
 *   t_cur_abs           <--  Current time value
 *----------------------------------------------------------------------------*/

static void
_cs_multigrid_post_function(cs_int_t   hierarchy_id,
                            cs_int_t   nt_cur_abs,
                            cs_real_t  t_cur_abs)
{
  int ii;
  size_t name_len;
  char *var_name = NULL;
  cs_multigrid_t *mg = NULL;
  const char name_prefix[] = "mg";
  const char *base_name = NULL;

  /* Return if necessary structures inconsistant or have been destroyed */

  if (hierarchy_id < cs_glob_multigrid_n_systems)
    mg = cs_glob_multigrid_systems[hierarchy_id];

  if (mg == NULL)
    return;

  if (mg->post_cell_num == NULL || cs_post_mesh_exists(-1) != true)
    return;

  /* Allocate name buffer */

  base_name = mg->info.name;
  name_len = strlen(name_prefix) + 1 + strlen(base_name) + 1 + 3 + 1 + 4 + 1;
  BFT_MALLOC(var_name, name_len, char);

  /* Loop on grid levels */

  for (ii = 0; ii < mg->n_levels_post; ii++) {

    sprintf(var_name, "%s %s %3d %2d",
            name_prefix, base_name, (int)(ii+1), (int)nt_cur_abs);

    cs_post_write_var(-1,
                      var_name,
                      1,
                      false,
                      true,
                      CS_POST_TYPE_int,
                      -1,
                      0.0,
                      mg->post_cell_num[ii],
                      NULL,
                      NULL);

    BFT_FREE(mg->post_cell_num[ii]);

  }
  mg->n_levels_post = 0;

  BFT_FREE(var_name);
}

/*----------------------------------------------------------------------------
 * Compute dot product, summing result over all ranks.
 *
 * parameters:
 *   n_elts <-- Local number of elements
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
  double s = cblas_ddot(n_elts, x, 1, y, 1);

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
 *   var_name      <-- Variable name
 *   n_f_cells     <-- Number of cells on fine mesh
 *   n_max_cycles  <-- Maximum number of cycles
 *   cycle_id      <-- Number of current cycle
 *
 *   verbosity     <-- Verbosity level
 *   n_iters       <-- Number of iterations
 *   precision     <-- Precision limit
 *   r_norm        <-- Residue normalization
 *   residue       <-> Residue
 *   rhs           --> Right-hand side
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   cycle number reached.
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
                  cs_real_t          *rhs)
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

  /* Compute residue */

  *residue = sqrt(_dot_product(n_f_cells, rhs, rhs));

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

  else if (verbosity > 2)
    bft_printf(_(cycle_fmt), cycle_id, n_iters, *residue/r_norm);

  return 0;
}

/*----------------------------------------------------------------------------
 * Sparse linear system resolution using multigrid.
 *
 * parameters:
 *   mg                    <-- Multigrid system
 *   descent_smoother_type <-- Type of smoother for descent (PCG, Jacobi, ...)
 *   ascent_smoother_type  <-- Type of smoother for ascent (PCG, Jacobi, ...)
 *   coarse_solver_type    <-- Type of solver (PCG, Jacobi, ...)
 *   symmetric             <-- Symmetric coefficients indicator
 *   poly_degree           <-- Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode         <-- Halo update option for rotational periodicity
 *   verbosity             <-- Verbosity level
 *   cycle_id              <-- Id of currect cycle
 *   n_max_cycles          <-- Maximum number of cycles
 *   n_max_iter            <-- Maximum number of iterations per grid level
 *                             n_max_iter[level * 2]     for descent
 *                             n_max_iter[level * 2 + 1] for ascent
 *   precision             <-- Precision limit
 *   r_norm                <-- Residue normalization
 *   n_level_iter          <-> Number of iterations per level
 *   residue               <-> Residue
 *   rhs                   <-- Right hand side
 *   vx                    --> System solution
 *   aux_size              <-- Number of elements in aux_vectors
 *   aux_vectors           --- Optional working area (allocation otherwise)
 *
 * Returns
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   cycle number reached.
 *----------------------------------------------------------------------------*/

static int
_multigrid_cycle(cs_multigrid_t     *mg,
                 cs_sles_type_t      descent_smoother_type,
                 cs_sles_type_t      ascent_smoother_type,
                 cs_sles_type_t      coarse_solver_type,
                 cs_bool_t           symmetric,
                 int                 poly_degree,
                 cs_perio_rota_t     rotation_mode,
                 int                 verbosity,
                 int                 cycle_id,
                 int                 n_max_cycles,
                 const int           n_max_iter[],
                 double              precision,
                 double              r_norm,
                 int                 n_level_iter[],
                 double             *residue,
                 const cs_real_t    *rhs,
                 cs_real_t          *vx,
                 size_t              aux_size,
                 void               *aux_vectors)
{
  int level, coarsest_level;
  fvm_lnum_t ii;

  int cvg = 0;
  int n_iter = 0;
  size_t alloc_size = 0;
  cs_real_t c_precision = precision;
  cs_real_t _residue = -1.;

  size_t _aux_size = aux_size;
  fvm_lnum_t n_cells = 0, n_cells_ext = 0;
  cs_real_t r_norm_l = r_norm;

  char _var_lv_name[33];
  char *var_lv_name = _var_lv_name;

  cs_real_t *_aux_vectors = aux_vectors;
  cs_real_t *wr = NULL;
  cs_real_t *_rhs_vx_val = NULL;
  cs_real_t **_rhs_vx = NULL, **_rhs = NULL, **_vx = NULL;
  cs_matrix_t  *_matrix;
  cs_multigrid_info_t *mg_info = NULL;

  const char *var_name = NULL;
  const cs_real_t *_rhs_level = NULL;
  const cs_real_t *_da = NULL, *_xa = NULL;
  const cs_grid_t *f = NULL, *c= NULL;

  cs_bool_t end_cycle = false;

  /* Initialization */

  mg_info = &(mg->info);
  var_name = mg_info->name;

  coarsest_level = mg->n_levels - 1;

  /* In theory, one should increase precision on coarsest mesh,
     but in practice, it more efficient to have a lower precision */
  /* c_precision = precision * 0.01; */
  c_precision = precision;

  f = mg->grid_hierarchy[0];

  cs_grid_get_info(f,
                   NULL,
                   NULL,
                   &n_cells,
                   &n_cells_ext,
                   NULL,
                   NULL);

  if (strlen(var_name) + 5 > 32)
    BFT_MALLOC(var_lv_name, strlen(var_name) + 5 + 1, char);

  /* Allocate wr or use working area */

  if (aux_size >= (size_t)n_cells_ext) {
    wr = aux_vectors;
    _aux_vectors = wr + n_cells_ext;
    _aux_size = aux_size - n_cells_ext;
  }
  else
    BFT_MALLOC(wr, n_cells_ext, cs_real_t);

  /* reserve memory for rhs and vx;
     for the finest level, simply point to input and output arrays */

  BFT_MALLOC(_rhs_vx, mg->n_levels*2, cs_real_t *);
  _rhs = _rhs_vx;
  _vx = _rhs_vx + mg->n_levels;

  _rhs[0] = NULL; /* Use _rhs_level when necessary to avoid const warning */
  _vx[0] = vx;

  /* Reserve memory for corrections and residues for coarse levels */

  if (mg->n_levels > 1) {

    alloc_size = 0;

    for (level = 1; level < mg->n_levels; level++)
      alloc_size += cs_grid_get_n_cells_ext(mg->grid_hierarchy[level]);

    BFT_MALLOC(_rhs_vx_val, alloc_size*2, cs_real_t);

    _rhs[1] = _rhs_vx_val;
    _vx[1] = _rhs_vx_val + alloc_size;

    for (level = 2; level < mg->n_levels; level++) {
      fvm_lnum_t _n_cells_ext_prev
        = cs_grid_get_n_cells_ext(mg->grid_hierarchy[level-1]);
      _rhs[level] = _rhs[level - 1] + _n_cells_ext_prev;
      _vx[level] = _vx[level - 1] + _n_cells_ext_prev;
    }
  }

  /* Descent */
  /*---------*/

  if (verbosity > 2)
    bft_printf(_("  Multigrid cycle: descent\n"));

  for (level = 0; level < coarsest_level; level++) {

    _rhs_level = (level == 0) ?  rhs : _rhs[level];

    sprintf(var_lv_name, "%s:%04d", var_name, level);

    c = mg->grid_hierarchy[level+1];

    /* Smoother pass */

    if (verbosity > 2)
      bft_printf(_("    level %3d: smoother\n"), level);

    cs_grid_get_matrix(f, &_da, &_xa, &_matrix);

    cs_sles_solve(var_lv_name,
                  descent_smoother_type,
                  false, /* Stats not updated here */
                  symmetric,
                  _da,
                  _xa,
                  _matrix,
                  NULL,
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

    n_level_iter[level] += n_iter;

    /* Restrict residue
       TODO: get residue from cs_sles_solve(). This optimisation would
       require adding an argument and exercising caution to ensure the
       correct sign and meaning of the residue. */

    cs_matrix_set_coefficients(_matrix, symmetric, _da, _xa);

    cs_matrix_vector_multiply(rotation_mode,
                              _matrix,
                              _vx[level],
                              wr);

    for (ii = 0; ii < n_cells; ii++)
      wr[ii] = _rhs_level[ii] - wr[ii];

    n_level_iter[level] += 1;

    /* Convergence test in beginning of cycle (fine mesh) */

    if (level == 0) {

      cvg = _convergence_test(var_name,
                              n_cells,
                              n_max_cycles,
                              cycle_id,
                              verbosity,
                              n_level_iter[0],
                              precision,
                              r_norm,
                              residue,
                              wr);

      /* If converged or cycle limit reached, break from descent loop */

      if (cvg != 0) {
        end_cycle = true;
        break;
      }

    }

    /* Prepare for next level */

    cs_grid_restrict_cell_var(f, c, wr, _rhs[level+1]);

    cs_grid_get_info(c,
                     NULL,
                     NULL,
                     &n_cells,
                     &n_cells_ext,
                     NULL,
                     NULL);

    f = c;

    for (ii = 0; ii < n_cells; ii++) /* Initialize correction */
      _vx[level+1][ii] = 0.0;

  } /* End of loop on levels (descent) */

  if (end_cycle == false) {

    /* Resolve coarsest level to convergence */
    /*---------------------------------------*/

    if (verbosity > 2)
      bft_printf(_("  Resolution on coarsest level\n"));

    assert(level = coarsest_level);
    assert(c == mg->grid_hierarchy[coarsest_level]);

    /* coarsest level == 0 should never happen, but we play it safe */
    _rhs_level = (level == 0) ?  rhs : _rhs[coarsest_level];

    sprintf(var_lv_name, "%s:%04d", var_name, coarsest_level);

    cs_grid_get_matrix(c, &_da, &_xa, &_matrix);

    cs_sles_solve(var_lv_name,
                  coarse_solver_type,
                  false, /* Stats not updated here */
                  symmetric,
                  _da,
                  _xa,
                  _matrix,
                  NULL,
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

    n_level_iter[level] += n_iter;

    /* Ascent */
    /*--------*/

    if (verbosity > 2)
      bft_printf(_("  Multigrid cycle: ascent\n"));

    for (level = coarsest_level - 1; level > -1; level--) {

      cs_real_t *_f_vx = _vx[level];

      c = mg->grid_hierarchy[level+1];
      f = mg->grid_hierarchy[level];

      cs_grid_get_info(f,
                       NULL,
                       NULL,
                       &n_cells,
                       &n_cells_ext,
                       NULL,
                       NULL);

      /* Prolong correction */

      cs_grid_prolong_cell_var(c, f, _vx[level+1], wr);

      for (ii = 0; ii < n_cells; ii++)
        _f_vx[ii] += wr[ii];

      /* Smoother pass if level > 0
         (smoother not called for finest mesh, as it will be called in
         descent phase of the next cycle, before the convergence test). */

      if (level > 0) {

        if (verbosity > 2)
          bft_printf(_("    level %3d: smoother\n"), level);

        sprintf(var_lv_name, "%s:%04d", var_name, level);

        cs_grid_get_matrix(f, &_da, &_xa, &_matrix);

        cs_sles_solve(var_lv_name,
                      ascent_smoother_type,
                      false, /* Stats not updated here */
                      symmetric,
                      _da,
                      _xa,
                      _matrix,
                      NULL,
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

        n_level_iter[level] += n_iter;
      }

    } /* End loop on levels (ascent) */

  } /* End of test on end_cycle */

  /* Free memory */

  if (var_lv_name != _var_lv_name)
    BFT_FREE(var_lv_name);

  if (wr != aux_vectors)
    BFT_FREE(wr);

  BFT_FREE(_rhs_vx);
  BFT_FREE(_rhs_vx_val);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Sparse linear system resolution using multigrid.
 *
 * parameters:
 *   var_name              <-- Variable name
 *   descent_smoother_type <-- Type of smoother for descent (PCG, Jacobi, ...)
 *   ascent_smoother_type  <-- Type of smoother for ascent (PCG, Jacobi, ...)
 *   coarse_solver_type    <-- Type of solver (PCG, Jacobi, ...)
 *   symmetric             <-- Symmetric coefficients indicator
 *   poly_degree           <-- Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode         <-- Halo update option for rotational periodicity
 *   verbosity             <-- Verbosity level
 *   n_max_cycles          <-- Maximum number of cycles
 *   n_max_iter_descent    <-- Maximum nb. of iterations for descent phases
 *   n_max_iter_ascent     <-- Maximum nb. of iterations for ascent phases
 *   n_max_iter_coarse     <-- Maximum nb. of iterations for coarsest solution
 *   precision             <-- Precision limit
 *   r_norm                <-- Residue normalization
 *   n_cycles              --> Number of cycles
 *   n_iter                --> Number of iterations
 *   residue               <-> Residue
 *   rhs                   <-- Right hand side
 *   vx                    --> System solution
 *   aux_size              <-- Number of elements in aux_vectors
 *   aux_vectors           --- Optional working area (allocation otherwise)
 *----------------------------------------------------------------------------*/

static void
_multigrid_solve(const char         *var_name,
                 cs_sles_type_t      descent_smoother_type,
                 cs_sles_type_t      ascent_smoother_type,
                 cs_sles_type_t      coarse_solver_type,
                 cs_bool_t           symmetric,
                 int                 poly_degree,
                 cs_perio_rota_t     rotation_mode,
                 int                 verbosity,
                 int                 n_max_cycles,
                 int                 n_max_iter_descent,
                 int                 n_max_iter_ascent,
                 int                 n_max_iter_coarse,
                 double              precision,
                 double              r_norm,
                 int                *n_cycles,
                 int                *n_iter,
                 double             *residue,
                 const cs_real_t    *rhs,
                 cs_real_t          *vx,
                 size_t              aux_size,
                 void               *aux_vectors)
{
  int ii;
  unsigned _n_iter[4] = {0, 0, 0, 0};

  fvm_lnum_t n_cells = 0;

  cs_multigrid_t *mg = NULL;
  cs_multigrid_info_t *mg_info = NULL;
  double  wt_start = 0.0, wt_stop = 0.0;
  double  cpu_start = 0.0, cpu_stop = 0.0;

  wt_start =bft_timer_wtime();
  cpu_start =bft_timer_cpu_time();
  mg = _find_or_add_system(var_name);
  mg_info = &(mg->info);

  cs_grid_get_info(mg->grid_hierarchy[0],
                   NULL,
                   NULL,
                   &n_cells,
                   NULL,
                   NULL,
                   NULL);

  /* Initialize number of iterations and residue,
     check for immediate return,
     solve sparse linear system using multigrid algorithm. */

  *n_cycles = 0;
  *n_iter = 0;

  if (cs_sles_needs_solving(var_name,
                            _("Multigrid"),
                            n_cells,
                            verbosity,
                            r_norm,
                            residue,
                            rhs) != 0) {

    int cycle_id = 1, cvg = 0;
    double it_count_num = 0.0;

    int *n_max_iter = NULL;
    int *n_level_iter = NULL;
    size_t  _aux_size = n_cells * 6;
    cs_real_t *_aux_vectors = aux_vectors;

    BFT_MALLOC(n_max_iter, mg->n_levels * 2, int);
    BFT_MALLOC(n_level_iter, mg->n_levels, int);

    if (_aux_size <= aux_size)
      BFT_MALLOC(_aux_vectors, _aux_size, cs_real_t);
    else
      _aux_size = aux_size;

    for (ii = 0; ii < mg->n_levels; ii++) {
      n_max_iter[ii*2] = n_max_iter_descent;
      n_max_iter[ii*2 + 1] = n_max_iter_ascent;
      n_level_iter[ii] = 0;
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
                             symmetric,
                             poly_degree,
                             rotation_mode,
                             verbosity,
                             cycle_id,
                             n_max_cycles,
                             n_max_iter,
                             precision,
                             r_norm,
                             n_level_iter,
                             residue,
                             rhs,
                             vx,
                             aux_size,
                             _aux_vectors);

      cycle_id++;
      *n_cycles += 1;
    }

    _n_iter[0] = n_level_iter[0];
    _n_iter[1] = n_level_iter[mg->n_levels - 1];

    for (ii = 0; ii < mg->n_levels; ii++)
      _n_iter[2] += n_level_iter[ii];

    /* Estimate "equivalent" iterations */

    for (ii = 0; ii < mg->n_levels; ii++) {
      fvm_gnum_t n_g_cells = cs_grid_get_n_g_cells(mg->grid_hierarchy[ii]);
      it_count_num += n_g_cells * n_level_iter[ii];
      _n_iter[2] += n_level_iter[ii];
    }

    _n_iter[3] = (int)(  it_count_num
                       / cs_grid_get_n_g_cells(mg->grid_hierarchy[0]));
    *n_iter = _n_iter[3];

    if (_aux_vectors != aux_vectors)
      BFT_FREE(_aux_vectors);
    BFT_FREE(n_level_iter);
    BFT_FREE(n_max_iter);
  }

  /* Update statistics */

  wt_stop =bft_timer_wtime();
  cpu_stop =bft_timer_cpu_time();

  /* Update stats on number of iterations (last, min, max, total) */

  mg_info->type[0] = descent_smoother_type;
  mg_info->type[1] = ascent_smoother_type;
  mg_info->type[2] = coarse_solver_type;

  for (ii = 0; ii < 4; ii++)
    mg_info->n_iterations[0][ii] = _n_iter[ii];

  if (mg_info->n_solves > 0) {
    for (ii = 0; ii < 4; ii++) {
      if (mg_info->n_iterations[1][ii] > _n_iter[ii])
        mg_info->n_iterations[1][ii] = _n_iter[ii];
      if (mg_info->n_iterations[2][ii] < _n_iter[ii])
        mg_info->n_iterations[2][ii] = _n_iter[ii];
    }
  }
  else {
    for (ii = 0; ii < 4; ii++) {
      mg_info->n_iterations[1][ii] = _n_iter[ii];
      mg_info->n_iterations[2][ii] = _n_iter[ii];
    }
  }

  for (ii = 0; ii < 4; ii++)
    mg_info->n_iterations_tot[ii] += _n_iter[ii];

  /* Update number of resolutions and timing data */

  mg_info->n_solves += 1;

  mg_info->wt_tot[1] += (wt_stop - wt_start);
  mg_info->cpu_tot[1] += (cpu_stop - cpu_start);
}

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Build a hierarchy of meshes starting from a fine mesh, for an
 * ACM (Additive Corrective Multigrid) method, grouping cells at
 * most 2 by 2.
 *----------------------------------------------------------------------------*/

void CS_PROCF(clmlga, CLMLGA)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname,     /* <-- variable name length */
 const cs_int_t   *ncelet,    /* <-- Number of cells, halo included */
 const cs_int_t   *ncel,      /* <-- Number of local cells */
 const cs_int_t   *nfac,      /* <-- Number of internal faces */
 const cs_int_t   *isym,      /* <-- Symmetry indicator:
                                     1: symmetric; 2: not symmetric */
 cs_int_t         *iagmax,    /* <-> Maximum agglomeration count */
 const cs_int_t   *nagmax,    /* <-- Agglomeration count limit */
 const cs_int_t   *ncpost,    /* <-- If > 0, postprocess coarsening, using
                                     coarse cell numbers modulo ncpost */
 const cs_int_t   *iwarnp,    /* <-- Verbosity level */
 const cs_int_t   *ngrmax,    /* <-- Maximum number of grid levels */
 const cs_int_t   *ncegrm,    /* <-- Maximum local number of cells on
                                     coarsest grid */
 const cs_real_t  *dam,       /* <-- Matrix diagonal */
 const cs_real_t  *xam        /* <-- Matrix extra-diagonal terms */
)
{
  char *var_name;
  double  wt_start, wt_stop, cpu_start, cpu_stop;

  fvm_lnum_t n_cells = 0;
  fvm_lnum_t n_faces = 0;
  fvm_gnum_t n_g_cells = 0;
  fvm_gnum_t n_g_cells_prev = 0;

  cs_int_t grid_lv = 0;
  cs_multigrid_t *mg = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;

  cs_grid_t *g = NULL;

  cs_bool_t symmetric = (*isym == 1) ? true : false;

  assert(*ncelet >= *ncel);
  assert(*nfac > 0);

  /* Initialization */

  wt_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

  var_name = cs_base_string_f_to_c_create(cname, *lname);

  mg = _find_or_add_system(var_name);

  if (*iwarnp > 1)
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
                                 mesh->i_face_cells,
                                 mesh->halo,
                                 mesh->i_face_numbering,
                                 mq->cell_cen,
                                 mq->cell_vol,
                                 mq->i_face_normal,
                                 dam,
                                 xam);

  _multigrid_add_level(mg, g); /* Assign to hierarchy */

  n_cells = mesh->n_cells;
  n_faces = mesh->n_i_faces;
  n_g_cells = mesh->n_g_cells;

  while (true) {

    n_g_cells_prev = n_g_cells;

    /* Recursion test */

    if (n_g_cells <= (fvm_gnum_t)(*ncegrm))
      break;

    else if (grid_lv >= *ngrmax) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_(" CLMLGA: maximum number of coarse grids (%d)\n"
                   "         reached for \"%s\".\n"),
                 (int)(*ngrmax), var_name);
      break;
    }

    /* Build coarser grid from previous grid */

    grid_lv += 1;

    if (*iwarnp > 2)
      bft_printf(_("\n   building level %2d grid\n"), grid_lv);

    g = cs_grid_coarsen(g, *iwarnp, *nagmax, iagmax);

    cs_grid_get_info(g,
                     &grid_lv,
                     &symmetric,
                     &n_cells,
                     NULL,
                     &n_faces,
                     &n_g_cells);

    _multigrid_add_level(mg, g); /* Assign to hierarchy */

    /* Print coarse mesh stats */

    if (*iwarnp > 2) {

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
             "     number of cells:     %12lu %10d %10d\n"
             "     number of faces:                  %10d %10d\n"),
           (unsigned long)n_g_cells, n_c_min, n_c_max, n_f_min, n_f_max);
      }

#endif

      if (cs_glob_n_ranks == 1)
        bft_printf(_("     number of cells:     %10d\n"
                     "     number of faces:     %10d\n"),
                   (int)n_cells, (int)n_faces);

    }

    /* If too few cells were grouped, we stop at this level */

    if (   n_g_cells > (0.8 * n_g_cells_prev)
        || n_g_cells < (1.5 * cs_glob_n_ranks))
      break;
  }

  cs_base_string_f_to_c_free(&var_name);

  /* Print final info */

  if (*iwarnp > 1)
    bft_printf
      (_("   number of coarse grids:           %d\n"
         "   number of cells in coarsest grid: %lu\n\n"),
       grid_lv, (unsigned long)n_g_cells);

  /* Prepare preprocessing info if necessary */

  if (*ncpost > 0) {

    if (mg->post_cell_max == 0) {
      int mg_id = _multigrid_id(mg);
      if (mg_id > -1)
        cs_post_add_time_dep_var(_cs_multigrid_post_function, mg_id);
      mg->post_cell_max = *ncpost;
    }

    _multigrid_add_post(mg, mesh->n_cells);

  }

  /* Update info */

  mg->info.n_levels_tot += grid_lv;

  mg->info.n_levels[0] = grid_lv;

  if (mg->info.n_builds > 0) {
    if (mg->info.n_levels[0] < mg->info.n_levels[1])
      mg->info.n_levels[1] = mg->info.n_levels[0];
    if (mg->info.n_levels[0] > mg->info.n_levels[2])
      mg->info.n_levels[2] = mg->info.n_levels[0];
  }
  else {
    mg->info.n_levels[1] = mg->info.n_levels[0];
    mg->info.n_levels[2] = mg->info.n_levels[0];
  }

  mg->info.n_builds += 1;

  /* Update timers */

  wt_stop = bft_timer_wtime();
  cpu_stop = bft_timer_cpu_time();

  mg->info.wt_tot[0] += (wt_stop - wt_start);
  mg->info.cpu_tot[0] += (cpu_stop - cpu_start);
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
  int ii;
  double  wt_start, wt_stop, cpu_start, cpu_stop;
  cs_multigrid_t *mg = NULL;

  /* Initialization */

  wt_start =bft_timer_wtime();
  cpu_start =bft_timer_cpu_time();

  var_name = cs_base_string_f_to_c_create(cname, *lname);

  mg = _find_or_add_system(var_name);

  cs_base_string_f_to_c_free(&var_name);

  /* Destroy grid hierarchy */

  if (mg->n_levels > 0) {
    for (ii = mg->n_levels - 1; ii > -1; ii--)
      cs_grid_destroy(mg->grid_hierarchy + ii);
    mg->n_levels = 0;
  }

  /* Update timers */

  wt_stop = bft_timer_wtime();
  cpu_stop = bft_timer_cpu_time();

  mg->info.wt_tot[0] += (wt_stop - wt_start);
  mg->info.cpu_tot[0] += (cpu_stop - cpu_start);
}

/*----------------------------------------------------------------------------
 * General sparse linear system resolution
 *----------------------------------------------------------------------------*/

void CS_PROCF(resmgr, RESMGR)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname,     /* <-- variable name length */
 const cs_int_t   *ncelet,    /* <-- Number of cells, halo included */
 const cs_int_t   *ncel,      /* <-- Number of local cells */
 const cs_int_t   *nfac,      /* <-- Number of faces */
 const cs_int_t   *isym,      /* <-- Symmetry indicator:
                                     1: symmetric; 2: not symmetric */
 const cs_int_t   *iresds,    /* <-- Descent smoother type:
                                     0: pcg; 1: Jacobi; 2: cg-stab */
 const cs_int_t   *iresas,    /* <-- Ascent smoother type:
                                     0: pcg; 1: Jacobi; 2: cg-stab */
 const cs_int_t   *ireslp,    /* <-- Coarse Resolution type:
                                     0: pcg; 1: Jacobi; 2: cg-stab */
 const cs_int_t   *ipol,      /* <-- Preconditioning polynomial degree
                                     (0: diagonal) */
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
 const cs_int_t   *ifacel,    /* <-- Face -> cell connectivity  */
 const cs_real_t  *rhs,       /* <-- System right-hand side */
 cs_real_t        *vx         /* <-> System solution */
)
{
  char *var_name;
  cs_sles_type_t type[4] = {CS_SLES_PCG,
                            CS_SLES_JACOBI,
                            CS_SLES_BICGSTAB,
                            CS_SLES_N_TYPES};

  int _iresds = *iresds;
  int _iresas = *iresas;
  int _ireslp = *ireslp;

  cs_bool_t symmetric = (*isym == 1) ? true : false;
  cs_perio_rota_t rotation_mode = CS_PERIO_ROTA_COPY;

  assert(*ncelet >= *ncel);
  assert(*nfac > 0);
  assert(ifacel != NULL);

  if (*iinvpe == 2)
    rotation_mode = CS_PERIO_ROTA_RESET;
  else if (*iinvpe == 3)
    rotation_mode = CS_PERIO_ROTA_IGNORE;

  var_name = cs_base_string_f_to_c_create(cname, *lname);

  assert(*iresds > -1 && *iresds < 3);
  assert(*iresas > -1 && *iresas < 3);
  assert(*ireslp > -1 && *ireslp < 3);

  if (_iresds < 0 || _iresds > 2)
    _iresds = 3;
  if (_iresas < 0 || _iresas > 2)
    _iresas = 3;
  if (_ireslp < 0 || _ireslp > 2)
    _ireslp = 3;

  _multigrid_solve(var_name,
                   type[_iresds],
                   type[_iresas],
                   type[_ireslp],
                   symmetric,
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
 *   post_cell_max <-- If > 0, activates postprocessing of coarseninsg,
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
    _multigrid_info_dump(&((cs_glob_multigrid_systems[ii])->info));

  /* Free multigrid structures */

  for (ii = 0; ii < cs_glob_multigrid_n_systems; ii++)
    _multigrid_destroy(cs_glob_multigrid_systems + ii);

  BFT_FREE(cs_glob_multigrid_systems);

  cs_glob_multigrid_n_systems = 0;
  cs_glob_multigrid_n_max_systems = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
