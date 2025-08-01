/*============================================================================
 * Multigrid solver.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_base_accel.h"
#include "alge/cs_blas.h"
#include "base/cs_dispatch.h"
#include "base/cs_file.h"
#include "alge/cs_grid.h"
#include "base/cs_halo.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "alge/cs_matrix.h"
#include "alge/cs_matrix_default.h"
#include "alge/cs_matrix_spmv_cuda.h"
#include "alge/cs_matrix_util.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "alge/cs_multigrid_smoother.h"
#include "base/cs_post.h"
#include "base/cs_profiling.h"
#include "alge/cs_sles.h"
#include "alge/cs_sles_it.h"
#include "alge/cs_sles_pc.h"
#include "base/cs_timer.h"
#include "base/cs_time_plot.h"
#include "base/cs_time_step.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_multigrid.h"

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

/* SIMD unit size to ensure SIMD alignement (2 to 4 required on most
 * current architectures, so 16 should be enough on most architectures
 * through at least 2012) */

#define CS_SIMD_SIZE(s) (((s-1)/16+1)*16)

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Multigrid types
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_MULTIGRID_MAIN,           /*!< Main multigrid type */
  CS_MULTIGRID_COARSE,         /*!< Coarse level solver */
  CS_MULTIGRID_BOTTOM,         /*!< Bottom level solver */
  CS_MULTIGRID_BOTTOM_SMOOTHE  /*!< Bottom level smoother */

} cs_multigrid_subtype_t;

/*----------------------------------------------------------------------------
 * Local Structure Definitions
 *----------------------------------------------------------------------------*/

/* Basic linear system solver or smoother context and functions */
/*--------------------------------------------------------------*/

typedef struct _cs_mg_sles_t {

  void                     *context;       /* solver context
                                              (options, state, logging) */

  cs_sles_setup_t          *setup_func;    /* solver setup function */
  cs_sles_solve_t          *solve_func;    /* solve function */
  cs_sles_destroy_t        *destroy_func;  /* destruction function */

} cs_mg_sles_t;

/* Basic per linear system options and logging */
/*---------------------------------------------*/

typedef struct _cs_multigrid_info_t {

  /* Settings */

  cs_sles_it_type_t    type[6];             /* Descent/ascent smoothers,
                                               Coarse solver, on host
                                               than device */

  bool                 is_pc;               /* True if used as preconditioner */
  int                  n_max_cycles;        /* Maximum allowed cycles */

  int                  n_max_iter[6];       /* maximum iterations allowed
                                               (host descent/ascent/coarse,
                                               device descent/ascent/coarse) */
  int                  poly_degree[6];      /* polynomial preconditioning degree
                                               (host descent/ascent/coarse,
                                               device descent/ascent/coarse) */

  double               precision_mult[3];   /* solver precision multiplier
                                               (descent/ascent/coarse) */

  /* Logging */

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
  unsigned long long   n_g_rows[4];        /* Global number of rows
                                               (last, min, max, total) */
  unsigned long long   n_elts[3][4];        /* Mean number of rows,
                                               rows + ghosts, and entries
                                               across ranks (last, min, max,
                                               total) */
  double               imbalance[3][4];     /* Imbalance for rows, rows
                                               + ghosts, and entries
                                               (last, min, max, total) */

  unsigned long long   n_it_solve[4];       /* Number of iterations for
                                               solving [last, min, max, total] */
  unsigned long long   n_it_ds_smoothe[4];  /* Number of iterations for
                                               descent smoothing:
                                                 [last, min, max, total] */
  unsigned long long   n_it_as_smoothe[4];  /* Number of iterations for
                                               ascent smoothing:
                                                 [last, min, max, total] */

  unsigned             n_calls[7];          /* Total number of calls:
                                               build, solve, descent smoothe,
                                               ascent smoothe, restrict from
                                               finer, prolong to finer,
                                               BLAS */

  cs_timer_counter_t   t_tot[7];            /* Total timers count:
                                               [build, solve, descent smoothe,
                                               ascent smoothe, restrict from
                                               finer, prolong to finer,
                                               BLAS] */

} cs_multigrid_level_info_t;

/* Grid hierarchy */
/*----------------*/

typedef struct _cs_multigrid_setup_data_t {

  /* Setup */

  unsigned        n_levels;           /* Current number of grid levels */
  unsigned        n_levels_alloc;     /* Allocated number of grid levels */

  cs_grid_t     **grid_hierarchy;     /* Array of grid pointers */
  cs_mg_sles_t   *sles_hierarchy;     /* Contexts for  associated smoothers
                                         and solvers (i*2: descent, coarse;
                                         i*2+1: ascent) */

  /* Level names */

  const char   **lv_names;            /* Names for logging */
  char          *lv_names_buffer;     /* Buffer for level names */

  /* Arrays used only for solving, but maintained until free,
     so as to be usable by convergence error handler. */

  double         exit_initial_residual;  /* Last level initial residual */
  double         exit_residual;          /* Last residual */
  int            exit_level;             /* Last level during solve */
  int            exit_cycle_id;          /* Last cycle id during solve */

  cs_real_t     *rhs_vx_buf;             /* Coarse grid "right hand sides"
                                            and corrections buffer */
  cs_real_t    **rhs_vx;                 /* Coarse grid "right hand sides"
                                            and corrections */

  /* Options used only when used as a preconditioner */

  char          *pc_name;                /* name of preconditioning system */

  int            pc_verbosity;           /* preconditioner verbosity */

  cs_real_t     *pc_aux;                 /* preconditioner auxiliary array */

} cs_multigrid_setup_data_t;

/* Grid hierarchy */
/*----------------*/

struct _cs_multigrid_t {

  /* Settings */

  cs_multigrid_type_t     type;     /* Multigrid type */
  cs_multigrid_subtype_t  subtype;  /* Multigrid subtype */

  int                     f_settings_threshold;  /* Maximum level considered
                                                    as a fine grid for specific
                                                    settings */
  int                     aggregation_limit[2];  /* Maximum allowed fine rows
                                                    per coarse cell (for fine
                                                    and coarse grids) */
  cs_grid_coarsening_t    coarsening_type[2];    /* Coarsening traversal type
                                                   (for fine and coarse grids) */

  int        n_levels_max;       /* Maximum number of grid levels */
  cs_gnum_t  n_g_rows_min;       /* Global number of rows on coarse grids
                                    under which no coarsening occurs */

  int        post_row_max;       /* If > 0, activates postprocessing of
                                    coarsening, projecting coarse cell
                                    numbers (modulo post_row_max)
                                    on the base grid */

  double     p0p1_relax;         /* p0/p1 relaxation_parameter */
  double     k_cycle_threshold;  /* threshold for k cycle */

  /* Setting for use as a preconditioner */

  double     pc_precision;       /* preconditioner precision */
  double     pc_r_norm;          /* preconditioner residual normalization */

  /* Data for postprocessing callback */

  int        post_location;      /* associated mesh location */
  int        n_levels_post;      /* Current number of postprocessed levels */

  int      **post_row_num;       /* If post_row_max > 0, array of
                                    (n_levels - 1) arrays of projected
                                    coarse cell numbers on the base grid */

  int      **post_row_rank;      /* If post_row_max > 0 and grid merging
                                    is active, array of (n_levels - 1) arrays
                                    of projected coarse cell ranks on the
                                    base grid */
  char      *post_name;          /* Name for postprocessing */

  /* Options and maintained state (statistics) */

  cs_multigrid_level_info_t  *lv_info;      /* Info for each level */
  cs_multigrid_t             *lv_mg[3];     /* Optional recursive multigrid
                                               descent, ascent, or null */
  cs_multigrid_t             *p_mg;         /* Optional parent multigrid,
                                               or null */

  cs_multigrid_info_t         info;         /* Base multigrid info */

# if defined(HAVE_MPI)

  /* Communicator used for reduction operations
     (if left at nullptr, main communicator will be used) */

  MPI_Comm comm;
  MPI_Comm caller_comm;

# endif

  int      caller_n_ranks;

  /* Coarse grid rank merging options */

  int      merge_stride;
  int      merge_bottom_n_max_ranks;
  float    merge_bottom_max_row_factor;

  cs_gnum_t merge_mean_threshold;
  cs_gnum_t merge_glob_threshold;

  /* Data available between "setup" and "solve" states */

  cs_multigrid_setup_data_t  *setup_data;   /* setup data */

  cs_time_plot_t             *cycle_plot;       /* plotting of cycles */
  int                         plot_time_stamp;  /* plotting time stamp;
                                                   if < 0, use wall clock */
};

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Multigrid type names */

const char *cs_multigrid_type_name[]
  = {N_("V-cycle"),
     N_("K-cycle"),
     N_("K-cycle, HPC variant"),
     N_("K-cycle, HPC coarse level")};

/* Tunable settings */

static int _k_cycle_hpc_merge_stride = 512;
static int _k_cycle_hpc_recurse_threshold = 256; /* under this size, coarsest
                                                    level solver does not
                                                    use k-cycle preconditioning */

static unsigned int _grid_max_level_for_device = 25; /* grids over this level are
                                                        solved on host only */

/* Force GPU settings (for CPU/GPU comparisons) */

constexpr bool _force_gpu_settings = false;

/*============================================================================
 * Private function prototypes for recursive
 *============================================================================*/

static void
_setup_k_cycle_hpc_sub(cs_multigrid_t     *mg,
                       const char         *name,
                       const cs_matrix_t  *a,
                       int                 n_ranks,
                       int                 verbosity);

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
_multigrid_info_init(cs_multigrid_info_t *info)
{
  int i;

  /* Options */

  info->type[0] = CS_SLES_PCG;
  info->type[1] = CS_SLES_PCG;
  info->type[2] = CS_SLES_PCG;
  info->type[3] = CS_SLES_N_SMOOTHER_TYPES;
  info->type[4] = CS_SLES_N_SMOOTHER_TYPES;
  info->type[5] = CS_SLES_N_SMOOTHER_TYPES;

  info->is_pc        = false;
  info->n_max_cycles = 100;

  info->n_max_iter[0] = 2;
  info->n_max_iter[1] = 10;
  info->n_max_iter[2] = 10000;
  info->n_max_iter[3] = 0;
  info->n_max_iter[4] = 0;
  info->n_max_iter[5] = 0;

  if (cs_get_device_id() > -1 || _force_gpu_settings) {
    info->n_max_iter[3] = 2;
    info->n_max_iter[4] = 10;
    info->n_max_iter[5] = 10000;
  }

  info->poly_degree[0] = -1;
  info->poly_degree[1] = -1;
  info->poly_degree[2] = -1;
  info->poly_degree[3] = -1;
  info->poly_degree[4] = -1;
  info->poly_degree[5] = -1;

  /* In theory, one should increase precision on coarsest mesh,
     but in practice, it is more efficient to have a lower precision,
     so we choose coarse_precision = global_precision; */

  info->precision_mult[0] = -1.;
  info->precision_mult[1] = -1.;
  info->precision_mult[2] = 1.;

  /* Counting and timing */

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
    info->imbalance[i][0] = HUGE_VALF;
    info->imbalance[i][1] = 0.;
  }

  for (i = 0; i < 7; i++)
    CS_TIMER_COUNTER_INIT(info->t_tot[i]);
}

/*----------------------------------------------------------------------------
 * Output information regarding multigrid options.
 *
 * parameters:
 *   mg <-> pointer to multigrid structure
 *----------------------------------------------------------------------------*/

static void
_multigrid_setup_log(const cs_multigrid_t *mg)
{
  if (mg->info.is_pc == false)
    cs_log_printf(CS_LOG_SETUP,
                  _("  Solver type:                       multigrid\n"));
  else
    cs_log_printf(CS_LOG_SETUP,
                  _("  Multigrid preconditioner parameters:\n"));

  if (mg->f_settings_threshold < 1) {
    cs_log_printf(CS_LOG_SETUP,
                  _("  Coarsening type:                   %s\n"
                    "    Max fine rows per coarse row:    %d\n"
                    "    Maximum number of levels :       %d\n"
                    "    Minimum number of coarse rows:   %llu\n"
                    "    P0/P1 relaxation parameter:      %g\n"
                    "  Maximum number of cycles:          %d\n"),
                  _(cs_grid_coarsening_type_name[mg->coarsening_type[1]]),
                  mg->aggregation_limit[1],
                  mg->n_levels_max, (unsigned long long)(mg->n_g_rows_min),
                  mg->p0p1_relax, mg->info.n_max_cycles);
  }
  else {
    cs_log_printf(CS_LOG_SETUP,
                  _("  Coarsening:\n"
                    "    Levels 0 - %d:\n"
                    "      Coarsening type:               %s\n"
                    "      Max fine rows per coarse row:  %d\n"
                    "    Levels %d and coarser:\n"
                    "      Coarsening type:               %s\n"
                    "      Max fine rows per coarse row:  %d\n"
                    "    Maximum number of levels :       %d\n"
                    "    Minimum number of coarse rows:   %llu\n"
                    "    P0/P1 relaxation parameter:      %g\n"
                    "  Maximum number of cycles:          %d\n"),
                  mg->f_settings_threshold-1,
                  _(cs_grid_coarsening_type_name[mg->coarsening_type[0]]),
                  mg->aggregation_limit[0],
                  mg->f_settings_threshold,
                  _(cs_grid_coarsening_type_name[mg->coarsening_type[1]]),
                  mg->aggregation_limit[1],
                  mg->n_levels_max, (unsigned long long)(mg->n_g_rows_min),
                  mg->p0p1_relax, mg->info.n_max_cycles);
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    "  Rank merge parameters:\n"
                    "    merge rank stride:               %d\n"
                    "    mean  coarse rows threshold:    %d\n"
                    "    total coarse rows threshold:    %llu\n"),
                  mg->merge_stride,
                  (int)(mg->merge_mean_threshold),
                  (unsigned long long)(mg->merge_glob_threshold));

    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    "  Bottom rank merge parameters:\n"
                    "    max ranks:                      %d\n"
                    "    max row factor                  %f\n\n"),
                  mg->merge_bottom_n_max_ranks,
                  (double)(mg->merge_bottom_max_row_factor));
  }
#endif

  cs_log_printf(CS_LOG_SETUP,
                _("  Cycle type:                        %s\n"),
                _(cs_multigrid_type_name[mg->type]));

  const char *stage_name[] = {"Descent smoother",
                              "Ascent smoother",
                              "Coarsest level solver"};

  for (int j = 0; j < 3; j++) {
    cs_log_printf(CS_LOG_SETUP, _("  %s:\n"), stage_name[j]);

    for (int k = 1; k > -1; k--) {
      int i = j + k*3;

      if (   mg->info.type[i] != CS_SLES_N_IT_TYPES
          && mg->info.type[i] < CS_SLES_N_SMOOTHER_TYPES) {
        if (k == 1)
          cs_log_printf(CS_LOG_SETUP,
                        _("    Type (device):                   %s\n"),
                        _(cs_sles_it_type_name[mg->info.type[i]]));
        else
          cs_log_printf(CS_LOG_SETUP,
                        _("    Type:                            %s\n"),
                        _(cs_sles_it_type_name[mg->info.type[i]]));

        if (mg->info.poly_degree[i] > -1) {
          cs_log_printf(CS_LOG_SETUP,
                        _("      Preconditioning:               "));
          if (mg->info.poly_degree[i] == 0)
            cs_log_printf(CS_LOG_SETUP, _("Jacobi\n"));
          else if (mg->info.poly_degree[i] < 0) {
            if (mg->lv_mg[i] != nullptr) {
              cs_log_printf(CS_LOG_SETUP, "%s\n",
                            _(cs_multigrid_type_name[mg->lv_mg[i]->type]));
            }
            else
              cs_log_printf(CS_LOG_SETUP, _("None\n"));
          }
          else
            cs_log_printf(CS_LOG_SETUP, _("polynomial, degree %d\n"),
                          mg->info.poly_degree[i]);
        }
        cs_log_printf(CS_LOG_SETUP,
                      _("      Maximum number of iterations:  %d\n"
                        "      Precision multiplier:          %g\n"),
                      mg->info.n_max_iter[i],
                      mg->info.precision_mult[i%3]);
      }
      else if (mg->lv_mg[j] != nullptr) {
        cs_log_printf(CS_LOG_SETUP, "  %s:\n", _(stage_name[j]));
        _multigrid_setup_log(mg->lv_mg[j]);
      }

    }
  }

  cs_log_printf(CS_LOG_SETUP,
                _("  Postprocess coarsening:            %d\n"),
                mg->post_row_max);
}

/*----------------------------------------------------------------------------
 * Output information regarding multigrid resolution.
 *
 * parameters:
 *   mg <-> pointer to multigrid structure
 *----------------------------------------------------------------------------*/

static void
_multigrid_performance_log(const cs_multigrid_t *mg)
{
  unsigned long long n_builds_denom = cs::max(mg->info.n_calls[0], 1u);
  unsigned long long n_solves_denom = cs::max(mg->info.n_calls[1], 1u);
  int n_lv_min = mg->info.n_levels[1];
  int n_lv_max = mg->info.n_levels[2];
  int n_lv_mean = (int)(mg->info.n_levels_tot / n_builds_denom);
  int n_cy_mean = (int)(mg->info.n_cycles[2] / n_solves_denom);

  char tmp_s[7][64] =  {"", "", "", "", "", "", ""};
  const char *stage_name[2] = {N_("Construction:"), N_("Resolution:")};
  const char *lv_stage_name[7] = {N_("build:"), N_("solve:"),
                                  N_("descent smoothe:"), N_("ascent smoothe:"),
                                  N_("restrict:"), N_("prolong:"),
                                  N_("BLAS")};

  if (mg->f_settings_threshold < 1) {
    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("\n"
         "  Multigrid:\n"
         "    %s\n"
         "    Coarsening                     : %s, r %d\n"),
       _(cs_multigrid_type_name[mg->type]),
       _(cs_grid_coarsening_type_name[mg->coarsening_type[1]]),
       mg->aggregation_limit[1]);
  }
  else {
    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("\n"
         "  Multigrid:\n"
         "    %s\n"
         "    Coarsening (0 - %d)             : %s, r %d\n"
         "    Coarsening (%d +)               : %s, r %d\n"),
       _(cs_multigrid_type_name[mg->type]),
       mg->f_settings_threshold-1,
       _(cs_grid_coarsening_type_name[mg->coarsening_type[0]]),
       mg->aggregation_limit[0],
       mg->f_settings_threshold,
       _(cs_grid_coarsening_type_name[mg->coarsening_type[1]]),
       mg->aggregation_limit[1]);
  }

  const char *stage_type_name[] =
    {"Descent smoother              ",
     "Ascent smoother               ",
     "Coarsest level solver         ",
     "Descent smoother (device)     ",
     "Ascent smoother (device)      ",
     "Coarsest level solver (device)"};

  for (int j = 0; j < 6; j++) {
    int i = j/2 + 3*((j+1)%2);

    if (   mg->info.type[i] != CS_SLES_N_IT_TYPES
        && mg->info.type[i] < CS_SLES_N_SMOOTHER_TYPES) {
      cs_log_printf(CS_LOG_PERFORMANCE,
                  _("    %s : %s\n"),
                  _(stage_type_name[i]),
                    _(cs_sles_it_type_name[mg->info.type[i]]));
    }
  }

  sprintf(tmp_s[0], "%-36s", "");
  cs_log_strpadl(tmp_s[1], _(" mean"), 12, 64);
  cs_log_strpadl(tmp_s[2], _("minimum"), 12, 64);
  cs_log_strpadl(tmp_s[3], _("maximum"), 12, 64);

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n  %s %s %s %s\n",
                tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

  cs_log_strpad(tmp_s[0], _("Number of levels:"), 36, 64);
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

  for (unsigned i = 0; i < mg->info.n_levels[2]; i++) {

    const cs_multigrid_level_info_t *lv_info = mg->lv_info + i;
    unsigned long long n_lv_builds = lv_info->n_calls[0];

    if (n_lv_builds < 1)
      continue;

    cs_log_strpad(tmp_s[0], _("Number of rows:"), 34, 64);
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  Grid level %d:\n"
                    "    %s %12llu %12llu %12llu\n"),
                  i, tmp_s[0],
                  lv_info->n_g_rows[3] / n_lv_builds,
                  lv_info->n_g_rows[1], lv_info->n_g_rows[2]);

    if (mg->caller_n_ranks == 1) {
      cs_log_strpad(tmp_s[1], _("Number of entries:"), 34, 64);
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "    %s %12llu %12llu %12llu\n",
                    tmp_s[1],
                    lv_info->n_elts[2][3] / n_lv_builds,
                    lv_info->n_elts[2][1], lv_info->n_elts[2][2]);
    }

#if defined(HAVE_MPI)

    if (mg->caller_n_ranks > 1) {
      cs_log_strpad(tmp_s[0], _("Number of active ranks:"), 34, 64);
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "    %s %12llu %12llu %12llu\n",
                    tmp_s[0],
                    lv_info->n_ranks[3] / n_lv_builds,
                    lv_info->n_ranks[1], lv_info->n_ranks[2]);
      cs_log_strpad(tmp_s[0], _("Mean local rows:"), 34, 64);
      cs_log_strpad(tmp_s[1], _("Mean local columns + ghosts:"), 34, 64);
      cs_log_strpad(tmp_s[2], _("Mean local entries:"), 34, 64);
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
      if (lv_info->n_ranks[2] > 1) {
        cs_log_strpad(tmp_s[0], _("Rows imbalance:"), 34, 64);
        cs_log_strpad(tmp_s[1], _("Columns + ghosts imbalance:"), 34, 64);
        cs_log_strpad(tmp_s[2], _("entries imbalance"), 34, 64);
        cs_log_printf(CS_LOG_PERFORMANCE,
                      "    %-34s %12.3f %12.3f %12.3f\n"
                      "    %-34s %12.3f %12.3f %12.3f\n"
                      "    %-34s %12.3f %12.3f %12.3f\n",
                      tmp_s[0],
                      lv_info->imbalance[0][3] / n_lv_builds,
                      lv_info->imbalance[0][1], lv_info->imbalance[0][2],
                      tmp_s[1],
                      lv_info->imbalance[1][3] / n_lv_builds,
                      lv_info->imbalance[1][1], lv_info->imbalance[1][2],
                      tmp_s[2],
                      lv_info->imbalance[2][3] / n_lv_builds,
                      lv_info->imbalance[2][1], lv_info->imbalance[2][2]);
      }
    }

#endif

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

  for (unsigned i = 0; i < mg->info.n_levels[2]; i++) {

    const cs_multigrid_level_info_t *lv_info = mg->lv_info + i;

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  Grid level %d:\n"), i);

    cs_log_timer_array(CS_LOG_PERFORMANCE,
                       4,                  /* indent, */
                       7,                  /* n_lines */
                       lv_stage_name,
                       lv_info->n_calls,
                       lv_info->t_tot);

  }

  {
    const char *names[]
      = {N_("coarse level descent smoother"),
         N_("coarse level ascent smoother"),
         N_("bottom level solver")};

    for (int i = 0; i < 3; i++) {
      if (mg->lv_mg[i] != nullptr) {

        cs_log_printf(CS_LOG_PERFORMANCE,
                      _("\n"
                        "  Nested %s:\n"), _(names[i]));
        cs_multigrid_log(mg->lv_mg[i], CS_LOG_PERFORMANCE);
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Create empty structure used to maintain setup data
 * (between cs_sles_setup and cs_sles_free type calls.
 *
 * returns:
 *   pointer to multigrid setup data structure
 *----------------------------------------------------------------------------*/

static cs_multigrid_setup_data_t *
_multigrid_setup_data_create(void)
{
  cs_multigrid_setup_data_t *mgd;

  CS_MALLOC(mgd, 1, cs_multigrid_setup_data_t);

  mgd->n_levels = 0;
  mgd->n_levels_alloc = 0;

  mgd->grid_hierarchy = nullptr;
  mgd->sles_hierarchy = nullptr;

  mgd->lv_names = nullptr;
  mgd->lv_names_buffer = nullptr;

  mgd->exit_initial_residual = -1.;
  mgd->exit_residual = -1.;
  mgd->exit_level = -1.;
  mgd->exit_cycle_id = -1.;

  mgd->rhs_vx_buf = nullptr;
  mgd->rhs_vx = nullptr;

  mgd->pc_name = nullptr;
  mgd->pc_aux = nullptr;
  mgd->pc_verbosity = 0;

  return mgd;
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
  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  unsigned ii;

  /* Reallocate arrays if necessary */

  if (mgd->n_levels == mgd->n_levels_alloc) {

    /* Max previous */
    unsigned int n_lv_max_prev = cs::max(mg->info.n_levels[2],
                                         mgd->n_levels);

    if (mgd->n_levels_alloc == 0) {
      mgd->n_levels_alloc = n_lv_max_prev;
      if (mgd->n_levels_alloc == 0)
        mgd->n_levels_alloc = 10;
    }
    else
      mgd->n_levels_alloc *= 2;

    CS_REALLOC(mgd->grid_hierarchy, mgd->n_levels_alloc, cs_grid_t *);
    CS_REALLOC(mgd->sles_hierarchy, mgd->n_levels_alloc*2, cs_mg_sles_t);

    for (ii = mgd->n_levels; ii < mgd->n_levels_alloc*2; ii++) {
      mgd->sles_hierarchy[ii].context = nullptr;
      mgd->sles_hierarchy[ii].setup_func = nullptr;
      mgd->sles_hierarchy[ii].solve_func = nullptr;
      mgd->sles_hierarchy[ii].destroy_func = nullptr;
    }

    if (n_lv_max_prev < mgd->n_levels_alloc) {
      CS_REALLOC(mg->lv_info, mgd->n_levels_alloc, cs_multigrid_level_info_t);
      for (ii = n_lv_max_prev; ii < mgd->n_levels_alloc; ii++)
        _multigrid_level_info_init(mg->lv_info + ii);
    }

  }

  /* Add new grid to hierarchy */

  mgd->grid_hierarchy[mgd->n_levels] = grid;
  mg->info.n_levels[0] = mgd->n_levels + 1;

  if (mg->post_row_num != nullptr) {
    int n_max_post_levels = (int)(mg->info.n_levels[2]) - 1;
    CS_REALLOC(mg->post_row_num, mgd->n_levels_alloc, int *);
    for (ii = n_max_post_levels + 1; ii < mgd->n_levels_alloc; ii++)
      mg->post_row_num[ii] = nullptr;
  }

  if (mg->post_row_rank != nullptr) {
    int n_max_post_levels = (int)(mg->info.n_levels[2]) - 1;
    CS_REALLOC(mg->post_row_rank, mgd->n_levels_alloc, int *);
    for (ii = n_max_post_levels + 1; ii < mgd->n_levels_alloc; ii++)
      mg->post_row_rank[ii] = nullptr;
  }

  /* Ready for next level */

  mgd->n_levels += 1;
}

/*----------------------------------------------------------------------------
 * Update information after adding grids to multigrid structure hierarchy.
 *
 * parameters:
 *   mg   <-- multigrid structure
 *----------------------------------------------------------------------------*/

static void
_update_level_info(cs_multigrid_t  *mg)
{
  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  for (unsigned lv = 0; lv < mgd->n_levels; lv++) {
    const cs_grid_t *grid = mgd->grid_hierarchy[lv];

    int  n_ranks;
    cs_lnum_t  n_rows, n_rows_with_ghosts, n_entries;
    cs_gnum_t  n_g_rows;
    cs_multigrid_level_info_t  *lv_info = mg->lv_info + lv;

    cs_grid_get_info(grid,
                     nullptr,
                     nullptr,
                     nullptr,
                     nullptr,
                     &n_ranks,
                     &n_rows,
                     &n_rows_with_ghosts,
                     &n_entries,
                     &n_g_rows);

    lv_info->n_ranks[0] = n_ranks;
    if (lv_info->n_ranks[1] > (unsigned)n_ranks)
      lv_info->n_ranks[1] = n_ranks;
    else if (lv_info->n_ranks[2] < (unsigned)n_ranks)
      lv_info->n_ranks[2] = n_ranks;
    lv_info->n_ranks[3] += n_ranks;

    lv_info->n_g_rows[0] = n_g_rows;
    if (lv_info->n_g_rows[1] > n_g_rows)
      lv_info->n_g_rows[1] = n_g_rows;
    else if (lv_info->n_g_rows[2] < n_g_rows)
      lv_info->n_g_rows[2] = n_g_rows;
    lv_info->n_g_rows[3] += n_g_rows;

    lv_info->n_elts[0][0] = n_rows;
    lv_info->n_elts[1][0] = n_rows_with_ghosts;
    lv_info->n_elts[2][0] = n_entries;

    for (unsigned ii = 0; ii < 3; ii++) {
      if (lv_info->n_elts[ii][1] > lv_info->n_elts[ii][0])
        lv_info->n_elts[ii][1] = lv_info->n_elts[ii][0];
      else if (lv_info->n_elts[ii][2] < lv_info->n_elts[ii][0])
        lv_info->n_elts[ii][2] = lv_info->n_elts[ii][0];
      lv_info->n_elts[ii][3] += lv_info->n_elts[ii][0];
    }

#if defined(HAVE_MPI)

    if (mg->caller_n_ranks > 1) {
      cs_gnum_t tot_sizes[3], max_sizes[3];
      cs_gnum_t loc_sizes[3] = {
        (cs_gnum_t)n_rows, (cs_gnum_t)n_rows_with_ghosts, (cs_gnum_t)n_entries};
      MPI_Allreduce(loc_sizes, tot_sizes, 3, CS_MPI_GNUM, MPI_SUM,
                    mg->caller_comm);
      MPI_Allreduce(loc_sizes, max_sizes, 3, CS_MPI_GNUM, MPI_MAX,
                    mg->caller_comm);
      for (unsigned ii = 0; ii < 3; ii++) {
        if (tot_sizes[ii] > 0)
          lv_info->imbalance[ii][0] = (  max_sizes[ii]
                                       / (tot_sizes[ii]*1.0/n_ranks)) - 1.0;
        else
          lv_info->imbalance[ii][0] = 0;

        if (lv_info->imbalance[ii][1] > lv_info->imbalance[ii][0])
          lv_info->imbalance[ii][1] = lv_info->imbalance[ii][0];
        else if (lv_info->imbalance[ii][2] < lv_info->imbalance[ii][0])
          lv_info->imbalance[ii][2] = lv_info->imbalance[ii][0];
        lv_info->imbalance[ii][3] += lv_info->imbalance[ii][0];
      }
    }

#endif /* defined(HAVE_MPI) */

    if (lv_info->n_calls[0] == 0) {
      lv_info->n_ranks[1] = n_ranks;
      lv_info->n_g_rows[1] = n_g_rows;
      for (unsigned ii = 0; ii < 3; ii++) {
        lv_info->n_elts[ii][1] = lv_info->n_elts[ii][0];
#if defined(HAVE_MPI)
        lv_info->imbalance[ii][1] = lv_info->imbalance[ii][0];
#endif
      }
    }

    lv_info->n_calls[0] += 1;
  }
}

/*----------------------------------------------------------------------------
 * Add postprocessing info to multigrid hierarchy
 *
 * parameters:
 *   mg            <-> multigrid structure
 *   name          <-- postprocessing name
 *   post_location <-- where to postprocess (cells or vertices)
 *   n_base_cells  <-- number of cells in base grid
 *----------------------------------------------------------------------------*/

static void
_multigrid_add_post(cs_multigrid_t  *mg,
                    const char      *name,
                    int              post_location,
                    cs_lnum_t        n_base_rows)
{
  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  int ii;

  assert(mg != nullptr);

  if (mg->post_row_max < 1)
    return;

  mg->post_location = post_location;
  mg->n_levels_post = mgd->n_levels - 1;

  CS_REALLOC(mg->post_name, strlen(name) + 1, char);
  strcpy(mg->post_name, name);

  assert(mg->n_levels_post <= mg->n_levels_max);

  /* Reallocate arrays if necessary */

  if (mg->post_row_num == nullptr) {
    CS_MALLOC(mg->post_row_num, mg->n_levels_max, int *);
    for (ii = 0; ii < mg->n_levels_max; ii++)
      mg->post_row_num[ii] = nullptr;
  }

  if (mg->post_row_rank == nullptr && mg->merge_stride > 1) {
    CS_MALLOC(mg->post_row_rank, mg->n_levels_max, int *);
    for (ii = 0; ii < mg->n_levels_max; ii++)
      mg->post_row_rank[ii] = nullptr;
  }

  for (ii = 0; ii < mg->n_levels_post; ii++) {
    CS_REALLOC(mg->post_row_num[ii], n_base_rows, int);
    cs_grid_project_row_num(mgd->grid_hierarchy[ii+1],
                            n_base_rows,
                            mg->post_row_max,
                            mg->post_row_num[ii]);
  }

  if (mg->post_row_rank != nullptr) {
    for (ii = 0; ii < mg->n_levels_post; ii++) {
      CS_REALLOC(mg->post_row_rank[ii], n_base_rows, int);
      cs_grid_project_row_rank(mgd->grid_hierarchy[ii+1],
                               n_base_rows,
                               mg->post_row_rank[ii]);
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
  CS_UNUSED(ts);

  int ii;
  size_t name_len;
  char *var_name = nullptr;
  cs_multigrid_t *mg = (cs_multigrid_t *)mgh;
  const char *base_name = nullptr;

  /* Return if necessary structures inconsistent or have been destroyed */

  if (mg == nullptr)
    return;

  if (mg->post_row_num == nullptr || cs_post_mesh_exists(-1) != true)
    return;

  int *s_num = nullptr;
  const cs_range_set_t *rs = nullptr;
  if (mg->post_location == CS_MESH_LOCATION_VERTICES) {
    CS_MALLOC(s_num, cs_glob_mesh->n_vertices, int);
    rs = cs_glob_mesh->vtx_range_set;
  }

  /* Allocate name buffer */

  base_name = mg->post_name;
  name_len = 3 + strlen(base_name) + 1 + 3 + 1 + 4 + 1;
  CS_MALLOC(var_name, name_len, char);

  /* Loop on grid levels */

  for (ii = 0; ii < mg->n_levels_post; ii++) {

    sprintf(var_name, "mg %s %2d", base_name, (ii+1));

    if (mg->post_location == CS_MESH_LOCATION_CELLS)
      cs_post_write_var(CS_POST_MESH_VOLUME,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        var_name,
                        1,
                        false,
                        true,
                        CS_POST_TYPE_int,
                        mg->post_row_num[ii],
                        nullptr,
                        nullptr,
                        cs_glob_time_step);

    else if (mg->post_location == CS_MESH_LOCATION_VERTICES) {
      cs_range_set_scatter(rs, CS_INT_TYPE, 1, mg->post_row_num[ii], s_num);
      cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                               CS_POST_WRITER_ALL_ASSOCIATED,
                               var_name,
                               1,
                               false,
                               true,
                               CS_POST_TYPE_int,
                               s_num,
                               cs_glob_time_step);
    }

    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid location for post-processing.\n", __func__);

    CS_FREE(mg->post_row_num[ii]);

    if (mg->post_row_rank != nullptr) {

      sprintf(var_name, "rk %s %2d",
              base_name, (ii+1));

      if (mg->post_location == CS_MESH_LOCATION_CELLS)
        cs_post_write_var(CS_POST_MESH_VOLUME,
                          CS_POST_WRITER_ALL_ASSOCIATED,
                          var_name,
                          1,
                          false,
                          true,
                          CS_POST_TYPE_int,
                          mg->post_row_rank[ii],
                          nullptr,
                          nullptr,
                          cs_glob_time_step);
      else if (mg->post_location == CS_MESH_LOCATION_VERTICES) {
        cs_range_set_scatter(rs, CS_INT_TYPE, 1, mg->post_row_rank[ii], s_num);
        cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                                 CS_POST_WRITER_ALL_ASSOCIATED,
                                 var_name,
                                 1,
                                 false,
                                 true,
                                 CS_POST_TYPE_int,
                                 s_num,
                                 cs_glob_time_step);
      }

      CS_FREE(mg->post_row_rank[ii]);

    }

  }
  mg->n_levels_post = 0;

  CS_FREE(s_num);
  CS_FREE(var_name);
}

/*----------------------------------------------------------------------------
 * Function returning the type name of the multigrid preconditioner context.
 *
 * parameters:
 *   context   <-- pointer to preconditioner context
 *   logging   <-- if true, logging description; if false, canonical name
 *----------------------------------------------------------------------------*/

static const char *
_multigrid_pc_get_type(const void  *context,
                       bool         logging)
{
  CS_UNUSED(context);

  if (logging == false) {
    static const char t[] = "multigrid";
    return t;
  }
  else {
    static const char t[] = N_("Multigrid");
    return _(t);
  }
}

/*----------------------------------------------------------------------------
 * Function for setup of a multigrid preconditioner.
 *
 * parameters:
 *   context   <-> pointer to preconditioner context
 *   name      <-- pointer to name of associated linear system
 *   a         <-- matrix
 *   accel     <-- use accelerator version ?
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

static void
_multigrid_pc_setup(void               *context,
                    const char         *name,
                    const cs_matrix_t  *a,
                    bool                accel,
                    int                 verbosity)
{
  CS_UNUSED(accel);

  cs_multigrid_setup(context, name, a, verbosity);

  cs_multigrid_t  *mg = (cs_multigrid_t  *)context;
  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  CS_REALLOC(mgd->pc_name, strlen(name) + 1, char);
  strcpy(mgd->pc_name, name);
}

/*----------------------------------------------------------------------------
 * Function for setup of a multigrid preconditioner.
 *
 * parameters:
 *   context   <-> pointer to preconditioner context
 *   name      <-- pointer to name of associated linear system
 *   a         <-- matrix
 *   accel     <-- use accelerator version ?
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

static void
_multigrid_pc_setup_k_sub(void               *context,
                          const char         *name,
                          const cs_matrix_t  *a,
                          bool                accel,
                          int                 verbosity)
{
  CS_UNUSED(accel);

  cs_multigrid_t  *mg = (cs_multigrid_t *)context;
  cs_multigrid_t  *parent = mg->p_mg;
  cs_multigrid_setup_data_t *p_mgd = parent->setup_data;

  int n_ranks = 1;
#if defined(HAVE_MPI)
  {
    int p_lv = p_mgd->n_levels-1;
    MPI_Comm lv_comm = cs_grid_get_comm(p_mgd->grid_hierarchy[p_lv]);
    if (lv_comm != MPI_COMM_NULL)
      MPI_Comm_size(lv_comm, &n_ranks);
  }
#endif

  _setup_k_cycle_hpc_sub(mg, name, a, n_ranks, verbosity);

  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  CS_REALLOC(mgd->pc_name, strlen(name) + 1, char);
  strcpy(mgd->pc_name, name);
}

/*----------------------------------------------------------------------------
 * Function for setup of a multigrid (local smoother)
 *
 * parameters:
 *   context   <-> pointer to preconditioner context
 *   name      <-- pointer to name of associated linear system
 *   a         <-- matrix
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

static void
_multigrid_setup_k_local_smoothe(void               *context,
                                 const char         *name,
                                 const cs_matrix_t  *a,
                                 int                 verbosity)
{
  cs_multigrid_t  *mg = (cs_multigrid_t *)context;

  _setup_k_cycle_hpc_sub(mg, name, a, 1, verbosity);

  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  CS_REALLOC(mgd->pc_name, strlen(name) + 1, char);
  strcpy(mgd->pc_name, name);
}

/*----------------------------------------------------------------------------
 * Function or setting of the required tolerance for multigrid as
 * as a preconditioner.
 *
 * The preconditioner is considered to have converged when
 * residual/r_norm <= precision, residual being the L2 norm of a.vx-rhs.
 *
 * parameters:
 *   context       <-> pointer to multigrid context
 *   precision     <-- preconditioner precision
 *   r_norm        <-- residual normalization
 *----------------------------------------------------------------------------*/

static void
_multigrid_pc_tolerance_t(void    *context,
                          double   precision,
                          double   r_norm)
{
  cs_multigrid_t  *mg = (cs_multigrid_t *)context;

  mg->pc_precision = precision;
  mg->pc_r_norm = r_norm;
}

/*----------------------------------------------------------------------------
 * Function for application of a Multigrid preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to nullptr, and x_out contains the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

static cs_sles_pc_state_t
_multigrid_pc_apply(void                *context,
                    const cs_real_t     *x_in,
                    cs_real_t           *x_out)
{
  CS_PROFILE_FUNC_RANGE();

  int     n_iter;
  double  residual;

  cs_multigrid_t  *mg = (cs_multigrid_t *)context;
  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  const cs_matrix_t  *a = cs_grid_get_matrix(mgd->grid_hierarchy[0]);

  const cs_real_t *rhs = x_in;

  /* Special case where no coarse grid is available
     (may happen with zero extra-diagonal or penalization everywhere):
     no preconditioning (should be an easy solve anyways) */

  if (mgd->n_levels == 1) {

    /* If preconditioner is "in-place", no-op */

    if (x_in != nullptr && x_in != x_out) {
      const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
      const cs_lnum_t n_rows = cs_matrix_get_n_rows(a) * db_size;

      cs_alloc_mode_t amode = cs_matrix_get_alloc_mode(a);
      bool use_gpu = (amode > CS_ALLOC_HOST) ? true : false;

      cs_dispatch_context ctx;
      ctx.set_use_gpu(use_gpu);

      ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
        x_out[ii] = x_in[ii];
      });
      ctx.wait();
    }

    mg->info.n_calls[1] += 1;

    return CS_SLES_PC_CONVERGED;
  }

  /* If preconditioner is "in-place", use additional buffer */

  if (x_in == nullptr) {
    const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
    const cs_lnum_t n_rows = cs_matrix_get_n_rows(a) * db_size;

    cs_alloc_mode_t amode = cs_matrix_get_alloc_mode(a);
    bool use_gpu = (amode > CS_ALLOC_HOST) ? true : false;

    cs_dispatch_context ctx;
    ctx.set_use_gpu(use_gpu);

    if (mgd->pc_aux == nullptr) {
      const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * db_size;
      CS_MALLOC_HD(mgd->pc_aux, n_cols, cs_real_t, amode);
    }
    cs_real_t *restrict _rhs = mgd->pc_aux;

    ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      _rhs[ii] = x_out[ii];
    });
    ctx.wait();

    rhs = _rhs;
  }

  cs_sles_convergence_state_t  cvg = cs_multigrid_solve(context,
                                                        mgd->pc_name,
                                                        a,
                                                        mgd->pc_verbosity,
                                                        mg->pc_precision,
                                                        mg->pc_r_norm,
                                                        &n_iter,
                                                        &residual,
                                                        rhs,
                                                        nullptr,
                                                        x_out,
                                                        0,
                                                        nullptr);

  cs_sles_pc_state_t state;

  switch(cvg) {
  case CS_SLES_DIVERGED:
    state = CS_SLES_PC_DIVERGED;
    break;
  case CS_SLES_BREAKDOWN:
    state = CS_SLES_PC_BREAKDOWN;
    break;
  case CS_SLES_CONVERGED:
    state = CS_SLES_PC_CONVERGED;
    break;
  default:
    state = CS_SLES_PC_MAX_ITERATION;
  }

  return state;
}

/*----------------------------------------------------------------------------
 * Create a Preconditioner structure using multigrid.
 *
 * parameters:
 *   mg_type  <-- type of multigrid algorithm to use
 *
 * returns:
 *   pointer to newly created preconditioner object.
 *----------------------------------------------------------------------------*/

static cs_multigrid_t *
_multigrid_pc_create(cs_multigrid_type_t  mg_type)
{
  cs_multigrid_t *mg = cs_multigrid_create(mg_type);

  switch(mg_type) {
  case CS_MULTIGRID_V_CYCLE:
    cs_multigrid_set_solver_options
      (mg,
       CS_SLES_P_SYM_GAUSS_SEIDEL, /* descent smoothe */
       CS_SLES_P_SYM_GAUSS_SEIDEL, /* ascent smoothe */
       CS_SLES_PCG,                /* coarse smoothe */
       1,                          /* n_max_cycles */
       1,                          /* n_max_iter_descent, */
       1,                          /* n_max_iter_ascent */
       500,                        /* n_max_iter_coarse */
       0, 0, -1,                   /* precond poly_degree */
       -1, -1, 1);                 /* precision_multiplier */
    break;

  case CS_MULTIGRID_K_CYCLE:
    cs_multigrid_set_solver_options
      (mg,
       CS_SLES_TS_F_GAUSS_SEIDEL,
       CS_SLES_TS_B_GAUSS_SEIDEL,
       CS_SLES_PCG,                /* coarse smoothe */
       1,    /* n max cycles */
       1,    /* n max iter for descent */
       1,    /* n max iter for ascent */
       500,  /* n max iter for coarse solve */
       0, 0, 0,     /* precond degree */
       -1, -1, 10); /* precision multiplier */
    break;

  case CS_MULTIGRID_K_CYCLE_HPC:
    cs_multigrid_set_solver_options
      (mg,
       CS_SLES_TS_F_GAUSS_SEIDEL,
       CS_SLES_TS_B_GAUSS_SEIDEL,
       CS_SLES_FCG,                /* coarse smoothe */
       1,   /* n max cycles */
       1,   /* n max iter for descent */
       1,   /* n max iter for ascent */
       500, /* n max iter for coarse solve */
       0, 0, 0,    /* precond degree */
       -1, -1, 1); /* precision multiplier */
    break;

  default:
    assert(0);
  }

  return mg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a multigrid smoother using a local K cycle for coarse-level
 * solvers used by high-performance system variant of K-cycle.
 *
 * \param[in]  parent  parent-level multigrid, or nullptr

 * \return  pointer to newly created mutigrid object.
 */
/*----------------------------------------------------------------------------*/

static cs_multigrid_t *
_multigrid_create_k_cycle_bottom_smoother(cs_multigrid_t  *parent)
{
  cs_multigrid_t *mg = _multigrid_pc_create(CS_MULTIGRID_K_CYCLE);

  mg->subtype = CS_MULTIGRID_BOTTOM_SMOOTHE;
  mg->info.is_pc = true;
  mg->f_settings_threshold = -1;
  mg->aggregation_limit[1] = 8;
  mg->k_cycle_threshold = -1;

#if defined(HAVE_MPI)
  mg->comm = MPI_COMM_NULL;
  if (parent != nullptr)
    mg->caller_comm = parent->comm;
  else
    mg->caller_comm = cs_glob_mpi_comm;
#endif

  return mg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destructor for preconditioner using persistent multigrid context.
 *
 * Using this destructor, the multigrid context is only freed, not destroyed.
 *
 * \param[in, out]  context  pointer to multigrid linear solver info
 *                           (actual type: cs_multigrid_t  **)
 */
/*----------------------------------------------------------------------------*/

static void
_pc_from_mg_destroy(void  **context)
{
  cs_multigrid_t *mg = (cs_multigrid_t *)(*context);

  if (mg != nullptr)
    cs_multigrid_free(mg);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a multigrid preconditioner based onusing a K cycle for coarse-level
 * solvers used by high-performance system variant of K-cycle.
 *
 * \param[in]  mg  associated coarse level multigrid structure
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

static cs_sles_pc_t *
_pc_create_from_mg_sub(cs_multigrid_t  *mg)
{
  cs_sles_pc_t *pc = cs_sles_pc_define(mg,
                                       _multigrid_pc_get_type,
                                       _multigrid_pc_setup_k_sub,
                                       _multigrid_pc_tolerance_t,
                                       _multigrid_pc_apply,
                                       cs_multigrid_free,
                                       cs_multigrid_log,
                                       cs_multigrid_copy,
                                       _pc_from_mg_destroy);

  return pc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create recursive high-performance system variant of K-cycle
 *        for coarse-solver bottom-level.
 *
 * Returns nullptr if the number of associated ranks is small enough that
 * no multigrid solver is required.
 *
 * \param[in]  parent  parent-level multigrid, or nullptr

 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

static cs_multigrid_t *
_multigrid_create_k_cycle_bottom_coarsest(cs_multigrid_t  *parent)
{
  cs_multigrid_t *mg = nullptr;

#if defined(HAVE_MPI)
  if (parent != nullptr) {
    int size = 1;
    if (parent->comm != MPI_COMM_NULL)
      MPI_Comm_size(parent->comm, &size);

    if (size >= _k_cycle_hpc_merge_stride) {
      if (size > _k_cycle_hpc_merge_stride)
        mg = _multigrid_pc_create(CS_MULTIGRID_K_CYCLE_HPC); /* recursive */
      else {
        if (size > _k_cycle_hpc_recurse_threshold) {
          mg = _multigrid_pc_create(CS_MULTIGRID_K_CYCLE);
          mg->f_settings_threshold = -1;
          mg->aggregation_limit[1] = 8;
          mg->k_cycle_threshold = -1;
        }
      }
      if (mg != nullptr) {
        mg->caller_comm = parent->comm;
        mg->comm = cs_grid_get_comm_merge(parent->comm,
                                          _k_cycle_hpc_merge_stride);
        mg->subtype = CS_MULTIGRID_COARSE;
      }
    }

  }
#endif /* defined(HAVE_MPI) */

  return mg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a 2 level coarse-level solver for use
 *        by high-performance system variant of K-cycle.
 *
 * \param[in]  parent  parent-level multigrid, or nullptr

 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

static cs_multigrid_t *
_multigrid_create_k_cycle_bottom(cs_multigrid_t  *parent)
{
  cs_multigrid_t *mg = cs_multigrid_create(CS_MULTIGRID_V_CYCLE);

  mg->subtype = CS_MULTIGRID_BOTTOM;
  mg->p_mg = parent;

  mg->f_settings_threshold = -1;
  mg->aggregation_limit[1] = 0;
  mg->coarsening_type[1] = CS_GRID_COARSENING_DEFAULT; /* Not used here */
  mg->n_levels_max = 2;
  mg->n_g_rows_min = 1;

  mg->pc_precision = 0.0;
  mg->pc_r_norm = 0.0;

  mg->n_levels_post = 0;
  mg->setup_data = nullptr;

  CS_REALLOC(mg->lv_info, mg->n_levels_max, cs_multigrid_level_info_t);

  for (int ii = 0; ii < mg->n_levels_max; ii++)
    _multigrid_level_info_init(mg->lv_info + ii);

  mg->caller_n_ranks = cs_glob_n_ranks;

#if defined(HAVE_MPI)
  mg->merge_stride = _k_cycle_hpc_merge_stride;
  mg->merge_mean_threshold = 0;
  mg->merge_glob_threshold = 1;
#endif

#if defined(HAVE_MPI)
  if (parent != nullptr) {
    mg->comm = parent->comm;
    mg->caller_comm = parent->comm;
  }
#endif /* defined(HAVE_MPI) */

  /* local k-cycle for level 0 smoother,
     k-cycle preconditioning for level 1 solver (such as FCG) */

  mg->lv_mg[0] = _multigrid_create_k_cycle_bottom_smoother(mg);
  mg->lv_mg[2] = _multigrid_create_k_cycle_bottom_coarsest(mg);

  mg->lv_mg[0]->p_mg = parent;
  if (mg->lv_mg[2] != nullptr)
    mg->lv_mg[2]->p_mg = mg;

  int pc_degree_bottom = (mg->lv_mg[2] != nullptr) ? -1 : 1;

  cs_multigrid_set_solver_options
    (mg,
     CS_SLES_N_SMOOTHER_TYPES, CS_SLES_N_SMOOTHER_TYPES, CS_SLES_FCG,
     1,   /* n max cycles */
     1,   /* n max iter for descent */
     1,   /* n max iter for ascent */
     500, /* n max iter coarse */
     -1, -1, pc_degree_bottom,   /* precond degree */
     1, 1, 1);  /* precision multiplier */

  return mg;
}

/*----------------------------------------------------------------------------
 * Allocate working array for coarse right hand sides and corrections.
 *
 * parameters:
 *   mg      <-> pointer to multigrid solver info and context
 *   stride  <-- matrix fill stride
 *----------------------------------------------------------------------------*/

static void
_multigrid_setup_sles_work_arrays(cs_multigrid_t  *mg,
                                  cs_lnum_t        stride)
{
  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  unsigned n0 = 0, n1 = 2;
  if (   mg->type >= CS_MULTIGRID_K_CYCLE
      && mg->type <= CS_MULTIGRID_K_CYCLE_HPC) {
    n0 = 4;
    n1 = n0 + 6;
  }

  CS_MALLOC(mgd->rhs_vx, mgd->n_levels*n1, cs_real_t *);

  for (unsigned i = 0; i < n1; i++)
    mgd->rhs_vx[i] = nullptr;

  if (mgd->n_levels > 1) {

    size_t wr_size0 = 0, wr_size1 = 0;
    for (unsigned i = 0; i < mgd->n_levels; i++) {
      size_t block_size
        = cs_grid_get_n_cols_max(mgd->grid_hierarchy[i])*stride;
      block_size = CS_SIMD_SIZE(block_size);
      if (i == 0)
        wr_size0 += block_size;
      else
        wr_size1 += block_size;
    }

    CS_MALLOC_HD(mgd->rhs_vx_buf, wr_size0*n0 + wr_size1*n1, cs_real_t,
                 cs_alloc_mode);

    size_t block_size_shift = 0;

    for (unsigned i = 0; i < mgd->n_levels; i++) {
      size_t block_size
        = cs_grid_get_n_cols_max(mgd->grid_hierarchy[i])*stride;
      block_size = CS_SIMD_SIZE(block_size);
      cs_lnum_t n = (i == 0) ? n0 : n1;
      for (int j = 0; j < n; j++) {
        mgd->rhs_vx[i*n1 + j] = mgd->rhs_vx_buf+ block_size_shift;
        block_size_shift += block_size;
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Setup multigrid sparse linear equation solvers on existing hierarchy.
 *
 * parameters:
 *   mg        <-> pointer to multigrid solver info and context
 *   name      <-- linear system name
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

static void
_multigrid_setup_sles_k_cycle_bottom(cs_multigrid_t  *mg,
                                     const char      *name,
                                     int              verbosity)
{
  cs_timer_t t0, t1;

  cs_multigrid_level_info_t *mg_lv_info;
  const cs_grid_t *g;
  const cs_matrix_t *m;

  size_t l = strlen(name) + 32;
  char *_name;
  CS_MALLOC(_name, l, char);

  /* Initialization */

  t0 = cs_timer_time();

  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  cs_lnum_t stride = 1; /* For diagonal blocks */

  /* Prepare solver context */

  assert(mgd->n_levels == 2);

  unsigned i = 0;

  assert(mg->subtype == CS_MULTIGRID_BOTTOM);
  assert(mg->type == CS_MULTIGRID_V_CYCLE);

  g = mgd->grid_hierarchy[i];
  m = cs_grid_get_matrix(g);

  mg_lv_info = mg->lv_info + i;

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t0, &t1);

  /* Intermediate grids */

  t0 = t1;

  g = mgd->grid_hierarchy[0];
  m = cs_grid_get_matrix(g);

  mg_lv_info = mg->lv_info;

  cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[0]);

  if (mg->info.type[0] < CS_SLES_N_SMOOTHER_TYPES) {
    mg_sles->context
      = cs_multigrid_smoother_create(mg->info.type[0],
                                     mg->info.poly_degree[0],
                                     mg->info.n_max_iter[0]);
    mg_sles->setup_func = cs_multigrid_smoother_setup;
    mg_sles->solve_func = cs_multigrid_smoother_solve;
    mg_sles->destroy_func = cs_sles_it_destroy;
  }
  else if (i == 0 && mg->lv_mg[0] != nullptr) {
    assert(mg_sles->context == nullptr);
    mg_sles->context = mg->lv_mg[0];
    mg_sles->setup_func = _multigrid_setup_k_local_smoothe;
    mg_sles->solve_func = cs_multigrid_solve;
    mg_sles->destroy_func = nullptr;
  }

  snprintf(_name, l-1, "%s:smoother:%d", name, i);
  _name[l-1] = '\0';

  mg_sles->setup_func(mg_sles->context, _name, m, verbosity - 2);
#if defined(HAVE_MPI)
  if (mg_sles->solve_func == cs_sles_it_solve) {
    cs_sles_it_t  *context = (cs_sles_it_t  *)(mg_sles->context);
    cs_sles_it_set_mpi_reduce_comm(context,
                                   MPI_COMM_NULL,
                                   MPI_COMM_NULL);
  }
#endif

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t0, &t1);

  /* Coarsest grid */

  t0 = t1;

  g = mgd->grid_hierarchy[1];
  m = cs_grid_get_matrix(g);

  mg_lv_info = mg->lv_info + 1;

  mg_sles = &(mgd->sles_hierarchy[2]);
  mg_sles->context
    = cs_sles_it_create(mg->info.type[2],
                        mg->info.poly_degree[2],
                        mg->info.n_max_iter[2],
                        false); /* stats not updated here */
  mg_sles->setup_func = cs_sles_it_setup;
  mg_sles->solve_func = cs_sles_it_solve;
  mg_sles->destroy_func = cs_sles_it_destroy;

  if (mg->lv_mg[2] != nullptr) {
    cs_sles_pc_t *pc = _pc_create_from_mg_sub(mg->lv_mg[2]);
    cs_sles_it_transfer_pc((cs_sles_it_t *)mg_sles->context, &pc);
  }

#if defined(HAVE_MPI)
  {
    cs_sles_it_t  *context = (cs_sles_it_t *)(mg_sles->context);
    cs_sles_it_set_mpi_reduce_comm(context,
                                   cs_grid_get_comm(mgd->grid_hierarchy[1]),
                                   mg->comm);
  }
#endif

  snprintf(_name, l-1, "%s:coarse:%d", name, i);
  _name[l-1] = '\0';

  mg_sles->setup_func(mg_sles->context, _name, m, verbosity - 2);

  /* Diagonal block size is the same for all levels */

  stride = cs_matrix_get_diag_block_size(m);

  _multigrid_setup_sles_work_arrays(mg, stride);

  CS_FREE(_name);

  /* Timing */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t0, &t1);
}

/*----------------------------------------------------------------------------
 * Setup multigrid sparse linear equation solvers on existing hierarchy.
 *
 * parameters:
 *   mg        <-> pointer to multigrid solver info and context
 *   name      <-- linear system name
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

static void
_multigrid_setup_sles(cs_multigrid_t  *mg,
                      const char      *name,
                      int              verbosity)
{
  cs_timer_t t0, t1;

  cs_multigrid_level_info_t *mg_lv_info;
  const cs_grid_t *g;
  const cs_matrix_t *m;

  size_t l = strlen(name) + 32;
  char *_name;
  CS_MALLOC(_name, l, char);

  /* Initialization */

  t0 = cs_timer_time();

  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  cs_lnum_t stride = 1; /* For diagonal blocks */

  /* Prepare solver context */

  unsigned n_levels = mgd->n_levels;

  unsigned i = 0;

  g = mgd->grid_hierarchy[i];
  m = cs_grid_get_matrix(g);

  mg_lv_info = mg->lv_info + i;

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t0, &t1);

  /* Intermediate grids */

  for (i = 0; i < n_levels - 1; i++) {

    t0 = t1;

    g = mgd->grid_hierarchy[i];
    m = cs_grid_get_matrix(g);

    mg_lv_info = mg->lv_info + i;

    int n_ops = 2;
    if (i == 0) {
      if (mg->type == CS_MULTIGRID_V_CYCLE && mg->info.n_max_cycles > 1)
        n_ops = 1;
    }

    int k = 0;
    if (   (cs_get_device_id() > -1 || _force_gpu_settings)
        && i <= _grid_max_level_for_device)
      k = 3;

    for (int j = 0; j < n_ops; j++) {
      if (   mg->info.type[j+k] != CS_SLES_N_IT_TYPES
          && mg->info.type[j+k] < CS_SLES_N_SMOOTHER_TYPES) {
        cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[i*2 + j]);
        mg_sles->context = nullptr;
        if (mg->info.precision_mult[j] < 0) {
          mg_sles->context
            = cs_multigrid_smoother_create(mg->info.type[j+k],
                                           mg->info.poly_degree[j+k],
                                           mg->info.n_max_iter[j+k]);
          if (mg_sles->context != nullptr) {
            mg_sles->setup_func = cs_multigrid_smoother_setup;
            mg_sles->solve_func = cs_multigrid_smoother_solve;
            mg_sles->destroy_func = cs_sles_it_destroy;
          }
        }
        if (mg_sles->context == nullptr) {
          mg_sles->context
            = cs_sles_it_create(mg->info.type[j+k],
                                mg->info.poly_degree[j+k],
                                mg->info.n_max_iter[j+k],
                                false); /* stats not updated here */
          mg_sles->setup_func = cs_sles_it_setup;
          mg_sles->solve_func = cs_sles_it_solve;
          mg_sles->destroy_func = cs_sles_it_destroy;
        }

        /* Share context between descent and ascent smoothers if both
           are of the cs_sles_it type */
        if (j == 1) {
          if (   mg->info.type[0] != CS_SLES_N_IT_TYPES
              && mg->info.type[0] < CS_SLES_N_SMOOTHER_TYPES) {
            cs_sles_it_t *c1
              = (cs_sles_it_t *)(mgd->sles_hierarchy[i*2 + 1].context);
            cs_sles_it_t *c0
              = (cs_sles_it_t *)(mgd->sles_hierarchy[i*2].context);
            cs_sles_it_set_shareable(c1, c0);
          }
        }
      }
    }

    if (i == 0 && mg->lv_mg[0] != nullptr) {
      cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[0]);
      assert(mg->subtype == CS_MULTIGRID_BOTTOM);
      assert(mg->type == CS_MULTIGRID_V_CYCLE);
      assert(n_ops == 1);
      assert(mg_sles->context == nullptr);
      mg_sles->context = mg->lv_mg[0];
      mg_sles->setup_func = _multigrid_setup_k_local_smoothe;
      mg_sles->solve_func = cs_multigrid_solve;
      mg_sles->destroy_func = nullptr;
    }

    for (int j = 0; j < n_ops; j++) {
      if (j == 0)
        snprintf(_name, l-1, "%s:descent:%d", name, i);
      else
        snprintf(_name, l-1, "%s:ascent:%d", name, i);
      _name[l-1] = '\0';
      cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[i*2 + j]);
      mg_sles->setup_func(mg_sles->context, _name, m, verbosity - 2);
#if defined(HAVE_MPI)
      if (mg_sles->solve_func == cs_sles_it_solve) {
        cs_sles_it_t  *context = (cs_sles_it_t *)mg_sles->context;
        MPI_Comm lv_comm = cs_grid_get_comm(mgd->grid_hierarchy[i]);
        cs_sles_it_set_mpi_reduce_comm(context,
                                       lv_comm,
                                       mg->comm);
      }
#endif
    }

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t0, &t1);

  }

  /* Coarsest grid */

  if (n_levels > 1) {

    t0 = t1;

    i = n_levels - 1;

    int k = 0;
    if (   (cs_get_device_id() > -1 || _force_gpu_settings)
        && i <= _grid_max_level_for_device)
      k = 3;

    g = mgd->grid_hierarchy[i];
    m = cs_grid_get_matrix(g);

    mg_lv_info = mg->lv_info + i;

    cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[i*2]);

    mg_sles->context = nullptr;
    if (mg->info.precision_mult[2] < 0) {
      mg_sles->context
        = cs_multigrid_smoother_create(mg->info.type[2+k],
                                       mg->info.poly_degree[2+k],
                                       mg->info.n_max_iter[2+k]);
      if (mg_sles->context != nullptr) {
        mg_sles->setup_func = cs_multigrid_smoother_setup;
        mg_sles->solve_func = cs_multigrid_smoother_solve;
        mg_sles->destroy_func = cs_sles_it_destroy;
      }
    }
    if (mg_sles->context == nullptr) {
      mg_sles->context
        = cs_sles_it_create(mg->info.type[2+k],
                            mg->info.poly_degree[2+k],
                            mg->info.n_max_iter[2+k],
                            false); /* stats not updated here */
      mg_sles->setup_func = cs_sles_it_setup;
      mg_sles->solve_func = cs_sles_it_solve;
      mg_sles->destroy_func = cs_sles_it_destroy;
    }

    if (mg->lv_mg[2] != nullptr) {
      cs_sles_pc_t *pc = _pc_create_from_mg_sub(mg->lv_mg[2]);
      cs_sles_it_transfer_pc((cs_sles_it_t *)mg_sles->context, &pc);
    }

#if defined(HAVE_MPI)
    {
      cs_sles_it_t  *context = (cs_sles_it_t *)mg_sles->context;
      cs_sles_it_set_mpi_reduce_comm(context,
                                     cs_grid_get_comm(mgd->grid_hierarchy[i]),
                                     mg->comm);
    }
#endif

    snprintf(_name, l-1, "%s:coarse:%d", name, i);
    _name[l-1] = '\0';

    mg_sles->setup_func(mg_sles->context, _name, m, verbosity - 2);

    /* Diagonal block size is the same for all levels */

    stride = cs_matrix_get_diag_block_size(m);

  }

  _multigrid_setup_sles_work_arrays(mg, stride);

  CS_FREE(_name);

  /* Timing */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t0, &t1);
}

/*----------------------------------------------------------------------------
 * Compute buffer size required for level names
 *
 * parameters:
 *   name     <-- linear system name
 *   n_levels <-- number multigrid levels
 *
 * returns:
 *   buffer size needed for level names
 *----------------------------------------------------------------------------*/

static size_t
_level_names_size(const char  *name,
                  int          n_levels)
{
  /* Format name width */

  int w = 1;
  for (int i = n_levels/10; i > 0; i /=10)
    w += 1;

  /* First part: pointers */

  size_t retval = n_levels*sizeof(char *)*2;
  retval = CS_SIMD_SIZE(retval);

  /* Second part: buffers */
  size_t buf_size = 0;

  if (n_levels > 1)
    buf_size =   (strlen(name) + strlen(":descent:") + w + 1)
               * (n_levels)*2;
  retval += CS_SIMD_SIZE(buf_size);

  return retval;
}

/*----------------------------------------------------------------------------
 * Initialize level names
 *
 * parameters:
 *   name <-- linear system name
 *   mgd  <-- multigrid setup data
 *----------------------------------------------------------------------------*/

static void
_level_names_init(const char                 *name,
                  cs_multigrid_setup_data_t  *mgd)
{
  int n_levels = mgd->n_levels;
  size_t lv_names_size = _level_names_size(name, n_levels);
  CS_REALLOC(mgd->lv_names_buffer, lv_names_size, char);

  /* Format name width */

  int w = 1;
  for (int i = n_levels/10; i > 0; i /=10)
    w += 1;

  /* First part: pointers */

  size_t ptr_size = n_levels*sizeof(char *)*2;
  ptr_size = CS_SIMD_SIZE(ptr_size);

  char *_buffer = (char *)(mgd->lv_names_buffer);
  char **_lv_names = (char **)(mgd->lv_names_buffer);
  const size_t name_len = strlen(name) + strlen(":descent:") + w + 1;

  /* Second part: buffers */

  for (short int i = 0; i < (short)n_levels -1; i++) {
    _lv_names[i*2] = _buffer + ptr_size + i*2*name_len;
    _lv_names[i*2+1] = _lv_names[i*2] + name_len;
    snprintf(_lv_names[i*2], name_len, "%s:descent:%0*hd", name, w, i);
    snprintf(_lv_names[i*2+1], name_len, "%s:ascent:%0*hd", name, w, i);
  }

  if (n_levels > 1) {
    short int i = n_levels - 1;
    _lv_names[i*2] = _buffer + ptr_size + i*2*name_len;
    _lv_names[i*2+1] = nullptr;
    snprintf(_lv_names[i*2], name_len, "%s:coarse:%0*hd", name, w, i);
  }

  mgd->lv_names = const_cast<const char **>(_lv_names);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup multigrid sparse linear equation solver.
 *
 * \param[in, out]  context    pointer to multigrid solver info and context
 *                             (actual type: cs_multigrid_t  *)
 * \param[in]       name       pointer to name of linear system
 * \param[in]       mesh       associated mesh (for visualization), or nullptr
 * \param[in, out]  f          associated fine grid
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

static void
_setup_hierarchy(void             *context,
                 const char       *name,
                 const cs_mesh_t  *mesh,
                 cs_grid_t        *f,
                 int               verbosity)

{
  cs_multigrid_t  *mg = (cs_multigrid_t *)context;

  cs_timer_t t0, t1, t2;

  int n_coarse_ranks = -2, n_coarse_ranks_prev = -2; /* for comparison only */
  cs_lnum_t n_rows = 0, n_cols_ext = 0, n_entries = 0;
  cs_gnum_t n_g_rows = 0, n_g_rows_prev = 0;

  cs_multigrid_level_info_t *mg_lv_info = nullptr;

  cs_grid_t *g = f;

  t0 = cs_timer_time();

  /* Initialization */

  mg->setup_data = _multigrid_setup_data_create();

  _multigrid_add_level(mg, f); /* Assign to hierarchy */

  /* Add info */

  cs_grid_get_info(f,
                   nullptr,
                   nullptr,
                   nullptr,
                   nullptr,
                   nullptr,
                   &n_rows,
                   &n_cols_ext,
                   &n_entries,
                   &n_g_rows);

  mg_lv_info = mg->lv_info;

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t0, &t1);

  bool add_grid = true;

  while (add_grid) {

    n_g_rows_prev = n_g_rows;
    n_coarse_ranks_prev = n_coarse_ranks;

    /* Recursion test */

    if ((int)(mg->setup_data->n_levels) >= mg->n_levels_max)
      break;

    /* Build coarser grid from previous grid */

    if (verbosity > 2)
      bft_printf(_("\n   building level %2u grid\n"), mg->setup_data->n_levels);

    cs_alloc_mode_t amode = cs_grid_get_alloc_mode(g);
    if (mg->setup_data->n_levels > _grid_max_level_for_device)
      amode = CS_ALLOC_HOST;

    if (mg->subtype == CS_MULTIGRID_BOTTOM)
      g = cs_grid_coarsen_to_single(g, amode, mg->merge_stride, verbosity);

    else {

      int grid_lv = mg->setup_data->n_levels;
      int fg_i = (grid_lv < mg->f_settings_threshold) ? 0 : 1;

      cs_grid_coarsening_t coarsening_type = mg->coarsening_type[fg_i];
      int aggregation_limit = mg->aggregation_limit[fg_i];

      g = cs_grid_coarsen(g,
                          amode,
                          coarsening_type,
                          aggregation_limit,
                          verbosity,
                          mg->merge_stride,
                          mg->merge_mean_threshold,
                          mg->merge_glob_threshold,
                          mg->p0p1_relax);
    }

    bool symmetric = true;
    int grid_lv;

    cs_grid_get_info(g,
                     &grid_lv,
                     &symmetric,
                     nullptr,
                     nullptr,
                     &n_coarse_ranks,
                     &n_rows,
                     &n_cols_ext,
                     &n_entries,
                     &n_g_rows);

#if defined(HAVE_MPI)
    if ((n_coarse_ranks != mg->caller_n_ranks) && (mg->caller_n_ranks > 1)) {
      cs_gnum_t _n_g_rows = n_g_rows;
      MPI_Allreduce(&_n_g_rows, &n_g_rows, 1, CS_MPI_GNUM, MPI_MAX, mg->caller_comm);
    }
#endif

    assert((unsigned)grid_lv == mg->setup_data->n_levels);

    /* If too few rows were grouped, we stop at this level */

    if (n_g_rows == 0)
      add_grid = false;
    else if (mg->setup_data->n_levels > 1) {
      if (   (n_g_rows < mg->n_g_rows_min)
          || (   n_g_rows > (0.8 * n_g_rows_prev)
              && n_coarse_ranks == n_coarse_ranks_prev)
          || n_g_rows == 0) {
        add_grid = false;
      }
    }

    if (add_grid) {

      _multigrid_add_level(mg, g); /* Assign to hierarchy */

      /* Print coarse mesh stats */

      if (verbosity > 2) {

#if defined(HAVE_MPI)

        if (mg->caller_n_ranks > 1) {

          int lcount[4], gcount[4];
          int n_c_min, n_c_max, n_f_min, n_f_max;

          lcount[0] = n_rows; lcount[1] = n_entries;
          if (n_rows > 0) {
            lcount[2] = -n_rows; lcount[3] = -n_entries;
          }
          else {
            lcount[2] = -INT_MAX; lcount[3] = -INT_MAX;
          }
          MPI_Allreduce(lcount, gcount, 4, MPI_INT, MPI_MAX,
                        mg->caller_comm);
          n_c_max = gcount[0]; n_f_max = gcount[1];
          n_c_min = -gcount[2]; n_f_min = -gcount[3];

          bft_printf
            (_("                                  total       min        max\n"
               "     number of rows:      %12llu %10d %10d\n"
               "     number of entries:                %10d %10d\n"),
             (unsigned long long)n_g_rows, n_c_min, n_c_max, n_f_min, n_f_max);
      }

#endif

        if (mg->caller_n_ranks == 1)
          bft_printf(_("     number of rows:      %10d\n"
                       "     number of entries:   %10d\n"),
                     (int)n_rows, (int)n_entries);

      }

      mg_lv_info = mg->lv_info + grid_lv;
      mg_lv_info->n_ranks[0] = n_coarse_ranks;
      mg_lv_info->n_elts[0][0] = n_rows;
      mg_lv_info->n_elts[1][0] = n_cols_ext;
      mg_lv_info->n_elts[2][0] = n_entries;

    } /* end adding grid */

    else
      cs_grid_destroy(&g); /* too coarse, not added */

    t2 = cs_timer_time();
    cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t1, &t2);
    t1 = t2;

  }

  if (   mg->merge_bottom_n_max_ranks < cs_glob_n_ranks
      || cs_glob_n_ranks == 1) {
    cs_multigrid_setup_data_t *mgd = mg->setup_data;
    cs_grid_merge_bottom(mgd->grid_hierarchy[mgd->n_levels-1],
                         verbosity,
                         mg->merge_bottom_n_max_ranks,
                         mg->merge_bottom_max_row_factor);

    t2 = cs_timer_time();
    cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t1, &t2);
    t1 = t2;
  }

  _update_level_info(mg);

  /* Print final info */

  if (verbosity > 1)
    bft_printf
      (_("   number of grid levels:           %u\n"
         "   number of rows in coarsest grid: %llu\n\n"),
       mg->setup_data->n_levels, (unsigned long long)n_g_rows);

  /* Initialize names for later logging */

  _level_names_init(name, mg->setup_data);

  /* Prepare preprocessing info if necessary */

  if (mg->post_row_max > 0) {
    if (mg->info.n_calls[0] == 0) {
      int l_id = 0;
      const cs_matrix_t *a = cs_grid_get_matrix(f);
      const cs_lnum_t n_rows_a = cs_matrix_get_n_rows(a);
      if (n_rows_a == mesh->n_cells)
        l_id = CS_MESH_LOCATION_CELLS;
      else if (n_rows_a <= mesh->n_vertices) {
        l_id = CS_MESH_LOCATION_VERTICES;
        if (mesh->vtx_range_set != nullptr) {
          if (mesh->vtx_range_set->n_elts[0] != n_rows_a)
            l_id = CS_MESH_LOCATION_NONE;
        }
      }
#if defined(HAVE_MPI)
      if (mg->caller_n_ranks > 1) {
        int _l_id = l_id;
        MPI_Allreduce(&_l_id, &l_id, 1, MPI_INT, MPI_MAX, mg->caller_comm);
        if (l_id != _l_id)
          _l_id = CS_MESH_LOCATION_NONE;
        MPI_Allreduce(&_l_id, &l_id, 1, MPI_INT, MPI_MIN, mg->caller_comm);
      }
#endif
      if (l_id != CS_MESH_LOCATION_NONE) {
        cs_post_add_time_dep_output(_cs_multigrid_post_function, (void *)mg);
        _multigrid_add_post(mg, name, l_id, n_rows_a);
      }
    }
  }

  /* Update info */

#if defined(HAVE_MPI)

  /* In parallel, get global (average) values from local values */

  if (mg->caller_n_ranks > 1) {

    int i, j;
    cs_gnum_t *_n_elts_l = nullptr, *_n_elts_s = nullptr, *_n_elts_m = nullptr;
    int grid_lv = mg->setup_data->n_levels;

    CS_MALLOC(_n_elts_l, 3*grid_lv, cs_gnum_t);
    CS_MALLOC(_n_elts_s, 3*grid_lv, cs_gnum_t);
    CS_MALLOC(_n_elts_m, 3*grid_lv, cs_gnum_t);

    for (i = 0; i < grid_lv; i++) {
      cs_multigrid_level_info_t *mg_inf = mg->lv_info + i;
      for (j = 0; j < 3; j++)
        _n_elts_l[i*3 + j] = mg_inf->n_elts[j][0];
    }

    MPI_Allreduce(_n_elts_l, _n_elts_s, 3*grid_lv, CS_MPI_GNUM, MPI_SUM,
                  mg->caller_comm);
    MPI_Allreduce(_n_elts_l, _n_elts_m, 3*grid_lv, CS_MPI_GNUM, MPI_MAX,
                  mg->caller_comm);

    for (i = 0; i < grid_lv; i++) {
      cs_multigrid_level_info_t *mg_inf = mg->lv_info + i;
      cs_gnum_t n_g_ranks = mg_inf->n_ranks[0];
      for (j = 0; j < 3; j++) {
        cs_gnum_t tmp_max = n_g_ranks * _n_elts_m[i*3+j];
        mg_inf->n_elts[j][0] = (_n_elts_s[i*3+j] + n_g_ranks/2) / n_g_ranks;
        if (_n_elts_s[i*3+j] > 0)
          mg_inf->imbalance[j][0] = (float)(tmp_max*1.0/_n_elts_s[i*3+j]);
        else
          mg_inf->imbalance[j][0] = 0;
      }
    }

    CS_FREE(_n_elts_m);
    CS_FREE(_n_elts_s);
    CS_FREE(_n_elts_l);

  }

#endif

  mg->info.n_levels_tot += mg->setup_data->n_levels;

  mg->info.n_levels[0] = mg->setup_data->n_levels;

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

  /* Cleanup temporary interpolation arrays */

  for (unsigned i = 0; i < mg->setup_data->n_levels; i++)
    cs_grid_free_quantities(mg->setup_data->grid_hierarchy[i]);

  /* Setup solvers */

  if (mg->subtype == CS_MULTIGRID_BOTTOM)
    _multigrid_setup_sles_k_cycle_bottom(mg, name, verbosity);
  else
    _multigrid_setup_sles(mg, name, verbosity);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup coarse multigrid for k cycle HPC variant.
 *
 * \param[in, out]  mg         pointer to multigrid solver info and context
 * \param[in]       name       pointer to name of linear system
 * \param[in]       a          associated matrix
 * \param[in]       n_ranks    associated number of MPI ranks
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

static void
_setup_k_cycle_hpc_sub(cs_multigrid_t     *mg,
                       const char         *name,
                       const cs_matrix_t  *a,
                       int                 n_ranks,
                       int                 verbosity)
{
  /* Destroy previous hierarchy if necessary */

  if (mg->setup_data != nullptr)
    cs_multigrid_free(mg);

  /* Initialization */

  cs_timer_t t0 = cs_timer_time();

  if (verbosity > 1)
    bft_printf(_("\n Construction of grid hierarchy for \"%s\"\n"),
               name);

  /* Build coarse grids hierarchy */
  /*------------------------------*/

  cs_grid_t *f = cs_grid_create_from_parent(a, n_ranks);

  cs_multigrid_level_info_t *mg_lv_info = mg->lv_info;

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t0, &t1);

  _setup_hierarchy(mg, name, nullptr, f, verbosity); /* Assign to and build
                                                     hierarchy */

  /* Update timers */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg->info.t_tot[0]), &t0, &t1);
}

/*----------------------------------------------------------------------------
 * Compute dot product, summing result over all ranks.
 *
 * parameters:
 *   mg <-- pointer to solver context info
 *   n  <-- local number of elements
 *   x  <-- vector in s = x.x
 *
 * returns:
 *   result of s = x.x
 *----------------------------------------------------------------------------*/

inline static double
_dot_xx(const cs_multigrid_t  *mg,
        cs_lnum_t              n,
        const cs_real_t       *x)
{
  double s = cs_dot_xx(n, x);

#if defined(HAVE_MPI)

  if (mg->comm != MPI_COMM_NULL) {
    double _sum;
    MPI_Allreduce(&s, &_sum, 1, MPI_DOUBLE, MPI_SUM, mg->comm);
    s = _sum;
  }

#endif /* defined(HAVE_MPI) */

  return s;
}

/*----------------------------------------------------------------------------
 * Compute 2 dot products x.x and y.y, summing result over all ranks.
 *
 * parameters:
 *   mg <-- pointer to solver context info
 *   n  <-- number of associated values
 *   x  <-- vector in s1 = x.x
 *   y  <-- vector in s2 = y.y
 *   s1 --> result of s1 = x.x
 *   s2 --> result of s2 = y.y
 *----------------------------------------------------------------------------*/

inline static void
_dot_xx_yy(const cs_multigrid_t  *mg,
           cs_lnum_t              n,
           const cs_real_t       *x,
           const cs_real_t       *y,
           double                *s1,
           double                *s2)
{
  double s[2];

  s[0] = cs_dot_xx(n, x);
  s[1] = cs_dot_xx(n, y);

#if defined(HAVE_MPI)

  if (mg->comm != MPI_COMM_NULL) {
    double _sum[2];
    MPI_Allreduce(s, _sum, 2, MPI_DOUBLE, MPI_SUM, mg->comm);
    s[0] = _sum[0];
    s[1] = _sum[1];
  }

#endif /* defined(HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
}

/*----------------------------------------------------------------------------
 * Compute 2 dot products x.x and x.y, summing result over all ranks.
 *
 * parameters:
 *   mg <-- pointer to solver context info
 *   n  <-- number of associated values
 *   x  <-- vector in s1 = x.y
 *   y  <-- vector in s1 = x.y and s2 = y.z
 *   z  <-- vector in s2 = y.z
 *   s1 --> result of s1 = x.y
 *   s2 --> result of s2 = y.z
 *----------------------------------------------------------------------------*/

inline static void
_dot_xy_yz(const cs_multigrid_t  *mg,
           cs_lnum_t              n,
           const cs_real_t       *x,
           const cs_real_t       *y,
           const cs_real_t       *z,
           double                *s1,
           double                *s2)
{
  double s[2];

  cs_dot_xy_yz(n, x, y, z, s, s+1);

#if defined(HAVE_MPI)

  if (mg->comm != MPI_COMM_NULL) {
    double _sum[2];
    MPI_Allreduce(s, _sum, 2, MPI_DOUBLE, MPI_SUM, mg->comm);
    s[0] = _sum[0];
    s[1] = _sum[1];
  }

#endif /* defined(HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
}

/*----------------------------------------------------------------------------
 * Compute 3 dot products x.u, xv, and x.w, summing result over all ranks.
 *
 * parameters:
 *   mg <-- pointer to solver context info
 *   n  <-- number of associated values
 *   x  <-- vector in s1 = x.u, s2 = x.v, s3 = x.w
 *   u <-- vector in s1 = x.u
 *   v <-- vector in s2 = x.v
 *   w <-- vector in s2 = x.w
 *   s1 --> result of s1 = x.y
 *   s2 --> result of s2 = y.z
 *----------------------------------------------------------------------------*/

inline static void
_dot_xu_xv_xw(const cs_multigrid_t  *mg,
              cs_lnum_t              n,
              const cs_real_t       *x,
              const cs_real_t       *u,
              const cs_real_t       *v,
              const cs_real_t       *w,
              double                *s1,
              double                *s2,
              double                *s3)
{
  double s[3];

  /* Use two separate call as cs_blas.c does not yet hav matching call */
  cs_dot_xy_yz(n, u, x, v, s, s+1);
  s[2] = cs_dot(n, x, w);

#if defined(HAVE_MPI)

  if (mg->comm != MPI_COMM_NULL) {
    double _sum[3];
    MPI_Allreduce(s, _sum, 3, MPI_DOUBLE, MPI_SUM, mg->comm);
    s[0] = _sum[0];
    s[1] = _sum[1];
    s[2] = _sum[2];
  }

#endif /* defined(HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
  *s3 = s[2];
}

/*----------------------------------------------------------------------------
 * Test if convergence is attained.
 *
 * parameters:
 *   mg               <-- associated multigrid structure
 *   var_name         <-- variable name
 *   n_f_rows         <-- number of rows on fine mesh
 *   n_max_cycles     <-- maximum number of cycles
 *   cycle_id         <-- number of current cycle
 *
 *   verbosity        <-- verbosity level
 *   n_iters          <-- number of iterations
 *   precision        <-- precision limit
 *   r_norm           <-- residual normalization
 *   initial_residual <-- initial residual
 *   residual         <-> residual
 *   rhs              <-- right-hand side
 *
 * returns:
 *   convergence status
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_convergence_test(cs_multigrid_t        *mg,
                  const char            *var_name,
                  cs_lnum_t              n_f_rows,
                  int                    n_max_cycles,
                  int                    cycle_id,
                  int                    verbosity,
                  int                    n_iters,
                  double                 precision,
                  double                 r_norm,
                  double                 initial_residual,
                  double                *residual,
                  const cs_real_t        rhs[])
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

  /* Compute residual */

  *residual = sqrt(_dot_xx(mg, n_f_rows, rhs));

  if (cycle_id == 1)
    initial_residual = *residual;

  /* Plot convergence if requested */

  if (mg->cycle_plot != nullptr) {
    double vals = *residual;
    double wall_time = cs_timer_wtime();
    mg->plot_time_stamp += 1;
    cs_time_plot_vals_write(mg->cycle_plot,
                            mg->plot_time_stamp,
                            wall_time,
                            1,
                            &vals);
  }

  if (*residual < precision*r_norm) {

    if (verbosity == 2)
      bft_printf(_(cycle_fmt), cycle_id, n_iters, *residual/r_norm);
    else if (verbosity > 2) {
      bft_printf(_(cycle_h_fmt));
      bft_printf(_(cycle_cv_fmt),
                 cycle_id, n_iters, *residual/r_norm);
      bft_printf(_(cycle_t_fmt));
    }
    return CS_SLES_CONVERGED;
  }

  else if (cycle_id > n_max_cycles) {

    if (  (verbosity > -1 && !(mg->info.is_pc || mg->subtype != CS_MULTIGRID_MAIN))
        || verbosity > 0) {
      if (verbosity == 1)
        bft_printf(_(cycle_fmt), cycle_id, n_iters, *residual/r_norm);
      else if (verbosity > 1) {
        bft_printf(_(cycle_fmt),
                   cycle_id, n_iters, *residual/r_norm);
        bft_printf(_(cycle_t_fmt));
      }
      bft_printf(_(" @@ Warning: algebraic multigrid for [%s]\n"
                   "    ********\n"
                   "    Maximum number of cycles (%d) reached.\n"),
                 var_name, n_max_cycles);
    }
    return CS_SLES_MAX_ITERATION;
  }

  else {

    if (*residual > initial_residual * 10000.0 && *residual > 100.) {
      if (verbosity > 2)
        bft_printf(_(cycle_fmt), cycle_id, n_iters, *residual/r_norm);
      return CS_SLES_DIVERGED;
    }
    else if (verbosity > 2) {
      if (cycle_id == 1)
        bft_printf(_(cycle_h_fmt));
      bft_printf(_(cycle_cv_fmt), cycle_id, n_iters, *residual/r_norm);
    }

#if (__STDC_VERSION__ >= 199901L)
    if (isnan(*residual) || isinf(*residual))
      return CS_SLES_DIVERGED;
#endif
  }

  return CS_SLES_ITERATING;
}

/*----------------------------------------------------------------------------
 * Log residual A.vx - Rhs
 *
 * parameters:
 *   mg              <-- pointer to multigrid context info
 *   int cycle_id    <-- cycle id
 *   var_name        <-- variable name
 *   a               <-- matrix
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static void
_log_residual(const cs_multigrid_t   *mg,
              int                     cycle_id,
              const char             *var_name,
              const cs_matrix_t      *a,
              const cs_real_t        *rhs,
              cs_real_t              *restrict vx)
{
  const cs_lnum_t diag_block_size = cs_matrix_get_diag_block_size(a);
  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a) * diag_block_size;
  const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;

  cs_real_t  *r;
  CS_MALLOC(r, n_cols, cs_real_t);

  cs_matrix_vector_multiply(a, vx, r);

  for (cs_lnum_t i = 0; i < n_rows; i++)
    r[i] -= rhs[i];

  double s = cs_dot_xx(n_rows, r);

  CS_FREE(r);

#if defined(HAVE_MPI)

  if (mg->comm != MPI_COMM_NULL) {
    double _sum;
    MPI_Allreduce(&s, &_sum, 1, MPI_DOUBLE, MPI_SUM, mg->comm);
    s = _sum;
  }

#endif /* defined(HAVE_MPI) */

  cs_log_printf(CS_LOG_DEFAULT, "  mg cycle %d: %s residual: %.3g\n",
                cycle_id, var_name, s);
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
  if (n_iter < lv_info_it[1])
    lv_info_it[1] = n_iter;
  else if (n_iter > lv_info_it[2])
    lv_info_it[2] = n_iter;
  if (lv_info_it[1] == 0)
    lv_info_it[1] = n_iter;
  lv_info_it[3] += n_iter;
}

/*----------------------------------------------------------------------------
 * Sparse linear system resolution using multigrid.
 *
 * parameters:
 *   mg               <-- multigrid system
 *   verbosity        <-- verbosity level
 *   cycle_id         <-- id of currect cycle
 *   n_equiv_iter     <-> equivalent number of iterations
 *   precision        <-- solver precision
 *   r_norm           <-- residual normalization
 *   initial_residual <-> initial residual
 *   residual         <-> residual
 *   rhs              <-- right hand side
 *   vx_ini           <-- initial system solution
 *                        (vx if nonzero, nullptr if zero)
 *   vx               --> system solution
 *   aux_size         <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors      --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence status
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_multigrid_v_cycle(cs_multigrid_t       *mg,
                   int                   verbosity,
                   int                   cycle_id,
                   int                  *n_equiv_iter,
                   double                precision,
                   double                r_norm,
                   double               *initial_residual,
                   double               *residual,
                   const cs_real_t      *rhs,
                   cs_real_t            *vx_ini,
                   cs_real_t            *vx,
                   size_t                aux_size,
                   void                 *aux_vectors)
{
  int level, coarsest_level;
  cs_timer_t t0, t1;

  cs_lnum_t db_size = 1;
  cs_lnum_t eb_size = 1;
  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING, c_cvg = CS_SLES_ITERATING;
  int n_iter = 0;
  double _residual = -1.;
  double _initial_residual = 0.;

  size_t _aux_r_size = aux_size / sizeof(cs_real_t);
  cs_lnum_t n_rows = 0, n_cols_ext = 0;
  cs_lnum_t _n_rows = 0;
  cs_gnum_t n_g_rows = 0;
  cs_real_t r_norm_l = r_norm;

  double denom_n_g_rows_0 = 1.0;

  cs_multigrid_setup_data_t *mgd = mg->setup_data;
  cs_multigrid_level_info_t  *lv_info = nullptr;

  cs_real_t *_aux_vectors = (cs_real_t *)aux_vectors;
  cs_real_t *restrict wr = nullptr;
  cs_real_t *restrict vx_lv = nullptr;

  const char **lv_names = mgd->lv_names;

  const cs_real_t *restrict rhs_lv = nullptr;
  const cs_grid_t *f = nullptr, *c= nullptr;

  bool end_cycle = false;

  cs_dispatch_context ctx;

  /* Initialization */

  coarsest_level = mgd->n_levels - 1;

  f = mgd->grid_hierarchy[0];

  cs_grid_get_info(f,
                   nullptr,
                   nullptr,
                   &db_size,
                   &eb_size,
                   nullptr,
                   &n_rows,
                   &n_cols_ext,
                   nullptr,
                   &n_g_rows);

  denom_n_g_rows_0 = 1.0 / n_g_rows;

  const cs_matrix_t *_matrix = cs_grid_get_matrix(f);
  cs_alloc_mode_t amode = cs_matrix_get_alloc_mode(_matrix);

  /* Allocate wr or use working area
     (note the finest grid could have less elements than a coarser
     grid to wich rank merging has been applied, hence the test below) */

  size_t wr_size = n_cols_ext*db_size;
  for (level = 1; level < (int)(mgd->n_levels); level++) {
    cs_lnum_t n_cols_max
      = cs_grid_get_n_cols_max(mgd->grid_hierarchy[level]);
    wr_size = cs::max(wr_size, (size_t)(n_cols_max*db_size));
    wr_size = CS_SIMD_SIZE(wr_size);
  }

  if (_aux_r_size >= wr_size && amode <= CS_ALLOC_HOST) {
    wr = (cs_real_t *)aux_vectors;
    _aux_vectors = wr + wr_size;
    _aux_r_size -= wr_size;
  }
  else {
    if (amode <= CS_ALLOC_HOST)
      CS_MALLOC(wr, wr_size, cs_real_t);
    else
      CS_MALLOC_HD(wr, wr_size, cs_real_t, CS_ALLOC_HOST_DEVICE_SHARED);
  }

  /* map arrays for rhs and vx;
     for the finest level, simply point to input and output arrays */

  mgd->rhs_vx[0] = nullptr; /* Use _rhs_level when necessary to avoid
                               const warning */
  mgd->rhs_vx[1] = vx;

  /* Descent */
  /*---------*/

  for (level = 0; level < coarsest_level; level++) {

    lv_info = mg->lv_info + level;
    t0 = cs_timer_time();

    rhs_lv = (level == 0) ? rhs : mgd->rhs_vx[level*2];
    vx_lv = mgd->rhs_vx[level*2 + 1];
    cs_real_t *restrict vx_lv_ini
      = (level == 0) ? vx_ini : nullptr;

    c = mgd->grid_hierarchy[level+1];

    /* Smoother pass */

    _matrix = cs_grid_get_matrix(f);

#if defined(HAVE_ACCEL)
    cs_alloc_mode_t amode_l = cs_matrix_get_alloc_mode(_matrix);
#endif

    cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[level*2]);

    c_cvg = mg_sles->solve_func(mg_sles->context,
                                lv_names[level*2],
                                _matrix,
                                verbosity - 4, /* verbosity */
                                precision*mg->info.precision_mult[0],
                                r_norm_l,
                                &n_iter,
                                &_residual,
                                rhs_lv,
                                vx_lv_ini,
                                vx_lv,
                                _aux_r_size*sizeof(cs_real_t),
                                _aux_vectors);

    if (mg->plot_time_stamp > -1)
      mg->plot_time_stamp += n_iter+1;

    if (mg_sles->solve_func == cs_sles_it_solve)
      _initial_residual = cs_sles_it_get_last_initial_residual
                            ((const cs_sles_it_t *)mg_sles->context);
    else
      _initial_residual = HUGE_VAL;

    if (level == 0 && cycle_id == 1)
      *initial_residual = _initial_residual;

    if (verbosity > 1)
      _log_residual(mg, cycle_id, lv_names[level*2],
                    _matrix, rhs_lv, vx_lv);

    if (c_cvg < CS_SLES_BREAKDOWN) {
      end_cycle = true;
      break;
    }

    /* Restrict residual
       TODO: get residual from cs_sles_solve(). This optimisation would
       require adding an argument and exercising caution to ensure the
       correct sign and meaning of the residual
       (regarding timing, this stage is part of the descent smoother) */

#if defined(HAVE_ACCEL)
    bool use_gpu = (amode_l > CS_ALLOC_HOST) ? true : false;
    ctx.set_use_gpu(use_gpu);

    if (use_gpu) {
#if defined(HAVE_CUDA)
      cudaStream_t stream = cs_matrix_spmv_cuda_get_stream();
      if (stream != 0)
        ctx.set_cuda_stream(stream);
#endif
      cs_matrix_vector_multiply_d(_matrix, vx_lv, wr);
    }
    else
#endif
      cs_matrix_vector_multiply(_matrix, vx_lv, wr);

    _n_rows = n_rows*db_size;
    ctx.parallel_for(_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      wr[ii] = rhs_lv[ii] - wr[ii];
    });
    ctx.wait();

    /* Convergence test in beginning of cycle (fine mesh) */

    if (level == 0) {

      cvg = _convergence_test(mg,
                              lv_names[0],
                              _n_rows,
                              mg->info.n_max_cycles,
                              cycle_id,
                              verbosity,
                              lv_info->n_it_ds_smoothe[3],
                              precision,
                              r_norm,
                              *initial_residual,
                              residual,
                              wr);

      /* If converged or cycle limit reached, break from descent loop */

      if (cvg != 0) {
        c_cvg = cvg;
        end_cycle = true;
        t1 = cs_timer_time();
        cs_timer_counter_add_diff(&(lv_info->t_tot[2]), &t0, &t1);
        lv_info->n_calls[2] += 1;
        _lv_info_update_stage_iter(lv_info->n_it_ds_smoothe, n_iter);
        *n_equiv_iter += n_iter * n_g_rows * denom_n_g_rows_0;
        break;
      }

    }

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[2]), &t0, &t1);
    lv_info->n_calls[2] += 1;
    _lv_info_update_stage_iter(lv_info->n_it_ds_smoothe, n_iter);
    *n_equiv_iter += n_iter * n_g_rows * denom_n_g_rows_0;

    /* Prepare for next level */

    cs_grid_restrict_row_var(ctx, f, c, wr, mgd->rhs_vx[(level+1)*2]);

    cs_grid_get_info(c,
                     nullptr,
                     nullptr,
                     nullptr,
                     nullptr,
                     nullptr,
                     &n_rows,
                     &n_cols_ext,
                     nullptr,
                     &n_g_rows);

    f = c;
    _n_rows = n_rows*db_size;

    t0 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[4]), &t1, &t0);
    lv_info->n_calls[4] += 1;

  } /* End of loop on levels (descent) */

  if (end_cycle == false) {

    /* Resolve coarsest level to convergence */
    /*---------------------------------------*/

    assert(level == coarsest_level);
    assert(c == mgd->grid_hierarchy[coarsest_level]);

    /* coarsest level == 0 should never happen, but we play it safe */
    rhs_lv = (level == 0) ?  rhs : mgd->rhs_vx[coarsest_level*2];
    vx_lv = mgd->rhs_vx[level*2 + 1];

    _matrix = cs_grid_get_matrix(c);

    cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[level*2]);

    _initial_residual = _residual;

    lv_info = mg->lv_info + level;

    t0 = cs_timer_time();

    c_cvg = mg_sles->solve_func(mg_sles->context,
                                lv_names[level*2],
                                _matrix,
                                verbosity - 3,
                                precision*mg->info.precision_mult[2],
                                r_norm_l,
                                &n_iter,
                                &_residual,
                                rhs_lv,
                                nullptr,
                                vx_lv,
                                _aux_r_size*sizeof(cs_real_t),
                                _aux_vectors);

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[1]), &t0, &t1);
    lv_info->n_calls[1] += 1;
    _lv_info_update_stage_iter(lv_info->n_it_solve, n_iter);

    if (mg_sles->solve_func == cs_sles_it_solve)
      _initial_residual = cs_sles_it_get_last_initial_residual
                            ((const cs_sles_it_t *)mg_sles->context);
    else
      _initial_residual = HUGE_VAL;

    *n_equiv_iter += n_iter * n_g_rows * denom_n_g_rows_0;

    if (verbosity > 1)
      _log_residual(mg, cycle_id, lv_names[level*2],
                    _matrix, rhs_lv, vx_lv);

    if (c_cvg < CS_SLES_BREAKDOWN)
      end_cycle = true;

  }

  if (end_cycle == false) {

    /* Ascent */
    /*--------*/

    for (level = coarsest_level - 1; level > -1; level--) {

      vx_lv = mgd->rhs_vx[level*2 + 1];

      lv_info = mg->lv_info + level;

      c = mgd->grid_hierarchy[level+1];
      f = mgd->grid_hierarchy[level];

      cs_grid_get_info(f,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr,
                       &n_rows,
                       &n_cols_ext,
                       nullptr,
                       &n_g_rows);

      /* Prolong correction */

      t0 = cs_timer_time();

      cs_real_t *restrict vx_lv1 = mgd->rhs_vx[(level+1)*2 + 1];
      cs_grid_prolong_row_var(ctx,
                              c, f,
                              true,  /* increment */
                              vx_lv1, vx_lv);

      _n_rows = n_rows*db_size;

      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&(lv_info->t_tot[5]), &t0, &t1);
      lv_info->n_calls[5] += 1;

      /* Smoother pass if level > 0
         (smoother not called for finest mesh, as it will be called in
         descent phase of the next cycle, before the convergence test). */

      if (level > 0) {

        rhs_lv = mgd->rhs_vx[level*2];

        _matrix = cs_grid_get_matrix(f);

        cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[level*2 + 1]);

        c_cvg = mg_sles->solve_func(mg_sles->context,
                                    lv_names[level*2+1],
                                    _matrix,
                                    verbosity - 4, /* verbosity */
                                    precision*mg->info.precision_mult[1],
                                    r_norm_l,
                                    &n_iter,
                                    &_residual,
                                    rhs_lv,
                                    vx_lv,
                                    vx_lv,
                                    _aux_r_size*sizeof(cs_real_t),
                                    _aux_vectors);

        t0 = cs_timer_time();
        cs_timer_counter_add_diff(&(lv_info->t_tot[3]), &t1, &t0);
        lv_info->n_calls[3] += 1;
        _lv_info_update_stage_iter(lv_info->n_it_as_smoothe, n_iter);

        if (mg_sles->solve_func == cs_sles_it_solve)
          _initial_residual = cs_sles_it_get_last_initial_residual
                                ((const cs_sles_it_t *)mg_sles->context);
        else
          _initial_residual = HUGE_VAL;

        *n_equiv_iter += n_iter * n_g_rows * denom_n_g_rows_0;

        if (verbosity > 1)
          _log_residual(mg, cycle_id, lv_names[level*2+1],
                        _matrix, rhs_lv, vx_lv);

        if (c_cvg < CS_SLES_BREAKDOWN)
          break;
      }

    } /* End loop on levels (ascent) */

  } /* End of tests on end_cycle */

  mgd->exit_level = level;
  mgd->exit_residual = _residual;
  if (level == 0)
    mgd->exit_initial_residual = *initial_residual;
  else
    mgd->exit_initial_residual = _initial_residual;
  mgd->exit_cycle_id = cycle_id;

  /* Free memory */

  if (wr != aux_vectors)
    CS_FREE(wr);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Single multigrid V-cycle for preconditioning.
 *
 * parameters:
 *   mg               <-- multigrid system
 *   verbosity        <-- verbosity level
 *   n_equiv_iter     <-> equivalent number of iterations
 *   precision        <-- solver precision
 *   r_norm           <-- residual normalization
 *   initial_residual <-> initial residual
 *   residual         <-> residual
 *   rhs              <-- right hand side
 *   vx_ini           <-- initial system solution
 *                        (vx if nonzero, nullptr if zero)
 *   vx               --> system solution
 *   aux_size         <-- number of elements in aux_vectors_d (in bytes)
 *   aux_size_h       <-- number of elements in aux_vectors_h (in bytes)
 *   aux_vectors      --- optional working area (allocation otherwise)
 *   aux_vectors_h    --- optional host-only working area (allocation otherwise)
 *
 * returns:
 *   convergence status
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_multigrid_v_cycle_pc(cs_multigrid_t        *mg,
                      int                    verbosity,
                      int                   *n_equiv_iter,
                      double                 precision,
                      double                 r_norm,
                      double                *initial_residual,
                      double                *residual,
                      const cs_real_t       *rhs,
                      cs_real_t             *vx_ini,
                      cs_real_t             *vx,
                      size_t                 aux_size,
                      size_t                 aux_size_h,
                      [[maybe_unused]]void  *aux_vectors,
                      void                  *aux_vectors_h)
{
#if CS_PROFILING == CS_PROFILING_NVTX
  static int call_id = -1;
  char       nvtx_name[48];
  call_id++;
  if (call_id == 2) {
    nvtxMark("Enter multigrid cycle");
  }
#endif

  CS_PROFILE_FUNC_RANGE();

  int level, coarsest_level;
  cs_timer_t t0, t1;

  cs_lnum_t db_size = 1;
  cs_lnum_t eb_size = 1;
  cs_sles_convergence_state_t c_cvg = CS_SLES_ITERATING;
  int n_iter = 0;
  double _residual = -1.;
  double _initial_residual = 0.;

  cs_lnum_t n_rows = 0, n_cols_ext = 0;
  cs_lnum_t _n_rows = 0;
  cs_gnum_t n_g_rows = 0;
  cs_real_t r_norm_l = r_norm;

  double denom_n_g_rows_0 = 1.0;

  cs_multigrid_setup_data_t *mgd = mg->setup_data;
  cs_multigrid_level_info_t  *lv_info = nullptr;

  cs_real_t *restrict vx_lv = nullptr;

  const char **lv_names = mgd->lv_names;

  const int cycle_id = 1;
  const cs_real_t *restrict rhs_lv = nullptr;
  const cs_grid_t *f = nullptr, *c= nullptr;

  bool end_cycle = false;

  cs_dispatch_context ctx;

  /* Initialization */

  coarsest_level = mgd->n_levels - 1;

  f = mgd->grid_hierarchy[0];

  cs_grid_get_info(f,
                   nullptr,
                   nullptr,
                   &db_size,
                   &eb_size,
                   nullptr,
                   &n_rows,
                   &n_cols_ext,
                   nullptr,
                   &n_g_rows);

  denom_n_g_rows_0 = 1.0 / n_g_rows;

  const cs_matrix_t *_matrix = cs_grid_get_matrix(f);
#if defined(HAVE_ACCEL)
  cs_alloc_mode_t amode_f = cs_matrix_get_alloc_mode(_matrix);
#endif

  /* map arrays for rhs and vx;
     for the finest level, simply point to input and output arrays */

  mgd->rhs_vx[0] = nullptr; /* Use _rhs_level when necessary to avoid
                               const warning */
  mgd->rhs_vx[1] = vx;

  /* Descent */
  /*---------*/

  for (level = 0; level < coarsest_level; level++) {

    lv_info = mg->lv_info + level;
    t0 = cs_timer_time();

    rhs_lv = (level == 0) ? rhs : mgd->rhs_vx[level*2];
    vx_lv = mgd->rhs_vx[level*2 + 1];
    cs_real_t *restrict vx_lv_ini
      = (level == 0) ? vx_ini : nullptr;

    c = mgd->grid_hierarchy[level+1];

    /* Smoother pass */

    _matrix = cs_grid_get_matrix(f);

    _n_rows = n_rows*db_size;

    size_t wr_size = n_cols_ext*db_size;
    wr_size = CS_SIMD_SIZE(wr_size);
    size_t _aux_r_size = aux_size_h - wr_size*sizeof(cs_real_t);
    cs_real_t *restrict wr = (cs_real_t *)aux_vectors_h;

#if defined(HAVE_ACCEL)
    cs_alloc_mode_t amode_p = amode_f;
    amode_f = cs_matrix_get_alloc_mode(_matrix);
    if (amode_f > CS_ALLOC_HOST) {
      assert(aux_size >= wr_size*sizeof(cs_real_t));
      _aux_r_size = aux_size - wr_size*sizeof(cs_real_t);
      wr = (cs_real_t *)aux_vectors;
    }
    else {
      if (amode_p > CS_ALLOC_HOST)
        cs_prefetch_d2h(mgd->rhs_vx[level*2], _n_rows*sizeof(cs_real_t));
      assert(aux_size_h >= wr_size*sizeof(cs_real_t));
    }
#else
    assert(aux_size >= wr_size*sizeof(cs_real_t));
#endif

    cs_real_t *_aux_vectors = wr + wr_size;

    cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[level*2]);

#if CS_PROFILING == CS_PROFILING_NVTX
    sprintf(nvtx_name, "Descent smoothe %d", level);
    nvtxRangePushA(nvtx_name);
#endif

    c_cvg = mg_sles->solve_func(mg_sles->context,
                                lv_names[level*2],
                                _matrix,
                                verbosity - 4, /* verbosity */
                                precision*mg->info.precision_mult[0],
                                r_norm_l,
                                &n_iter,
                                &_residual,
                                rhs_lv,
                                vx_lv_ini,
                                vx_lv,
                                _aux_r_size,
                                _aux_vectors);

#if CS_PROFILING == CS_PROFILING_NVTX
    nvtxRangePop();
#endif

    if (mg->plot_time_stamp > -1)
      mg->plot_time_stamp += n_iter+1;

    if (mg_sles->solve_func == cs_sles_it_solve)
      _initial_residual = cs_sles_it_get_last_initial_residual
                            ((const cs_sles_it_t *)mg_sles->context);
    else
      _initial_residual = HUGE_VAL;

    if (level == 0)
      *initial_residual = _initial_residual;

    if (verbosity > 1)
      _log_residual(mg, cycle_id, lv_names[level*2],
                    _matrix, rhs_lv, vx_lv);

    if (c_cvg < CS_SLES_BREAKDOWN) {
      end_cycle = true;
      break;
    }

    /* Restrict residual
       TODO: get residual from cs_sles_solve(). This optimisation would
       require adding an argument and exercising caution to ensure the
       correct sign and meaning of the residual
       (regarding timing, this stage is part of the descent smoother) */

#if defined(HAVE_ACCEL)
    bool use_gpu = (amode_f > CS_ALLOC_HOST) ? true : false;
    ctx.set_use_gpu(use_gpu);

    if (use_gpu) {
#if defined(HAVE_CUDA)
      cudaStream_t stream = cs_matrix_spmv_cuda_get_stream();
      if (stream != 0)
        ctx.set_cuda_stream(stream);
#endif
      cs_matrix_vector_multiply_d(_matrix, vx_lv, wr);
    }
    else
#endif
      cs_matrix_vector_multiply(_matrix, vx_lv, wr);

    ctx.parallel_for(_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      wr[ii] = rhs_lv[ii] - wr[ii];
    });
    ctx.wait();

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[2]), &t0, &t1);
    lv_info->n_calls[2] += 1;
    _lv_info_update_stage_iter(lv_info->n_it_ds_smoothe, n_iter);
    *n_equiv_iter += n_iter * n_g_rows * denom_n_g_rows_0;

    /* Prepare for next level */

#if CS_PROFILING == CS_PROFILING_NVTX
    sprintf(nvtx_name, "Restrict %d to %d", level, level +1);
    nvtxRangePushA(nvtx_name);
#endif
    cs_grid_restrict_row_var(ctx, f, c, wr, mgd->rhs_vx[(level+1)*2]);

    cs_grid_get_info(c,
                     nullptr,
                     nullptr,
                     nullptr,
                     nullptr,
                     nullptr,
                     &n_rows,
                     &n_cols_ext,
                     nullptr,
                     &n_g_rows);

#if CS_PROFILING == CS_PROFILING_NVTX
    nvtxRangePop();
#endif

    f = c;
    _n_rows = n_rows*db_size;

    t0 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[4]), &t1, &t0);
    lv_info->n_calls[4] += 1;

  } /* End of loop on levels (descent) */

  if (end_cycle == false) {

    /* Resolve coarsest level to convergence */
    /*---------------------------------------*/

#if CS_PROFILING == CS_PROFILING_NVTX
    sprintf(nvtx_name, "Coarse solve (%d)", level);
    nvtxRangePushA(nvtx_name);
#endif

    assert(level == coarsest_level);
    assert(c == mgd->grid_hierarchy[coarsest_level]);

    /* coarsest level == 0 should never happen, but we play it safe */
    rhs_lv = (level == 0) ?  rhs : mgd->rhs_vx[coarsest_level*2];
    vx_lv = mgd->rhs_vx[level*2 + 1];

    _matrix = cs_grid_get_matrix(c);

    size_t _aux_r_size = aux_size_h;
    cs_real_t *_aux_vectors = reinterpret_cast<cs_real_t *>(aux_vectors_h);
#if defined(HAVE_ACCEL)
    if (cs_matrix_get_alloc_mode(_matrix) > CS_ALLOC_HOST) {
      _aux_r_size = aux_size;
      _aux_vectors = reinterpret_cast<cs_real_t *>(aux_vectors);
    }
#endif

    cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[level*2]);

    _initial_residual = _residual;

    lv_info = mg->lv_info + level;

    t0 = cs_timer_time();

    c_cvg = mg_sles->solve_func(mg_sles->context,
                                lv_names[level*2],
                                _matrix,
                                verbosity - 3,
                                precision*mg->info.precision_mult[2],
                                r_norm_l,
                                &n_iter,
                                &_residual,
                                rhs_lv,
                                nullptr,
                                vx_lv,
                                _aux_r_size,
                                _aux_vectors);

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[1]), &t0, &t1);
    lv_info->n_calls[1] += 1;
    _lv_info_update_stage_iter(lv_info->n_it_solve, n_iter);

    if (mg_sles->solve_func == cs_sles_it_solve)
      _initial_residual = cs_sles_it_get_last_initial_residual
                            ((const cs_sles_it_t *)mg_sles->context);
    else
      _initial_residual = HUGE_VAL;

    *n_equiv_iter += n_iter * n_g_rows * denom_n_g_rows_0;

    if (verbosity > 1)
      _log_residual(mg, cycle_id, lv_names[level*2],
                    _matrix, rhs_lv, vx_lv);

    if (c_cvg < CS_SLES_BREAKDOWN)
      end_cycle = true;

#if CS_PROFILING == CS_PROFILING_NVTX
    nvtxRangePop();
#endif
  }

  if (end_cycle == false) {

    /* Ascent */
    /*--------*/

    for (level = coarsest_level - 1; level > -1; level--) {

      vx_lv = mgd->rhs_vx[level*2 + 1];

      lv_info = mg->lv_info + level;

      c = mgd->grid_hierarchy[level+1];
      f = mgd->grid_hierarchy[level];

      cs_grid_get_info(f,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr,
                       &n_rows,
                       &n_cols_ext,
                       nullptr,
                       &n_g_rows);

      /* Prolong correction */

      t0 = cs_timer_time();

      cs_real_t *restrict vx_lv1 = mgd->rhs_vx[(level+1)*2 + 1];

      _matrix = cs_grid_get_matrix(f);

#if defined(HAVE_ACCEL)
      cs_alloc_mode_t amode_p = amode_f;
      amode_f = cs_matrix_get_alloc_mode(_matrix);
      if (amode_f > CS_ALLOC_HOST) {
        ctx.set_use_gpu(true);
#if defined(HAVE_CUDA)
        cudaStream_t stream = cs_matrix_spmv_cuda_get_stream();
        if (stream != 0)
          ctx.set_cuda_stream(stream);
#endif
        if (amode_p == CS_ALLOC_HOST)
          cs_prefetch_h2d(vx_lv1, _n_rows*sizeof(cs_real_t));
      }
#endif

#if CS_PROFILING == CS_PROFILING_NVTX
      sprintf(nvtx_name, "Prolong %d to %d", level+1, level);
      nvtxRangePushA(nvtx_name);
#endif

      cs_grid_prolong_row_var(ctx,
                              c, f,
                              true,  /* increment */
                              vx_lv1, vx_lv);

#if CS_PROFILING == CS_PROFILING_NVTX
      nvtxRangePop();
#endif
      _n_rows = n_rows*db_size;

      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&(lv_info->t_tot[5]), &t0, &t1);
      lv_info->n_calls[5] += 1;

      rhs_lv = (level == 0) ? rhs : mgd->rhs_vx[level*2];

      size_t _aux_r_size = aux_size_h;
      cs_real_t *_aux_vectors = reinterpret_cast<cs_real_t *>(aux_vectors_h);
#if defined(HAVE_ACCEL)
      if (amode_f > CS_ALLOC_HOST) {
        _aux_r_size = aux_size;
        _aux_vectors = reinterpret_cast<cs_real_t *>(aux_vectors);
      }
#endif

      cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[level*2 + 1]);

#if CS_PROFILING == CS_PROFILING_NVTX
      sprintf(nvtx_name, "Ascent smoothe %d", level);
      nvtxRangePushA(nvtx_name);
#endif

      c_cvg = mg_sles->solve_func(mg_sles->context,
                                  lv_names[level*2+1],
                                  _matrix,
                                  verbosity - 4, /* verbosity */
                                  precision*mg->info.precision_mult[1],
                                  r_norm_l,
                                  &n_iter,
                                  &_residual,
                                  rhs_lv,
                                  vx_lv,
                                  vx_lv,
                                  _aux_r_size,
                                  _aux_vectors);
#if CS_PROFILING == CS_PROFILING_NVTX
      nvtxRangePop();
#endif

      t0 = cs_timer_time();
      cs_timer_counter_add_diff(&(lv_info->t_tot[3]), &t1, &t0);
      lv_info->n_calls[3] += 1;
      _lv_info_update_stage_iter(lv_info->n_it_as_smoothe, n_iter);

      if (mg_sles->solve_func == cs_sles_it_solve)
        _initial_residual = cs_sles_it_get_last_initial_residual
                              ((const cs_sles_it_t *)mg_sles->context);
      else
        _initial_residual = HUGE_VAL;

      if (level == 0)
        *residual = _residual;

      *n_equiv_iter += n_iter * n_g_rows * denom_n_g_rows_0;

      if (verbosity > 1)
        _log_residual(mg, cycle_id, lv_names[level*2+1],
                      _matrix, rhs_lv, vx_lv);

      if (c_cvg < CS_SLES_BREAKDOWN)
        break;

    } /* End loop on levels (ascent) */

  } /* End of tests on end_cycle */

  mgd->exit_level = level;
  mgd->exit_residual = _residual;
  if (level == 0)
    mgd->exit_initial_residual = *initial_residual;
  else
    mgd->exit_initial_residual = _initial_residual;
  mgd->exit_cycle_id = cycle_id;

#if CS_PROFILING == CS_PROFILING_NVTX
  if (call_id == 2) {
    nvtxMark("Exit multigrid cycle");
  }
#endif

  return CS_SLES_ITERATING;
}

/*----------------------------------------------------------------------------
 * Apply a multigrid K-cycle on the current level.
 *
 * This version of the multigrid cycle is a K-cycle, described by Y. Notay in
 * "An aggregation-based algebraic multigrid method", Electronic Transactions on
 * Numerical Analysis, 2010, vol 37, pages 123-146.
 *
 * Used as a preconditioner, it is not constant (sometimes there is only one
 * visit of the coarse level, sometimes 2) and therefore, it has to be used with
 * IPCG (or FCG) instead of CG. Descent and ascent smoothers still have to be
 * symmetric for the IPCG to converge.
 *
 * parameters:
 *   mg               <-- multigrid system
 *   level            <-- level on which we apply the cycle
 *   verbosity        <-- verbosity level
 *   cycle_id         <-- id of currect cycle
 *   n_equiv_iter     <-> equivalent number of iterations
 *   precision        <-- solver precision
 *   r_norm           <-- residual normalization
 *   initial_residual <-> initial residual
 *   residual         <-> residual
 *   rhs              <-- right hand side
 *   vx               --> system solution
 *   aux_size         <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors      --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence status
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_multigrid_k_cycle(cs_multigrid_t       *mg,
                   int                   level,
                   int                   verbosity,
                   int                   cycle_id,
                   int                  *n_equiv_iter,
                   double                precision,
                   double                r_norm,
                   double               *initial_residual,
                   double               *residual,
                   const cs_real_t      *rhs,
                   cs_real_t            *vx,
                   size_t                aux_size,
                   void                 *aux_vectors)
{
  cs_timer_t t0, t1;

  bool end_cycle = false;
  cs_sles_convergence_state_t c_cvg = CS_SLES_ITERATING, cvg = CS_SLES_ITERATING;
  int n_iter = 0;
  double _residual = -1.;
  cs_real_t r_norm_l = r_norm;

  size_t _aux_r_size = aux_size / sizeof(cs_real_t);
  cs_lnum_t f_n_rows = 0, c_n_rows = 0;
  cs_lnum_t f_n_cols_ext = 0, c_n_cols_ext = 0;
  cs_gnum_t f_n_g_rows = 0, c_n_g_rows = 0;

  cs_multigrid_setup_data_t *mgd = mg->setup_data;
  int coarsest_level = mgd->n_levels-1;

  const char **lv_names = mgd->lv_names;

  cs_lnum_t db_size = 1;
  cs_lnum_t eb_size = 1;

  /* Recursion threshold */
  const cs_real_t trsh = mg->k_cycle_threshold;

  /* RHS and unknowns arrays */
  cs_real_t *restrict vx_lv = nullptr;
  const cs_real_t *restrict rhs_lv = nullptr;
  const cs_matrix_t  *f_matrix = nullptr, *c_matrix = nullptr;

  /* Initialization of the computing tools */

  cs_grid_t *f = mgd->grid_hierarchy[level];
  cs_grid_t *c = mgd->grid_hierarchy[level+1];

  f_matrix = cs_grid_get_matrix(f);
  c_matrix = cs_grid_get_matrix(c);

  cs_grid_get_info(f, nullptr, nullptr, &db_size, &eb_size, nullptr, &f_n_rows,
                   &f_n_cols_ext, nullptr, &f_n_g_rows);
  cs_grid_get_info(c, nullptr, nullptr, &db_size, &eb_size, nullptr, &c_n_rows,
                   &c_n_cols_ext, nullptr, &c_n_g_rows);

  cs_lnum_t _f_n_rows = f_n_rows*db_size;
  cs_lnum_t _c_n_rows = c_n_rows*db_size;

  static double denom_n_g_rows_0 = 1;
  if (level == 0)
    denom_n_g_rows_0 = 1.0 / f_n_g_rows;

  rhs_lv = rhs;
  vx_lv = vx;

  cs_dispatch_context ctx_f, ctx_c;

#if defined(HAVE_ACCEL)
  bool f_use_gpu = (cs_matrix_get_alloc_mode(f_matrix) > CS_ALLOC_HOST)
    ? true : false;
  bool c_use_gpu = (cs_matrix_get_alloc_mode(c_matrix) > CS_ALLOC_HOST)
    ? true : false;
  ctx_f.set_use_gpu(f_use_gpu);
  ctx_c.set_use_gpu(c_use_gpu);
# if defined(HAVE_CUDA)
  cudaStream_t stream = cs_matrix_spmv_cuda_get_stream();
  if (stream != 0) {
    ctx_f.set_cuda_stream(stream);
    ctx_c.set_cuda_stream(stream);
  }
# endif
#endif

  /* Descent
     ------- */

  const int na = 10; /* number of array slots in mgd->rhs_vx */

  /* Arrays that are needed for both descent and ascent phases */
  cs_real_t *restrict rt_lv = mgd->rhs_vx[level*na];
  assert(rt_lv != nullptr);

  cs_multigrid_level_info_t  *lv_info = mg->lv_info + level;
  t0 = cs_timer_time();

  /* Pre-smoothing */

  cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[level*2]);

  c_cvg = mg_sles->solve_func(mg_sles->context,
                              lv_names[level*2],
                              f_matrix,
                              0, /* verbosity */
                              precision*mg->info.precision_mult[0],
                              r_norm_l,
                              &n_iter,
                              &_residual,
                              rhs_lv,
                              nullptr,
                              vx_lv,
                              _aux_r_size*sizeof(cs_real_t),
                              aux_vectors);

  if (initial_residual != nullptr && cycle_id == 1)
    *initial_residual = HUGE_VAL;

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(lv_info->t_tot[2]), &t0, &t1);
  lv_info->n_calls[2] += 1;
  _lv_info_update_stage_iter(lv_info->n_it_ds_smoothe, n_iter);
  *n_equiv_iter += n_iter * f_n_g_rows * denom_n_g_rows_0;

  if (verbosity > 0)
    _log_residual(mg, cycle_id, lv_names[level*2],
                  f_matrix, rhs_lv, vx_lv);

  /* Compute new residual */

#if defined(HAVE_ACCEL)
  if (f_use_gpu) {
    cs_matrix_vector_multiply_d(f_matrix, vx_lv, rt_lv);
  }
  else
#endif
    cs_matrix_vector_multiply(f_matrix, vx_lv, rt_lv);

  ctx_f.parallel_for(_f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
    rt_lv[ii] = rhs_lv[ii] - rt_lv[ii];
  });
  ctx_f.wait();

  /* Convergence test in beginning of cycle (fine mesh) */

  if ((level == 0 && mg->info.n_max_cycles > 1) || c_cvg < CS_SLES_BREAKDOWN) {

    if (c_cvg >= CS_SLES_BREAKDOWN)
      cvg = _convergence_test(mg,
                              lv_names[0],
                              f_n_rows*db_size,
                              mg->info.n_max_cycles,
                              cycle_id,
                              verbosity,
                              lv_info->n_it_ds_smoothe[3],
                              precision,
                              r_norm,
                              *initial_residual,
                              residual,
                              rt_lv);
    else
      cvg = c_cvg;

    /* If converged or cycle limit reached, break from descent loop */

    if (cvg != CS_SLES_ITERATING)
      end_cycle = true;
  }

  t0 = cs_timer_time();
  cs_timer_counter_add_diff(&(lv_info->t_tot[6]), &t1, &t0);
  lv_info->n_calls[6] += 1; /* Single update for this cycle */

  if (end_cycle) {
    return cvg;
  }

  /* Arrays needed for the descent phase
     (mapped to preallocated mgd->rhs_vx to avoid extra allocations) */

  cs_real_t *restrict vx_lv1 = mgd->rhs_vx[(level+1)*na + 4];
  cs_real_t *restrict rhs_lv1 = mgd->rhs_vx[(level+1)*na + 5];

  /* Restriction of rt_lv into the next level rhs */
  cs_grid_restrict_row_var(ctx_f, f, c, rt_lv, rhs_lv1);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(lv_info->t_tot[4]), &t0, &t1);
  lv_info->n_calls[4] += 1;

  /* Approximate solution to the coarser problem
     ------------------------------------------- */

  /* If the next level is the coarsest one, then we solve it directly */

  if (level == coarsest_level-1) {

    lv_info = mg->lv_info + coarsest_level;

    mg_sles = &(mgd->sles_hierarchy[coarsest_level*2]);

    cvg = mg_sles->solve_func(mg_sles->context,
                              lv_names[coarsest_level*2],
                              c_matrix,
                              verbosity - 2,
                              precision*mg->info.precision_mult[2],
                              r_norm_l,
                              &n_iter,
                              residual,
                              rhs_lv1,
                              nullptr,
                              vx_lv1,
                              _aux_r_size*sizeof(cs_real_t),
                              aux_vectors);

    t0 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[1]), &t1, &t0);
    lv_info->n_calls[1] += 1;
    _lv_info_update_stage_iter(lv_info->n_it_solve, n_iter);

    if (verbosity > 0)
      _log_residual(mg, cycle_id, lv_names[coarsest_level*2],
                    c_matrix, rhs_lv1, vx_lv1);

    *n_equiv_iter += n_iter * c_n_g_rows * denom_n_g_rows_0;

  }

  /* Otherwise, we apply the K-cycle recursively */

  else {

    c_cvg = _multigrid_k_cycle(mg,
                               level + 1,
                               verbosity,
                               cycle_id,
                               &n_iter,
                               precision,
                               r_norm,
                               nullptr,
                               residual,
                               rhs_lv1,
                               vx_lv1,
                               aux_size,
                               aux_vectors);

    t1 = cs_timer_time();

    /* Arrays needed for the first Krylov iteration inside the cycle */
    cs_real_t *restrict v_lv1 = mgd->rhs_vx[(level+1)*na + 6];
    cs_real_t *restrict rt_lv1 = mgd->rhs_vx[(level+1)*na + 7];

#if defined(HAVE_ACCEL)
    if (c_use_gpu) {
      cs_matrix_vector_multiply_d(c_matrix, vx_lv1, rt_lv1);
    }
    else
#endif
      cs_matrix_vector_multiply(c_matrix, vx_lv1, v_lv1);

    /* Coefficients for the Krylov iteration */

    cs_real_t rho1 = 0, alpha1 = 0;
    _dot_xy_yz(mg, _c_n_rows, v_lv1, vx_lv1, rhs_lv1, &rho1, &alpha1);

    cs_real_t ar1 = alpha1 / rho1;

    /* New residual */
    for (cs_lnum_t i = 0; i < _c_n_rows; i++)
      rt_lv1[i] = rhs_lv1[i] - ar1 * v_lv1[i];

    cs_real_t  rt_lv1_norm = 1.0, r_lv1_norm = 0.0;

    if (trsh > 0)
      _dot_xx_yy(mg, _c_n_rows, rt_lv1, rhs_lv1, &rt_lv1_norm, &r_lv1_norm);

    /* Free (unmap) arrays that were needed only for the descent phase */
    rhs_lv1 = nullptr;

    /* Test for the second coarse resolution */
    if (rt_lv1_norm < trsh * trsh * r_lv1_norm) {
      ctx_c.parallel_for(_c_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
        vx_lv1[i] = ar1 * vx_lv1[i];
      });
      ctx_c.wait();
    }
    else {

      /* Arrays for the (optional) second Krylov iteration */
      cs_real_t *restrict vx2_lv1 = mgd->rhs_vx[(level+1)*na + 8];
      cs_real_t *restrict w_lv1 = mgd->rhs_vx[(level+1)*na + 9];

      ctx_c.parallel_for(_c_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
        w_lv1[i] = 0.0;
      });
      ctx_c.wait();

      t0 = cs_timer_time();
      cs_timer_counter_add_diff(&(lv_info->t_tot[6]), &t1, &t0);

      c_cvg = _multigrid_k_cycle(mg,
                                 level + 1,
                                 verbosity,
                                 cycle_id,
                                 &n_iter,
                                 precision,
                                 r_norm,
                                 initial_residual,
                                 residual,
                                 rt_lv1,
                                 vx2_lv1,
                                 aux_size,
                                 aux_vectors);

      t1 = cs_timer_time();

#if defined(HAVE_ACCEL)
      if (c_use_gpu) {
        cs_matrix_vector_multiply_d(c_matrix, vx2_lv1, w_lv1);
      }
      else
#endif
        cs_matrix_vector_multiply(c_matrix, vx2_lv1, w_lv1);

      /* Krylov iteration */

      cs_real_t gamma = 0, beta = 0, alpha2 = 0;

      _dot_xu_xv_xw(mg, _c_n_rows,
                    vx2_lv1, v_lv1, w_lv1, rt_lv1,
                    &gamma, &beta, &alpha2);

      cs_real_t rho2 = beta - (gamma * gamma) / rho1;
      cs_real_t ar2 = alpha2 / rho2;
      cs_real_t ar1_ar2 = ar1 - (gamma / rho1) * ar2;

      ctx_c.parallel_for(_c_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
        vx_lv1[i] = ar1_ar2 * vx_lv1[i] + ar2 * vx2_lv1[i];
      });
      ctx_c.wait();

      vx2_lv1 = nullptr;
      w_lv1 = nullptr;
    }

    t0 = cs_timer_time();
    cs_timer_counter_add_diff(&(lv_info->t_tot[6]), &t1, &t0);

    v_lv1 = nullptr;
    rt_lv1 = nullptr;
  }

  /* Ascent
     ------ */

  lv_info = mg->lv_info + level;
  t0 = cs_timer_time();

  /* Arrays that are needed for the ascent phase */
  cs_real_t *restrict rb_lv = mgd->rhs_vx[level*na + 1];
  cs_real_t *restrict z1_lv = mgd->rhs_vx[level*na + 2];
  cs_real_t *restrict z2_lv = mgd->rhs_vx[level*na + 3];

  ctx_f.parallel_for(_f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
    z2_lv[ii] = 0.0;
  });
  ctx_f.wait();

  /* Prolongation */
  cs_grid_prolong_row_var(ctx_f, c, f, false, vx_lv1, z1_lv);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(lv_info->t_tot[5]), &t0, &t1);
  lv_info->n_calls[5] += 1;

  /* New residual */
#if defined(HAVE_ACCEL)
  if (f_use_gpu) {
    cs_matrix_vector_multiply_d(f_matrix, z1_lv, rb_lv);
  }
  else
#endif
    cs_matrix_vector_multiply(f_matrix, z1_lv, rb_lv);

  ctx_f.parallel_for(_f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
    rb_lv[ii] = rt_lv[ii] - rb_lv[ii];
  });
  ctx_f.wait();

  t0 = cs_timer_time();
  cs_timer_counter_add_diff(&(lv_info->t_tot[5]), &t1, &t0);
  t1 = t0;

  /* Post-smoothing */

  mg_sles = &(mgd->sles_hierarchy[level*2+1]);

  cvg = mg_sles->solve_func(mg_sles->context,
                            lv_names[level*2+1],
                            f_matrix,
                            0, /* verbosity */
                            precision*mg->info.precision_mult[1],
                            r_norm_l,
                            &n_iter,
                            &_residual,
                            rb_lv,
                            z2_lv,
                            z2_lv,
                            _aux_r_size*sizeof(cs_real_t),
                            aux_vectors);

  t0 = cs_timer_time();
  cs_timer_counter_add_diff(&(lv_info->t_tot[3]), &t1, &t0);
  lv_info->n_calls[3] += 1;
  _lv_info_update_stage_iter(lv_info->n_it_as_smoothe, n_iter);
  *n_equiv_iter += n_iter * f_n_g_rows * denom_n_g_rows_0;

  if (verbosity > 0)
    _log_residual(mg, cycle_id, lv_names[level*2 + 1],
                  c_matrix, rb_lv, z2_lv);

  ctx_f.parallel_for(_f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
    vx_lv[ii] += z1_lv[ii] + z2_lv[ii];
  });
  ctx_f.wait();

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(lv_info->t_tot[6]), &t0, &t1);

  /* Free/unmap arrays of the ascent phase */
  rb_lv = nullptr;
  z1_lv = nullptr;
  z2_lv = nullptr;

  /* Free/unmap arrays that were needed for both phases */
  rt_lv = nullptr;
  vx_lv1 = nullptr;

  return cvg;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize multigrid solver API.
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_initialize(void)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize multigrid solver API.
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_finalize(void)
{
  cs_grid_finalize();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set maximum grid level which should run on device (i.e. GPU).
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_max_grid_level_for_device(int  level)
{
  _grid_max_level_for_device = level;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and associate a multigrid sparse linear system solver
 *        for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * \ref cs_sles_define and \ref cs_multigrid_create.
 *
 * Note that this function returns a pointer directly to the multigrid solver
 * management structure. This may be used to set further options, for
 * example calling \ref cs_multigrid_set_coarsening_options and
 * \ref cs_multigrid_set_solver_options.
 * If needed, \ref cs_sles_find may be used to obtain a pointer to the
 * matching \ref cs_sles_t container.
 *
 * \param[in]  f_id     associated field id, or < 0
 * \param[in]  name     associated name if f_id < 0, or nullptr
 * \param[in]  mg_type  type of multigrid algorithm to use
 *
 * \return  pointer to new multigrid info and context
 */
/*----------------------------------------------------------------------------*/

cs_multigrid_t *
cs_multigrid_define(int                   f_id,
                    const char           *name,
                    cs_multigrid_type_t   mg_type)
{
  cs_multigrid_t *
    mg = cs_multigrid_create(mg_type);

  cs_sles_t *sc = cs_sles_define(f_id,
                                 name,
                                 mg,
                                 "cs_multigrid_t",
                                 cs_multigrid_setup,
                                 cs_multigrid_solve,
                                 cs_multigrid_free,
                                 cs_multigrid_log,
                                 cs_multigrid_copy,
                                 cs_multigrid_destroy);

  cs_sles_set_error_handler(sc,
                            cs_multigrid_error_post_and_abort);

  return mg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create multigrid linear system solver info and context.
 *
 * The multigrid variant is an ACM (Additive Corrective Multigrid) method.
 *
 * \param[in]  mg_type  type of multigrid algorithm to use
 *
 * \return  pointer to new multigrid info and context
 */
/*----------------------------------------------------------------------------*/

cs_multigrid_t *
cs_multigrid_create(cs_multigrid_type_t  mg_type)
{
  int ii;
  cs_multigrid_t *mg;

  /* Increment number of setups */

  CS_MALLOC(mg, 1, cs_multigrid_t);

  mg->type = mg_type;
  mg->subtype = CS_MULTIGRID_MAIN;

  mg->caller_n_ranks = cs_glob_n_ranks;

#if defined(HAVE_MPI)
  mg->comm = cs_glob_mpi_comm;
  mg->caller_comm = cs_glob_mpi_comm;
  if (mg->caller_n_ranks < 2) {
    mg->comm = MPI_COMM_NULL;
  }
#endif

  mg->merge_stride = 1;
  mg->merge_bottom_n_max_ranks = cs_glob_n_ranks;
  mg->merge_bottom_max_row_factor = 1.0;

  mg->merge_mean_threshold = 300;
  mg->merge_glob_threshold = 500;

  mg->f_settings_threshold = -1;
  for (int i = 0; i < 2; i++) {
    mg->aggregation_limit[i] = 3;
    mg->coarsening_type[i] = CS_GRID_COARSENING_DEFAULT;
  }
  mg->n_levels_max = 25;
  mg->n_g_rows_min = 30;

  mg->post_row_max = 0;

  mg->p0p1_relax = 0.;
  mg->k_cycle_threshold = 0;

  _multigrid_info_init(&(mg->info));
  for (int i = 0; i < 3; i++)
    mg->lv_mg[i] = nullptr;
  mg->p_mg = nullptr;

  if (mg->type == CS_MULTIGRID_V_CYCLE) {
    mg->p0p1_relax = 0.95;
  }
  else if (mg->type == CS_MULTIGRID_K_CYCLE) {
    for (int i = 0; i < 2; i++) {
      mg->coarsening_type[i] = CS_GRID_COARSENING_SPD_PW;
      mg->aggregation_limit[i] = 4;
    }
    mg->n_levels_max = 10;
    mg->n_g_rows_min = 256;
    mg->k_cycle_threshold = 0.25;
  }
  else if (mg->type == CS_MULTIGRID_K_CYCLE_HPC) {
    for (int i = 0; i < 2; i++) {
      mg->coarsening_type[i] = CS_GRID_COARSENING_SPD_PW;
      mg->aggregation_limit[i] = 8;
    }
    mg->n_levels_max = 4;
    mg->k_cycle_threshold = -1;
    mg->lv_mg[2] = _multigrid_create_k_cycle_bottom(mg);
  }

  mg->pc_precision = 0.0;
  mg->pc_r_norm = 0.0;

  mg->n_levels_post = 0;

  mg->setup_data = nullptr;

  CS_MALLOC(mg->lv_info, mg->n_levels_max, cs_multigrid_level_info_t);

  for (ii = 0; ii < mg->n_levels_max; ii++)
    _multigrid_level_info_init(mg->lv_info + ii);

  mg->post_row_num = nullptr;
  mg->post_row_rank = nullptr;
  mg->post_name = nullptr;

  mg->cycle_plot = nullptr;
  mg->plot_time_stamp = -1;

  if (mg_type == CS_MULTIGRID_V_CYCLE)
    cs_multigrid_set_solver_options
      (mg,
       CS_SLES_PCG, CS_SLES_PCG, CS_SLES_PCG,
       100, /* n max cycles */
       2,   /* n max iter for descent */
       10,  /* n max iter for ascent */
       500,
       0, 0, 0,  /* precond degree */
       1, 1, 1); /* precision multiplier */

  else if (   mg->type >= CS_MULTIGRID_K_CYCLE
           && mg->type <= CS_MULTIGRID_K_CYCLE_HPC)
    cs_multigrid_set_solver_options
      (mg,
       CS_SLES_P_SYM_GAUSS_SEIDEL,
       CS_SLES_P_SYM_GAUSS_SEIDEL,
       CS_SLES_PCG,
       100,  /* n max cycles */
       1,    /* n max iter for descent */
       1,    /* n max iter for ascent */
       50,   /* n max iter for coarse solve */
       -1, -1, 0,   /* precond degree */
       -1, -1, 1);  /* precision multiplier */

  return mg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy multigrid linear system solver info and context.
 *
 * \param[in, out]  context  pointer to multigrid linear solver info
 *                           (actual type: cs_multigrid_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_destroy(void  **context)
{
  cs_multigrid_t *mg = (cs_multigrid_t *)(*context);

  if (mg == nullptr)
    return;

  CS_FREE(mg->lv_info);

  if (mg->post_row_num != nullptr) {
    int n_max_post_levels = (int)(mg->info.n_levels[2]) - 1;
    for (int i = 0; i < n_max_post_levels; i++)
      if (mg->post_row_num[i] != nullptr)
        CS_FREE(mg->post_row_num[i]);
    CS_FREE(mg->post_row_num);
  }

  if (mg->post_row_rank != nullptr) {
    int n_max_post_levels = (int)(mg->info.n_levels[2]) - 1;
    for (int i = 0; i < n_max_post_levels; i++)
      if (mg->post_row_rank[i] != nullptr)
        CS_FREE(mg->post_row_rank[i]);
    CS_FREE(mg->post_row_rank);
  }

  CS_FREE(mg->post_name);

  if (mg->cycle_plot != nullptr)
    cs_time_plot_finalize(&(mg->cycle_plot));

  for (int i = 0; i < 3; i++) {
    if (mg->lv_mg[i] != nullptr)
      cs_multigrid_destroy((void **)(&(mg->lv_mg[i])));
  }

  CS_FREE(mg);
  *context = (void *)mg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create multigrid sparse linear system solver info and context
 *        based on existing info and context.
 *
 * \param[in]  context  pointer to reference info and context
 *                      (actual type: cs_multigrid_t  *)
 *
 * \return  pointer to newly created solver info object
 *          (actual type: cs_multigrid_t  *)
 */
/*----------------------------------------------------------------------------*/

void *
cs_multigrid_copy(const void  *context)
{
  cs_multigrid_t *d = nullptr;

  if (context != nullptr) {

    const cs_multigrid_t *c = (const cs_multigrid_t *)context;
    d = cs_multigrid_create(c->type);

    /* Beginning of cs_multigrid_info_t contains settings, the rest logging */
    memcpy(&(d->info), &(c->info),
           offsetof(cs_multigrid_info_t, n_calls));

    /* Same here: settings at beginningof structure */
    memcpy(d, c, offsetof(cs_multigrid_t, n_levels_post));

  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log multigrid solver info.
 *
 * \param[in]  context   pointer to iterative solver info and context
 *                       (actual type: cs_multigrid_t  *)
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_log(const void  *context,
                 cs_log_t     log_type)
{
  const cs_multigrid_t *mg = (const cs_multigrid_t *)context;

  if (log_type == CS_LOG_SETUP)
    _multigrid_setup_log(mg);

  else if (log_type == CS_LOG_PERFORMANCE)
    _multigrid_performance_log(mg);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set multigrid coarsening parameters.
 *
 * \param[in, out]  mg                 pointer to multigrid info and context
 * \param[in]       aggregation_limit  maximum allowed fine rows
 *                                     per coarse row
 * \param[in]       coarsening_type    coarsening type;
 *                                     see \ref cs_grid_coarsening_t
 * \param[in]       n_max_levels       maximum number of grid levels
 * \param[in]       min_g_rows         global number of rows on coarse grids
 *                                     under which no coarsening occurs
 * \param[in]       p0p1_relax         p0/p1 relaxation_parameter
 * \param[in]       postprocess_block_size  if > 0, postprocess coarsening
 *                                     (uses coarse row numbers
 *                                      modulo this value)
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_coarsening_options
  (cs_multigrid_t       *mg,
   int                   aggregation_limit,
   cs_grid_coarsening_t  coarsening_type,
   int                   n_max_levels,
   cs_gnum_t             min_g_rows,
   double                p0p1_relax,
   int                   postprocess_block_size)
{
  if (mg == nullptr)
    return;

  mg->aggregation_limit[1] = aggregation_limit;
  mg->coarsening_type[1] = coarsening_type;
  mg->n_levels_max = n_max_levels;
  mg->n_g_rows_min = min_g_rows;

  mg->post_row_max = postprocess_block_size;

  mg->p0p1_relax = p0p1_relax;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set specific multigrid coarsening parameters for fine grids.
 *
 * \param[in, out]  mg                    pointer to multigrid info and context
 * \param[in]       f_settings_threshold  grids of this level or higher use
 *                                        standard (coarse grid) settings.
 * \param[in]       aggregation_limit     maximum allowed fine rows
 *                                        per coarse row
 * \param[in]       coarsening_type       coarsening type;
 *                                        see \ref cs_grid_coarsening_t
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_coarsening_options_fine_grid
  (cs_multigrid_t       *mg,
   int                   f_settings_threshold,
   int                   aggregation_limit,
   cs_grid_coarsening_t  coarsening_type)
{
  if (mg == nullptr)
    return;

  mg->f_settings_threshold = f_settings_threshold;
  mg->aggregation_limit[0] = aggregation_limit;
  mg->coarsening_type[0] = coarsening_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set multigrid parameters for associated iterative solvers.
 *
 * On a GPU, some parameters may be replaced by the closest GPU equivalents.
 * For finer control, use \ref cs_multigrid_set_solver_options_d after
 * calling this function to modify parameters for solvers running on device.
 *
 * \param[in, out]  mg                      pointer to multigrid info
 *                                          and context
 * \param[in]       descent_smoother_type   type of smoother for descent
 * \param[in]       ascent_smoother_type    type of smoother for ascent
 * \param[in]       coarse_solver_type      type of solver for coarsest grid
 * \param[in]       n_max_cycles            maximum number of cycles
 * \param[in]       n_max_iter_descent      maximum iterations
 *                                          per descent smoothing
 * \param[in]       n_max_iter_ascent       maximum iterations
 *                                          per ascent smoothing
 * \param[in]       n_max_iter_coarse       maximum iterations
 *                                          per coarsest solution
 * \param[in]       poly_degree_descent     preconditioning polynomial degree
 *                                          for descent phases (0: diagonal)
 * \param[in]       poly_degree_ascent      preconditioning polynomial degree
 *                                          for ascent phases (0: diagonal)
 * \param[in]       poly_degree_coarse      preconditioning polynomial degree
 *                                          for coarse solver (0: diagonal)
 * \param[in]       precision_mult_descent  precision multiplier
 *                                          for descent smoothers (levels >= 1)
 * \param[in]       precision_mult_ascent   precision multiplier
 *                                          for ascent smoothers
 * \param[in]       precision_mult_coarse   precision multiplier
 *                                          for coarsest grid
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_solver_options(cs_multigrid_t     *mg,
                                cs_sles_it_type_t   descent_smoother_type,
                                cs_sles_it_type_t   ascent_smoother_type,
                                cs_sles_it_type_t   coarse_solver_type,
                                int                 n_max_cycles,
                                int                 n_max_iter_descent,
                                int                 n_max_iter_ascent,
                                int                 n_max_iter_coarse,
                                int                 poly_degree_descent,
                                int                 poly_degree_ascent,
                                int                 poly_degree_coarse,
                                double              precision_mult_descent,
                                double              precision_mult_ascent,
                                double              precision_mult_coarse)
{
  if (mg == nullptr)
    return;

  cs_multigrid_info_t  *info = &(mg->info);

  info->type[0] = descent_smoother_type;
  info->type[1] = ascent_smoother_type;
  info->type[2] = coarse_solver_type;

  info->n_max_cycles = n_max_cycles;

  info->n_max_iter[0] = n_max_iter_descent;
  info->n_max_iter[1] = n_max_iter_ascent;
  info->n_max_iter[2] = n_max_iter_coarse;

  info->poly_degree[0] = poly_degree_descent;
  info->poly_degree[1] = poly_degree_ascent;
  info->poly_degree[2] = poly_degree_coarse;

  info->precision_mult[0] = precision_mult_descent;
  info->precision_mult[1] = precision_mult_ascent;
  info->precision_mult[2] = precision_mult_coarse;

  if (cs_get_device_id() > -1 || _force_gpu_settings) {
    for (int i = 0; i < 3; i++) {
      info->type[i+3] = info->type[i];
      info->n_max_iter[i+3] = info->n_max_iter[i];
      info->poly_degree[i+3] = info->poly_degree[i];
      switch (info->type[i+3]) {
      case CS_SLES_P_GAUSS_SEIDEL:
        info->type[i+3] = CS_SLES_JACOBI;
        info->n_max_iter[i+3] *= 2;
        break;
      case CS_SLES_P_SYM_GAUSS_SEIDEL:
        info->type[i+3] = CS_SLES_JACOBI;
        if (i == 0 && info->n_max_iter[i+3] == 1)
          info->n_max_iter[i+3] = 2;
        else
          info->n_max_iter[i+3] *= 3;
        break;
      case CS_SLES_TS_F_GAUSS_SEIDEL:
        info->type[i+3] = CS_SLES_JACOBI;
        info->n_max_iter[i+3] = 2;
        break;
      case CS_SLES_TS_B_GAUSS_SEIDEL:
        info->type[i+3] = CS_SLES_JACOBI;
        info->n_max_iter[i+3] = 2;
        break;
      case CS_SLES_PCG:
        info->type[i+3] = CS_SLES_FCG;
        break;
      case CS_SLES_BICGSTAB:
        [[fallthrough]];
      case CS_SLES_BICGSTAB2:
        [[fallthrough]];
      case CS_SLES_GMRES:
        [[fallthrough]];
      case CS_SLES_PCR3:
        info->type[i+3] = CS_SLES_GCR;
        break;
      default:
        break;
      }
    }
  }

  for (int i = 0; i < 6; i++) {
    switch(info->type[i]) {
    case CS_SLES_JACOBI:
    case CS_SLES_P_GAUSS_SEIDEL:
    case CS_SLES_P_SYM_GAUSS_SEIDEL:
      info->poly_degree[i] = -1;
      break;
    default:
      break;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set multigrid device solver parameters for associated
 *        iterative solvers.
 *
 * \param[in, out]  mg                      pointer to multigrid info
 *                                          and context
 * \param[in]       descent_smoother_type   type of smoother for descent
 * \param[in]       ascent_smoother_type    type of smoother for ascent
 * \param[in]       coarse_solver_type      type of solver for coarsest grid
 * \param[in]       n_max_iter_descent      maximum iterations
 *                                          per descent smoothing
 * \param[in]       n_max_iter_ascent       maximum iterations
 *                                          per ascent smoothing
 * \param[in]       n_max_iter_coarse       maximum iterations
 *                                          per coarsest solution
 * \param[in]       poly_degree_descent     preconditioning polynomial degree
 *                                          for descent phases (0: diagonal)
 * \param[in]       poly_degree_ascent      preconditioning polynomial degree
 *                                          for ascent phases (0: diagonal)
 * \param[in]       poly_degree_coarse      preconditioning polynomial degree
 *                                          for coarse solver (0: diagonal)
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_solver_options_d(cs_multigrid_t     *mg,
                                  cs_sles_it_type_t   descent_smoother_type,
                                  cs_sles_it_type_t   ascent_smoother_type,
                                  cs_sles_it_type_t   coarse_solver_type,
                                  int                 n_max_iter_descent,
                                  int                 n_max_iter_ascent,
                                  int                 n_max_iter_coarse,
                                  int                 poly_degree_descent,
                                  int                 poly_degree_ascent,
                                  int                 poly_degree_coarse)
{
  if (mg == nullptr || (cs_get_device_id() < 0 && _force_gpu_settings == false))
    return;

  cs_multigrid_info_t  *info = &(mg->info);

  info->type[3] = descent_smoother_type;
  info->type[4] = ascent_smoother_type;
  info->type[5] = coarse_solver_type;

  info->n_max_iter[3] = n_max_iter_descent;
  info->n_max_iter[4] = n_max_iter_ascent;
  info->n_max_iter[5] = n_max_iter_coarse;

  info->poly_degree[3] = poly_degree_descent;
  info->poly_degree[4] = poly_degree_ascent;
  info->poly_degree[5] = poly_degree_coarse;

  for (int i = 0; i < 3; i++) {
    switch(info->type[i+3]) {
    case CS_SLES_JACOBI:
    case CS_SLES_P_GAUSS_SEIDEL:
    case CS_SLES_P_SYM_GAUSS_SEIDEL:
      info->poly_degree[i+3] = -1;
      break;
    default:
      break;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the max. number of cycles for a multigrid
 *
 * \param[in, out]  mg              pointer to multigrid info and context
 * \param[in]       n_max_cycles    maximum number of cycles
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_max_cycles(cs_multigrid_t     *mg,
                            int                 n_max_cycles)
{
  if (mg == nullptr)
    return;

  cs_multigrid_info_t  *info = &(mg->info);

  info->n_max_cycles = n_max_cycles;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if a multigrid solver requires an MSR matrix input.
 *
 * \param[in]  mg  pointer to multigrid info and context
 *
 * \return   true if MSR is needed, false otherwise.
 */
/*----------------------------------------------------------------------------*/

bool
cs_multigrid_need_msr(const cs_multigrid_t  *mg)
{
  bool retval = false;

  if (mg != nullptr) {
    const cs_multigrid_info_t  *info = &(mg->info);
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        int k = i + j*3;
        cs_sles_it_type_t fs_type = info->type[k];
        if (   fs_type >= CS_SLES_P_GAUSS_SEIDEL
            && fs_type <= CS_SLES_TS_B_GAUSS_SEIDEL)
          retval = true;
      }
      int fg_i = (mg->f_settings_threshold < 1) ? 1 : 0;
      if (mg->coarsening_type[fg_i] == CS_GRID_COARSENING_SPD_PW)
        retval = true;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup multigrid sparse linear equation solver.
 *
 * \param[in, out]  context    pointer to multigrid solver info and context
 *                             (actual type: cs_multigrid_t  *)
 * \param[in]       name       pointer to name of linear system
 * \param[in]       a          associated matrix
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_setup(void               *context,
                   const char         *name,
                   const cs_matrix_t  *a,
                   int                 verbosity)
{
  cs_multigrid_setup_conv_diff(context,
                               name,
                               a,
                               false,
                               verbosity);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup multigrid sparse linear equation solver.
 *
 * \param[in, out]  context    pointer to multigrid solver info and context
 *                             (actual type: cs_multigrid_t  *)
 * \param[in]       name       pointer to name of linear system
 * \param[in]       a          associated matrix
 * \param[in]       conv_diff  convection-diffusion mode
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_setup_conv_diff(void               *context,
                             const char         *name,
                             const cs_matrix_t  *a,
                             bool                conv_diff,
                             int                 verbosity)

{
  cs_multigrid_t  *mg = (cs_multigrid_t *)context;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Destroy previous hierarchy if necessary */

  if (mg->setup_data != nullptr)
    cs_multigrid_free(mg);

  /* Initialization */

  cs_timer_t t0 = cs_timer_time();

  if (verbosity > 1) {
    if (!(mg->info.is_pc)) {
      bft_printf(_("\n Setup of solver for linear system \"%s\"\n"),
                 name);
      cs_matrix_log_info(a, verbosity);
    }
    bft_printf(_("\n Construction of grid hierarchy for \"%s\"\n"),
               name);
  }

  /* Build coarse grids hierarchy */
  /*------------------------------*/

  const cs_lnum_t diag_block_size = cs_matrix_get_diag_block_size(a);
  const cs_lnum_t extra_diag_block_size
    = cs_matrix_get_extra_diag_block_size(a);

  cs_grid_t *f
    = cs_grid_create_from_shared(mesh->n_i_faces,
                                 diag_block_size,
                                 extra_diag_block_size,
                                 mesh->i_face_cells,
                                 a,
                                 conv_diff);

  cs_multigrid_level_info_t *mg_lv_info = mg->lv_info;

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg_lv_info->t_tot[0]), &t0, &t1);

  _setup_hierarchy(mg, name, mesh, f, verbosity); /* Assign to and build
                                                     hierarchy */

  /* Update timers */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg->info.t_tot[0]), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call multigrid sparse linear equation solver.
 *
 * \param[in, out]  context        pointer to multigrid solver info and context
 *                                 (actual type: cs_multigrid_t  *)
 * \param[in]       name           pointer to name of linear system
 * \param[in]       a              matrix
 * \param[in]       verbosity      associated verbosity
 * \param[in]       precision      solver precision
 * \param[in]       r_norm         residual normalization
 * \param[out]      n_iter         number of "equivalent" iterations
 * \param[out]      residual       residual
 * \param[in]       rhs            right hand side
 * \param[in]       vx_ini         initial system solution
 *                                 (vx if nonzero, nullptr if zero)
 * \param[in, out]  vx             system solution
 * \param[in]       aux_size       size of aux_vectors (in bytes)
 * \param           aux_vectors    optional working area
 *                                 (internal allocation if nullptr)
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_multigrid_solve(void                *context,
                   const char          *name,
                   const cs_matrix_t   *a,
                   int                  verbosity,
                   double               precision,
                   double               r_norm,
                   int                 *n_iter,
                   double              *residual,
                   const cs_real_t     *rhs,
                   cs_real_t           *vx_ini,
                   cs_real_t           *vx,
                   size_t               aux_size,
                   void                *aux_vectors)
{
  CS_PROFILE_FUNC_RANGE();

  cs_timer_t t0, t1;
  t0 = cs_timer_time();

  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  cs_multigrid_t *mg = (cs_multigrid_t *)context;
  cs_multigrid_info_t *mg_info = &(mg->info);

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);

  /* Initialize number of equivalent iterations and residual,
     check for immediate return,
     solve sparse linear system using multigrid algorithm. */

  *n_iter = 0;
  unsigned n_cycles = 0;

  if (mg->setup_data == nullptr) {
    /* Stop solve timer to switch to setup timer */
    t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(mg->info.t_tot[1]), &t0, &t1);

    /* Setup grid hierarchy */
    cs_multigrid_setup(context, name, a, verbosity);

    /* Restart solve timer */
    t0 = cs_timer_time();
  }

  if (mg_info->is_pc == false && verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" RHS norm:          %11.4e\n\n"), r_norm);

  /* Buffer size sufficient to avoid local reallocation for most solvers */
  size_t  _aux_size_d = 0;
  size_t  _aux_size_h = 0;

  const cs_multigrid_setup_data_t *mgd = mg->setup_data;
  for (int level = 0; level < (int)(mgd->n_levels); level++) {
    const cs_grid_t *g = mgd->grid_hierarchy[level];
    const cs_matrix_t *m = cs_grid_get_matrix(g);
    cs_alloc_mode_t amode_l = cs_matrix_get_alloc_mode(m);
    size_t wr_l =   cs_grid_get_n_cols_max(mgd->grid_hierarchy[level])
                  * db_size;
    wr_l = CS_SIMD_SIZE(wr_l);
    if (amode_l > CS_ALLOC_HOST) {
      if (wr_l > _aux_size_d)
        _aux_size_d = wr_l;
    }
    else if (wr_l > _aux_size_h)
      _aux_size_h = wr_l;
  }

  /* Buffer size sufficient to avoid local reallocation for most solvers */
  _aux_size_d = 6*_aux_size_d*sizeof(cs_real_t);
  _aux_size_h = 6*_aux_size_h*sizeof(cs_real_t);

  size_t _aux_size = aux_size;
  void *_aux_buf = nullptr;
  void *_aux_buf_h = nullptr;

  cs_alloc_mode_t amode_aux = cs_check_device_ptr(aux_vectors);
  if (amode_aux > CS_ALLOC_HOST) {
    if (_aux_size_d > aux_size) {
      CS_MALLOC_HD(_aux_buf, _aux_size_d,
                   unsigned char,
                   CS_ALLOC_HOST_DEVICE_SHARED);
      _aux_size = _aux_size_d;
    }
    else
      _aux_buf = aux_vectors;
    CS_MALLOC(_aux_buf_h, _aux_size_h, unsigned char);
  }
  else if (_aux_size_d > 0) {
    CS_MALLOC_HD(_aux_buf, _aux_size_d,
                 unsigned char,
                 CS_ALLOC_HOST_DEVICE_SHARED);
    _aux_size = _aux_size_d;
    if (_aux_size_h > aux_size)
      CS_MALLOC(_aux_buf_h, _aux_size_h, unsigned char);
    else {
      _aux_buf_h = aux_vectors;
      _aux_size_h = aux_size;
    }
  }
  else {
    if (_aux_size_h > aux_size) {
      CS_MALLOC(_aux_buf, _aux_size_h, unsigned char);
      _aux_size = _aux_size_h;
    }
    else {
      _aux_buf = aux_vectors;
      _aux_size_h = aux_size;
    }
    _aux_buf_h = _aux_buf;
  }

  if (verbosity == 2) /* More detailed headers later if > 2 */
    bft_printf(_("Multigrid [%s]:\n"), name);

  /* Initial residual should be improved, but this is consistent
     with the legacy (by increment) case */

  double initial_residual = -1;

  *residual = initial_residual; /* not known yet, so be safe */

  /* Cycle to solution */

  if (mg->type == CS_MULTIGRID_V_CYCLE) {
    if (mg->info.n_max_cycles == 1) {
      cvg = _multigrid_v_cycle_pc(mg,
                                  verbosity,
                                  n_iter,
                                  precision,
                                  r_norm,
                                  &initial_residual,
                                  residual,
                                  rhs,
                                  vx_ini,
                                  vx,
                                  _aux_size,
                                  _aux_size_h,
                                  _aux_buf,
                                  _aux_buf_h);
    }
    else {
      for (n_cycles = 0; cvg == CS_SLES_ITERATING; n_cycles++) {
        int cycle_id = n_cycles+1;
        cs_real_t *_vx_ini = (n_cycles > 0) ? vx : vx_ini;
        cvg = _multigrid_v_cycle(mg,
                                 verbosity,
                                 cycle_id,
                                 n_iter,
                                 precision,
                                 r_norm,
                                 &initial_residual,
                                 residual,
                                 rhs,
                                 _vx_ini,
                                 vx,
                                 _aux_size,
                                 _aux_buf);
      }
    }
  }
  else if (   mg->type >= CS_MULTIGRID_K_CYCLE
           && mg->type <= CS_MULTIGRID_K_CYCLE_HPC) {
    for (n_cycles = 0;
         n_cycles < (unsigned)mg->info.n_max_cycles && cvg == CS_SLES_ITERATING;
         n_cycles++) {
      int cycle_id = n_cycles+1;
      cvg = _multigrid_k_cycle(mg,
                               0,
                               verbosity,
                               cycle_id,
                               n_iter,
                               precision,
                               r_norm,
                               &initial_residual,
                               residual,
                               rhs,
                               vx,
                               _aux_size,
                               _aux_buf);
    }
  }

  if (_aux_buf_h == _aux_buf) {
    if (_aux_buf != aux_vectors) {
      CS_FREE(_aux_buf);
      _aux_buf_h = _aux_buf;
    }
  }
  else {
    if (_aux_buf != aux_vectors)
      CS_FREE_HD(_aux_buf);
    if (_aux_buf_h != aux_vectors)
      CS_FREE(_aux_buf_h);
  }

  if (verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "\n");

  /* Update statistics */

  t1 = cs_timer_time();

  /* Update stats on number of iterations (last, min, max, total) */

  mg_info->n_cycles[2] += n_cycles;

  if (mg_info->n_calls[1] > 0) {
    if (mg_info->n_cycles[0] > n_cycles)
      mg_info->n_cycles[0] = n_cycles;
    if (mg_info->n_cycles[1] < n_cycles)
      mg_info->n_cycles[1] = n_cycles;
  }
  else {
    mg_info->n_cycles[0] = n_cycles;
    mg_info->n_cycles[1] = n_cycles;
  }

  /* Update number of resolutions and timing data */

  mg_info->n_calls[1] += 1;
  cs_timer_counter_add_diff(&(mg->info.t_tot[1]), &t0, &t1);

  return cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free multigrid sparse linear equation solver setup context.
 *
 * This function frees resolution-related data, incuding the current
 * grid hierarchy, but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.
 *
 * \param[in, out]  context  pointer to multigrid solver info and context
 *                           (actual type: cs_multigrid_t  *)
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_free(void  *context)
{
  CS_PROFILE_FUNC_RANGE();

  cs_multigrid_t *mg = (cs_multigrid_t *)context;

  cs_timer_t t0, t1;

  /* Initialization */

  t0 = cs_timer_time();

  /* Free subgrid info first if needed */

  for (int i = 0; i < 3; i++) {
    if (mg->lv_mg[i] != nullptr)
      cs_multigrid_free(mg->lv_mg[i]);
  }

  if (mg->setup_data != nullptr) {

    cs_multigrid_setup_data_t *mgd = mg->setup_data;

    /* Free coarse solution data */

    CS_FREE(mgd->rhs_vx);
    CS_FREE_HD(mgd->rhs_vx_buf);

    /* Destroy solver hierarchy */

    for (int i = mgd->n_levels - 1; i > -1; i--) {
      for (int j = 0; j < 2; j++) {
        cs_mg_sles_t  *mg_sles = &(mgd->sles_hierarchy[i*2+j]);
        if (mg_sles->context != nullptr && mg_sles->destroy_func)
          mg_sles->destroy_func(&(mg_sles->context));
      }
    }
    CS_FREE(mgd->sles_hierarchy);

    /* Destroy grid hierarchy */

    for (int i = mgd->n_levels - 1; i > -1; i--)
      cs_grid_destroy(mgd->grid_hierarchy + i);
    CS_FREE(mgd->grid_hierarchy);

    CS_FREE(mgd->lv_names_buffer);

    /* Destroy peconditioning-only arrays */

    CS_FREE(mgd->pc_name);
    CS_FREE(mgd->pc_aux);

    CS_FREE(mg->setup_data);
  }

  /* Update timers */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(mg->info.t_tot[0]), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a multigrid preconditioner.
 *
 * \param[in]  mg_type  type of multigrid algorithm to use
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_multigrid_pc_create(cs_multigrid_type_t  mg_type)
{
  cs_multigrid_t *mg = _multigrid_pc_create(mg_type);

  mg->info.is_pc = true;

  cs_sles_pc_t *pc = cs_sles_pc_define(mg,
                                       _multigrid_pc_get_type,
                                       _multigrid_pc_setup,
                                       _multigrid_pc_tolerance_t,
                                       _multigrid_pc_apply,
                                       cs_multigrid_free,
                                       cs_multigrid_log,
                                       cs_multigrid_copy,
                                       cs_multigrid_destroy);

  return pc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Error handler for multigrid sparse linear equation solver.
 *
 * In case of divergence or breakdown, this error handler outputs
 * postprocessing data to assist debugging, then aborts the run.
 * It does nothing in case the maximum iteration count is reached.

 * \param[in, out]  sles           pointer to solver object
 * \param[in]       state          convergence state
 * \param[in]       a              matrix
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 *
 * \return  false (do not attempt new solve)
 */
/*----------------------------------------------------------------------------*/

bool
cs_multigrid_error_post_and_abort(cs_sles_t                    *sles,
                                  cs_sles_convergence_state_t   state,
                                  const cs_matrix_t            *a,
                                  const cs_real_t               rhs[],
                                  cs_real_t                     vx[])
{
  if (state >= CS_SLES_MAX_ITERATION)
    return false;

  const cs_multigrid_t  *mg = (cs_multigrid_t *)cs_sles_get_context(sles);
  const char *name = cs_sles_get_name(sles);

  cs_multigrid_setup_data_t *mgd = mg->setup_data;
  if (mgd == nullptr)
    return false;

  int level = mgd->exit_level;

  int mesh_id = cs_post_init_error_writer_cells();
  int location_id = CS_MESH_LOCATION_CELLS;

  const cs_range_set_t  *rs = nullptr;

  if (mesh_id != 0) {

    char var_name[32];

    int lv_id = 0;
    cs_real_t *var = nullptr, *da = nullptr;

    cs_lnum_t db_size = 1;
    cs_lnum_t eb_size = 1;

    const cs_grid_t *g = mgd->grid_hierarchy[0];
    const cs_lnum_t n_base_rows = cs_grid_get_n_rows(g);
    const cs_matrix_t  *_matrix = nullptr;

    CS_MALLOC(var, cs_grid_get_n_cols_ext(g), cs_real_t);
    CS_MALLOC(da, cs_grid_get_n_cols_ext(g), cs_real_t);

    /* Output info on main level */

    cs_sles_post_error_output_def(name,
                                  mesh_id,
                                  a,
                                  rhs,
                                  vx);

    /* Output diagonal and diagonal dominance for all coarse levels */

    for (lv_id = 1; lv_id < (int)(mgd->n_levels); lv_id++) {

      g = mgd->grid_hierarchy[lv_id];

      cs_grid_get_info(g,
                       nullptr,
                       nullptr,
                       &db_size,
                       &eb_size,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr);

      _matrix = cs_grid_get_matrix(g);

      cs_matrix_copy_diagonal(_matrix, da);
      cs_grid_project_var(g, n_base_rows, da, var);
      cs_range_set_scatter(rs, CS_REAL_TYPE, db_size, var, var);
      sprintf(var_name, "Diag_%04d", lv_id);
      cs_sles_post_output_var(var_name,
                              mesh_id,
                              location_id,
                              CS_POST_WRITER_ERRORS,
                              db_size,
                              var);

      cs_grid_project_diag_dom(g, n_base_rows, var);
      cs_range_set_scatter(rs, CS_REAL_TYPE, db_size, var, var);
      sprintf(var_name, "Diag_Dom_%04d", lv_id);
      cs_sles_post_output_var(var_name,
                              mesh_id,
                              location_id,
                              CS_POST_WRITER_ERRORS,
                              db_size,
                              var);
    }

    /* Output info on current level if > 0 */

    if (level > 0) {

      cs_lnum_t ii;
      cs_lnum_t n_cells = 0;
      cs_lnum_t n_cols_ext = 0;

      cs_real_t *c_res = nullptr;

      g = mgd->grid_hierarchy[level];

      cs_grid_get_info(g,
                       nullptr,
                       nullptr,
                       &db_size,
                       &eb_size,
                       nullptr,
                       &n_cells,
                       &n_cols_ext,
                       nullptr,
                       nullptr);

      cs_grid_project_var(g, n_base_rows, mgd->rhs_vx[level*2], var);
      cs_range_set_scatter(rs, CS_REAL_TYPE, db_size, var, var);
      sprintf(var_name, "RHS_%04d", level);
      cs_sles_post_output_var(var_name,
                              mesh_id,
                              location_id,
                              CS_POST_WRITER_ERRORS,
                              db_size,
                              var);

      cs_grid_project_var(g, n_base_rows, mgd->rhs_vx[level*2+1], var);
      cs_range_set_scatter(rs, CS_REAL_TYPE, db_size, var, var);
      sprintf(var_name, "X_%04d", level);
      cs_sles_post_output_var(var_name,
                              mesh_id,
                              location_id,
                              CS_POST_WRITER_ERRORS,
                              db_size,
                              var);

      /* Compute residual */

      CS_MALLOC(c_res, n_cols_ext*db_size, cs_real_t);

      _matrix = cs_grid_get_matrix(g);

      cs_matrix_vector_multiply(_matrix,
                                mgd->rhs_vx[level*2+1],
                                c_res);

      const cs_real_t *c_rhs_lv = mgd->rhs_vx[level*2];
      for (ii = 0; ii < n_cells; ii++) {
        for (cs_lnum_t i = 0; i < db_size; i++)
          c_res[ii*db_size + i]
            = fabs(c_res[ii*db_size + i] - c_rhs_lv[ii*db_size + i]);
      }

      cs_grid_project_var(g, n_base_rows, c_res, var);
      cs_range_set_scatter(rs, CS_REAL_TYPE, db_size, var, var);

      CS_FREE(c_res);

      sprintf(var_name, "Residual_%04d", level);
      cs_sles_post_output_var(var_name,
                              mesh_id,
                              location_id,
                              CS_POST_WRITER_ERRORS,
                              db_size,
                              var);
    }

    cs_post_finalize();

    CS_FREE(da);
    CS_FREE(var);
  }

  /* Now abort */

  const char *error_type[] = {N_("divergence"), N_("breakdown")};
  int err_id = (state == CS_SLES_BREAKDOWN) ? 1 : 0;

  if (level == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("algebraic multigrid [%s]: %s after %d cycles:\n"
                "  initial residual: %11.4e; current residual: %11.4e"),
              name, _(error_type[err_id]), mgd->exit_cycle_id,
              mgd->exit_initial_residual, mgd->exit_residual);
  else
    bft_error(__FILE__, __LINE__, 0,
              _("algebraic multigrid [%s]: %s after %d cycles\n"
                "  during resolution at level %d:\n"
                "  initial residual: %11.4e; current residual: %11.4e"),
              name, _(error_type[err_id]),
              mgd->exit_cycle_id, level,
              mgd->exit_initial_residual, mgd->exit_residual);

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set plotting options for multigrid.
 *
 * \param[in, out]  mg             pointer to multigrid info and context
 * \param[in]       base_name      base plot name to activate, nullptr otherwise
 * \param[in]       use_iteration  if true, use iteration as time stamp
 *                                 otherwise, use wall clock time
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_plot_options(cs_multigrid_t  *mg,
                              const char      *base_name,
                              bool             use_iteration)
{
  if (mg != nullptr) {
    if (cs_glob_rank_id < 1 && base_name != nullptr) {

      /* Destroy previous plot if options reset */
      if (mg->cycle_plot != nullptr)
        cs_time_plot_finalize(&(mg->cycle_plot));

      /* Create new plot */
      cs_file_mkdir_default("monitoring");
      const char *probe_names[] = {base_name};
      mg->cycle_plot = cs_time_plot_init_probe(base_name,
                                               "monitoring/residual_",
                                               CS_TIME_PLOT_CSV,
                                               use_iteration,
                                               -1.,      /* force flush */
                                               0,        /* no buffer */
                                               1,        /* n_probes */
                                               nullptr,  /* probe_list */
                                               nullptr,  /* probe_coords */
                                               probe_names);

      if (use_iteration)
        mg->plot_time_stamp = 0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query the global multigrid parameters for parallel grid merging.
 *
 * \param[in]   mg                   pointer to multigrid info and context
 * \param[out]  rank_stride          number of ranks over which merging
 *                                   takes place, or nullptr
 * \param[out]  rows_mean_threshold  mean number of rows under which merging
 *                                   should be applied, or nullptr
 * \param[out]  rows_glob_threshold  global number of rows under which
 *                                   merging should be applied, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_get_merge_options(const cs_multigrid_t  *mg,
                               int                   *rank_stride,
                               int                   *rows_mean_threshold,
                               cs_gnum_t             *rows_glob_threshold)
{
#if defined(HAVE_MPI)
  if (rank_stride != nullptr)
    *rank_stride = mg->merge_stride;
  if (rows_mean_threshold != nullptr)
    *rows_mean_threshold = mg->merge_mean_threshold;
  if (rows_glob_threshold != nullptr)
    *rows_glob_threshold = mg->merge_glob_threshold;
#else
  if (rank_stride != nullptr)
    *rank_stride = 0;
  if (rows_mean_threshold != nullptr)
    *rows_mean_threshold = 0;
  if (rows_glob_threshold != nullptr)
    *rows_glob_threshold = 0;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set global multigrid parameters for parallel grid merging behavior.
 *
 * \param[in, out]  mg                   pointer to multigrid info and context
 * \param[in]       rank_stride          number of ranks over which merging
 *                                       takes place
 * \param[in]       rows_mean_threshold  mean number of rows under which
 *                                       merging should be applied
 * \param[in]       rows_glob_threshold  global number of rows under which
 *                                       merging should be applied
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_merge_options(cs_multigrid_t  *mg,
                               int              rank_stride,
                               int              rows_mean_threshold,
                               cs_gnum_t        rows_glob_threshold)
{
#if defined(HAVE_MPI)
  mg->merge_stride = rank_stride;
  if (mg->merge_stride > cs_glob_n_ranks)
    mg->merge_stride = cs_glob_n_ranks;
  mg->merge_mean_threshold = rows_mean_threshold;
  mg->merge_glob_threshold = rows_glob_threshold;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query the global multigrid parameters for parallel grid merging.
 *
 * \param[in]   mg              pointer to multigrid info and context
 * \param[out]  n_max_ranks     maximum number of MPI ranks for bottom grid
 * \param[out]  max_row_factor  maximum acceptable mean ratio of merged rows
 *                              (per MPI rank) to finest rows.
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_get_merge_bottom_options(const cs_multigrid_t  *mg,
                                      int                   *n_max_ranks,
                                      float                 *max_row_factor)
{
#if defined(HAVE_MPI)
  if (n_max_ranks != nullptr)
    *n_max_ranks = mg->merge_bottom_n_max_ranks;
  if (max_row_factor != nullptr)
    *max_row_factor = mg->merge_bottom_max_row_factor;
#else
  if (n_max_ranks != nullptr)
    *n_max_ranks = 1;
  if (max_row_factor != nullptr)
    *max_row_factor = 1.;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set global multigrid parameters for parallel grid merging behavior.
 *
 * \param[in, out]  mg              pointer to multigrid info and context
 * \param[in]       n_max_ranks     maximum number of MPI ranks for bottom grid
 * \param[in]       max_row_factor  maximum acceptable mean ratio of merged rows
 *                                  (per MPI rank) to finest rows.
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_merge_bottom_options(cs_multigrid_t  *mg,
                                      int              n_max_ranks,
                                      float            max_row_factor)
{
#if defined(HAVE_MPI)
  mg->merge_bottom_n_max_ranks = cs::min(n_max_ranks, cs_glob_n_ranks);
  mg->merge_bottom_max_row_factor = max_row_factor;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a grid associated with a given multigrid
 *        setup and level.
 *
 * If the multigrid hierarchy is not set up, or a level coarser than the
 * coarsest level is requested, nullptr is returned.

 * \param[in]  mg     pointer to multigrid info and context
 * \param[in]  level  level of the requested grid (or -1 for coarsest)
 *
 * \return  pointer to grid of requested level (nullptr id not present)
 */
/*----------------------------------------------------------------------------*/

const cs_grid_t *
cs_multigrid_get_grid(const cs_multigrid_t  *mg,
                      int                    level)
{
  const cs_grid_t *g = nullptr;

  cs_multigrid_setup_data_t *mgd = mg->setup_data;

  if (mgd != nullptr) {

    if (level < 0)
      level = mgd->n_levels - 1;

    if ((unsigned)level < mgd->n_levels)
      g = mgd->grid_hierarchy[level];

  }

  return g;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
