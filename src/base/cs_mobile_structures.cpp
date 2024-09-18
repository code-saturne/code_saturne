/*============================================================================
 * Mobile structures management.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_ale.h"
#include "cs_ast_coupling.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_file.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_parameters_check.h"
#include "cs_prototypes.h"
#include "cs_time_plot.h"
#include "cs_time_step.h"
#include "cs_timer_stats.h"
#include "cs_turbulence_model.h"
#include "cs_velocity_pressure.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mobile_structures.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_mobile_structures.c
        Mobile structures management.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef cs_real_t  _cs_real_11_t[11];  /* Vector of 11 real values */

/*! Mobile_structures type */
/*-------------------------*/

typedef struct {

  /* Base structure definitions and input */

  int            n_structs;  /*!< number of structures */

  cs_real_t      aexxst;     /*!< coefficient for the predicted displacement */
  cs_real_t      bexxst;     /*!< coefficient for the predicted displacement */

  cs_real_t      cfopre;     /*!< coefficient for the predicted force */

  cs_real_t      alpnmk;     /*!< alpha coefficient for the Newmark hht method */
  cs_real_t      betnmk;     /*!< beta coefficient for the Newmark hht method */
  cs_real_t      gamnmk;     /*!< gamma coefficient for the Newmark hht method */

  cs_real_33_t  *xmstru;     /*!< mass matrices (kg)*/
  cs_real_33_t  *xcstru;     /*!< damping matrix coefficients (kg/s) */
  cs_real_33_t  *xkstru;     /*!< spring matrix constants (kg/s2 = N/m) */

  /* Output (plotting) control */

  int                plot;               /*!< monitoring format mask
                                           0: no plot
                                           1: plot to text (.dat) format
                                           2: plot to .csv format *
                                           3: plot to both formats */

  cs_time_control_t  plot_time_control;  /*!< time control for plotting */
  char              *plot_dir_name;      /*!< monitoring output directory */

  /* Computed structure values */

  cs_real_3_t   *xstr;       /*!< displacement vectors compared to structure
                              *   positions in the initial mesh (m) */
  cs_real_3_t  *xsta;        /*!< values of xstr at the previous time step */
  cs_real_3_t  *xstp;        /*!< predicted values of xstr */
  cs_real_3_t  *xstreq;      /*!< equilibrum positions of a structure (m) */

  cs_real_3_t  *xpstr;       /*!< velocity vectors (m/s) */
  cs_real_3_t  *xpsta;       /*!< xpstr at previous time step */

  cs_real_3_t  *xppstr;      /*!< acceleration vectors (m/s2) */
  cs_real_3_t  *xppsta;      /*!< acceleration vectors at previous
                              *   time step (m/s2) */

  cs_real_3_t  *forstr;      /*!< force vectors acting on the structure (N) */
  cs_real_3_t  *forsta;      /*!< forstr at previous time step (N) */
  cs_real_3_t  *forstp;      /*!< predicted force vectors (N) */

  cs_real_t  *dtstr;         /*!< time step used to solve structure movements */

  /* Association with mesh */

  int        *idfstr;        /*!< structure number associated to each
                              *   boundary face:
                              *   - 0 if face is not coupled to a structure
                              *   - if > 0, internal structure id + 1
                              *   - if < 0, - code_aster instance id  - 1 */

  /* Plotting */

  int            n_plots;    /*!< number of plots for format */

  cs_time_plot_t  **plot_files[2];  /*!< Associated plot files */

} cs_mobile_structures_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_mobile_structures_t  *_mobile_structures = nullptr;

/* Arrays allowing return to initial state at end of ALE iteration.
 * - Mass flux: save at the first call of cs_theta_scheme_update_var.
 * - Gradient BC's for P and U (since we use a reconstruction to compute
 *   the real BC's); this part might not really be necessary.
 * - The initial pressure (since the initial pressure is also overwritten when
 *   nterup > 1); maybe this could be optimized... */

static _cs_real_11_t *_bc_coeffs_save = nullptr;
static cs_real_t *_pr_save = nullptr;

static int _post_out_stat_id = -1;

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Number of code_aster couplings */
int cs_glob_ast_coupling_n_couplings = 0;

/*! Maximum number of implicitation iterations of the structure displacement */
int cs_glob_mobile_structures_i_max = 1;

/*! Relative precision of implicitation of the structure displacement */
double cs_glob_mobile_structures_i_eps = 1e-5;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new mobile structures handling structure.
 *
 * \return  pointer to newly created structure.
 */
/*----------------------------------------------------------------------------*/

static cs_mobile_structures_t *
_mobile_structures_create(void)
{
  cs_mobile_structures_t *ms;

  BFT_MALLOC(ms, 1, cs_mobile_structures_t);

  ms->n_structs = 0;

  ms->aexxst = -cs_math_big_r;
  ms->bexxst = -cs_math_big_r;
  ms->cfopre = -cs_math_big_r;

  ms->alpnmk = 0.;
  ms->betnmk = -cs_math_big_r;
  ms->gamnmk = -cs_math_big_r;

  ms->plot = 2;

  cs_time_control_init_by_time_step(&(ms->plot_time_control),
                                    -1,
                                    -1,
                                    1,       /* nt_interval */
                                    true,    /* at_start */
                                    false);  /* at end */

  ms->xmstru = nullptr;
  ms->xcstru = nullptr;
  ms->xkstru = nullptr;

  ms->xstr = nullptr;
  ms->xsta = nullptr;
  ms->xstp = nullptr;
  ms->xstreq = nullptr;

  ms->xpstr = nullptr;
  ms->xpsta = nullptr;

  ms->xppstr = nullptr;
  ms->xppsta = nullptr;

  ms->forstr = nullptr;
  ms->forsta = nullptr;
  ms->forstp = nullptr;

  ms->dtstr = nullptr;

  /* Plot info */

  ms->n_plots = 0;

  ms->plot_dir_name = nullptr;
  const char dir_name[] = "monitoring";
  BFT_MALLOC(ms->plot_dir_name, strlen(dir_name) + 1, char);
  strcpy(ms->plot_dir_name, dir_name);

  ms->plot_files[0] = nullptr;
  ms->plot_files[1] = nullptr;

  /* Face association */

  ms->idfstr = nullptr;

  return ms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a mobile structures handling structure.
 *
 * \param[in, out]  pointer to mobile structures object to destroy.
 */
/*----------------------------------------------------------------------------*/

static void
_mobile_structures_destroy(cs_mobile_structures_t  **ms)
{
  cs_mobile_structures_t  *_ms = *ms;

  BFT_FREE(_ms->xmstru);
  BFT_FREE(_ms->xcstru);
  BFT_FREE(_ms->xkstru);

  BFT_FREE(_ms->xstr);
  BFT_FREE(_ms->xsta);
  BFT_FREE(_ms->xstp);
  BFT_FREE(_ms->xstreq);

  BFT_FREE(_ms->xpstr);
  BFT_FREE(_ms->xpsta);

  BFT_FREE(_ms->xppstr);
  BFT_FREE(_ms->xppsta);

  BFT_FREE(_ms->forstr);
  BFT_FREE(_ms->forsta);
  BFT_FREE(_ms->forstp);

  BFT_FREE(_ms->dtstr);

  /* Plot info */

  for (int fmt = CS_TIME_PLOT_DAT; fmt <= CS_TIME_PLOT_CSV; fmt++) {
    if (_ms->plot_files[fmt] != nullptr) {

      for (int i = 0; i < _ms->n_plots; i++) {
        cs_time_plot_t *p = _ms->plot_files[fmt][i];

        if (p != nullptr) {
          cs_time_plot_finalize(&p);
          _ms->plot_files[fmt][i] = nullptr;
        }
      }

      BFT_FREE(_ms->plot_files[fmt]);

    }
  }

  BFT_FREE(_ms->plot_dir_name);

  /* Face association */

  BFT_FREE(_ms->idfstr);

  BFT_FREE(_ms);
  *ms = _ms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize mobile structures with ALE for internal coupling.
 *
 * In case of multiple calls of this function, only added structures
 * are modified (so that structures may be added by calling this function
 * with increasing values of n_structures).
 *
 * \param[in, out]  ms            pointer to mobile structure
 * \param[in]       n_structures  number of internal structures
 */
/*----------------------------------------------------------------------------*/

static void
_init_internal_structures(cs_mobile_structures_t *ms,
                          int                     n_structures)
{
  int n_structs_prev = ms->n_structs;
  ms->n_structs = n_structures;

  BFT_REALLOC(ms->xmstru, n_structures, cs_real_33_t);
  BFT_REALLOC(ms->xcstru, n_structures, cs_real_33_t);
  BFT_REALLOC(ms->xkstru, n_structures, cs_real_33_t);

  BFT_REALLOC(ms->xstr, n_structures, cs_real_3_t);
  BFT_REALLOC(ms->xsta, n_structures, cs_real_3_t);
  BFT_REALLOC(ms->xstp, n_structures, cs_real_3_t);
  BFT_REALLOC(ms->xstreq, n_structures, cs_real_3_t);

  BFT_REALLOC(ms->xpstr, n_structures, cs_real_3_t);
  BFT_REALLOC(ms->xpsta, n_structures, cs_real_3_t);

  BFT_REALLOC(ms->xppstr, n_structures, cs_real_3_t);
  BFT_REALLOC(ms->xppsta, n_structures, cs_real_3_t);

  BFT_REALLOC(ms->forstr, n_structures, cs_real_3_t);
  BFT_REALLOC(ms->forsta, n_structures, cs_real_3_t);
  BFT_REALLOC(ms->forstp, n_structures, cs_real_3_t);

  BFT_REALLOC(ms->dtstr, n_structures, cs_real_t);

  for (int i = n_structs_prev; i < n_structures; i++) {

    ms->dtstr[i] = 0;
    for (int j = 0; j < 3; j++) {
      ms->xstr[i][j]   = 0;
      ms->xpstr[i][j]  = 0;
      ms->xppstr[i][j] = 0;
      ms->xsta[i][j]   = 0;
      ms->xpsta[i][j]  = 0;
      ms->xppsta[i][j] = 0;
      ms->xstp[i][j]   = 0;
      ms->forstr[i][j] = 0;
      ms->forsta[i][j] = 0;
      ms->forstp[i][j] = 0;
      ms->xstreq[i][j] = 0;
      for (int k = 0; k < 3; k++) {
        ms->xmstru[i][j][k] = 0;
        ms->xcstru[i][j][k] = 0;
        ms->xkstru[i][j][k] = 0;
      }
    }

  }

  _post_out_stat_id = cs_timer_stats_id_by_name("postprocessing_output");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Newmark HHT method to solve second order linear differential
 *         equation.
 *
 * The equation id of the form:
 *    M.X'' + C.X' + K.(X+X0) = F
 *
 * Where X is a vector field, and M, C, and K are arbitrary 3x3 matrices.
 *
 * The solution is of displacement type.
 *
 * \param[in]   structure_id  id of matching structure
 * \param[in]   alpnmk        alpha coefficient
 * \param[in]   betnmk        beta coefficient
 * \param[in]   gamnmk        gamma coefficient
 * \param[in]   xm            system's mass matrix
 * \param[in]   xc            system's friction matrix
 * \param[in]   xk            system's stiffness matrix
 * \param[in]   xn0           initial displacement
 * \param[out]  xn            displacement at time n
 * \param[out]  xpn           velocity at time n
 * \param[out]  xppn          acceleration at time n
 * \param[in]   xnm1          displacement at time n-1
 * \param[in]   xpnm1         velocity at time n-1
 * \param[in]   xppnm1        acceleration at time n-1
 * \param[in]   xfn           force at time n
 * \param[in]   xfnm1         force at time n-1
 * \param[in]   dt            time step
 */
/*----------------------------------------------------------------------------*/

static void
_newmark(int        structure_id,
         cs_real_t  alpnmk,
         cs_real_t  betnmk,
         cs_real_t  gamnmk,
         cs_real_t  xm[3][3],
         cs_real_t  xc[3][3],
         cs_real_t  xk[3][3],
         cs_real_t  xn0[3],
         cs_real_t  xn[3],
         cs_real_t  xpn[3],
         cs_real_t  xppn[3],
         cs_real_t  xnm1[3],
         cs_real_t  xpnm1[3],
         cs_real_t  xppnm1[3],
         cs_real_t  xfn[3],
         cs_real_t  xfnm1[3],
         cs_real_t  dt)
{
  /* Null displacement criterion */
  const cs_real_t epsdet = 1e-12;

  /* Equation coefficients */
  cs_real_t a0 = 1.0/betnmk/cs_math_pow2(dt);
  cs_real_t a1 = (1.0+alpnmk)*gamnmk/betnmk/dt;
  cs_real_t a2 = 1.0/betnmk/dt;
  cs_real_t a3 = 1.0/2.0/betnmk - 1.0;
  cs_real_t a4 = (1.0+alpnmk)*gamnmk/betnmk - 1.0;
  cs_real_t a5 = (1.0+alpnmk)*dt*(gamnmk/2.0/betnmk - 1.0);
  cs_real_t a6 = dt*(1.0 - gamnmk);
  cs_real_t a7 = gamnmk*dt;

  double a[3][3], b1[3], b2[3], b[3];
  double det, det1, det2, det3;

  for (int ii = 0; ii < 3; ii++) {
    for (int jj = 0; jj < 3; jj++) {
      a[jj][ii] =  (1.0+alpnmk)*xk[jj][ii]
                  + a1*xc[jj][ii] + a0*xm[jj][ii];
    }
    b[ii]  = (1.0+alpnmk)*xfn[ii] - alpnmk*xfnm1[ii];
    b1[ii] = a0*xnm1[ii] + a2*xpnm1[ii] + a3*xppnm1[ii];
    b2[ii] = a1*xnm1[ii] + a4*xpnm1[ii] + a5*xppnm1[ii];
  }

  for (int ii = 0; ii < 3; ii++) {
    for (int jj = 0; jj < 3; jj++) {
      b[ii] +=   xm[jj][ii]*b1[jj] + xc[jj][ii]*b2[jj]
               + xk[jj][ii]*(alpnmk*xnm1[jj] + xn0[jj]);
    }
  }

  det =   a[0][0]*a[1][1]*a[2][2]
        + a[0][1]*a[1][2]*a[2][0]
        + a[0][2]*a[1][0]*a[2][1]
        - a[0][2]*a[1][1]*a[2][0]
        - a[0][1]*a[1][0]*a[2][2]
        - a[0][0]*a[1][2]*a[2][1];

  if (fabs(det) < epsdet) {
    cs_log_printf
      (CS_LOG_DEFAULT,
       _("@\n"
         "@ @@ Warning: ALE displacement of internal structures\n"
         "@    =======\n"
         "@  Structure: %d\n"
         "@  The absolute value of the discriminant of the\n"
         "@    displacement matrix is: %14.5e\n"
         "@  The matrix is considered not inversible\n"
         "@    (limit value fixed to %14.5e)\n"
         "@\n"
         "@  Calculation abort\n"),
       structure_id, fabs(det), epsdet);

    cs_time_step_define_nt_max(cs_glob_time_step->nt_cur);
  }

  det1 =   b[0]*a[1][1]*a[2][2]
         + b[1]*a[1][2]*a[2][0]
         + b[2]*a[1][0]*a[2][1]
         - b[2]*a[1][1]*a[2][0]
         - b[1]*a[1][0]*a[2][2]
         - b[0]*a[1][2]*a[2][1];

  det2 =   a[0][0]*b[1]*a[2][2]
         + a[0][1]*b[2]*a[2][0]
         + a[0][2]*b[0]*a[2][1]
         - a[0][2]*b[1]*a[2][0]
         - a[0][1]*b[0]*a[2][2]
         - a[0][0]*b[2]*a[2][1];

  det3 =   a[0][0]*a[1][1]*b[2]
         + a[0][1]*a[1][2]*b[0]
         + a[0][2]*a[1][0]*b[1]
         - a[0][2]*a[1][1]*b[0]
         - a[0][1]*a[1][0]*b[2]
         - a[0][0]*a[1][2]*b[1];

  xn[0] = det1/det;
  xn[1] = det2/det;
  xn[2] = det3/det;

  for (int ii = 0; ii < 3; ii++) {
    xppn[ii] = a0*(xn[ii]-xnm1[ii]) - a2*xpnm1[ii]  - a3*xppnm1[ii];
    xpn[ii]  = xpnm1[ii]            + a6*xppnm1[ii] + a7*xppn[ii];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize time plot for internal mobile structures.
 *
 * \param[in, out]  ms   pointer to mobile structures, or nullptr
 */
/*----------------------------------------------------------------------------*/

static void
_init_time_plot(cs_mobile_structures_t   *ms)
{
  if (ms == nullptr || cs_glob_rank_id > 0)
    return;

  bool use_iteration = (cs_glob_time_step->is_local) ? true : false;

  if (cs_file_mkdir_default(ms->plot_dir_name) != 0)
    bft_error(__FILE__, __LINE__, 0,
              _("The %s directory cannot be created"), ms->plot_dir_name);

  const int n_plots = 12;

  const char *name[12] = {"displacement x",
                          "displacement y",
                          "displacement z",
                          "velocity x",
                          "velocity y",
                          "velocity z",
                          "acceleration x",
                          "acceleration y",
                          "acceleration z",
                          "force x",
                          "force y",
                          "force z"};

  char *file_prefix;
  const char base_prefix[] = "structures_";
  size_t l = strlen(base_prefix);

  if (ms->plot_dir_name != nullptr) {
    l += strlen(ms->plot_dir_name) + 1;
    BFT_MALLOC(file_prefix, l+1, char);
    snprintf(file_prefix, l+1, "%s/%s", ms->plot_dir_name, base_prefix);
  }
  else {
    BFT_MALLOC(file_prefix, l+1, char);
    strncpy(file_prefix, base_prefix, l+1);
  }
  file_prefix[l] = '\0';

  float flush_wtime = -1;
  int   n_buffer_steps = -1;

  cs_time_plot_get_flush_default(&flush_wtime, &n_buffer_steps);

  for (int i = 0; i < n_plots; i++) {
    for (int fmt = CS_TIME_PLOT_DAT; fmt <= CS_TIME_PLOT_CSV; fmt++) {
      int fmt_mask = fmt + 1;

      if (ms->plot & fmt_mask) {
        if (i == 0) {
          BFT_REALLOC(ms->plot_files[fmt], n_plots, cs_time_plot_t *);
          for (int j = 0; j < n_plots; j++)
            ms->plot_files[fmt][j] = nullptr;
        }

        ms->plot_files[fmt][i] =
          cs_time_plot_init_struct(name[i],
                                   file_prefix,
                                   static_cast<cs_time_plot_format_t>(fmt),
                                   use_iteration,
                                   flush_wtime,
                                   n_buffer_steps,
                                   ms->n_structs,
                                   (cs_real_t *)ms->xmstru,
                                   (cs_real_t *)ms->xcstru,
                                   (cs_real_t *)ms->xkstru);
      }
    }
  }

  ms->n_plots = n_plots;

  BFT_FREE(file_prefix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output single time plot series for internal mobile structures.
 *
 * \param[in]  pl_it   plot id
 * \param[in]  ms      pointer to mobile structures, or nullptr
 * \param[in]  n       number of values
 * \param[in]  nt_cur  current time step value
 * \param[in]  t_cur   current time value
 * \param[in]  val     values to output (size: n)
 */
/*----------------------------------------------------------------------------*/

static void
_time_plot_write(int                      pl_id,
                 cs_mobile_structures_t  *ms,
                 int                      n,
                 int                      nt_cur,
                 cs_real_t                t_cur,
                 const cs_real_t           val[])
{
  for (int fmt = CS_TIME_PLOT_DAT; fmt <= CS_TIME_PLOT_CSV; fmt++) {
    assert(pl_id > -1 && pl_id < ms->n_plots);

    if (ms->plot_files[fmt] != nullptr) {
      cs_time_plot_t *p = ms->plot_files[fmt][pl_id];
      if (p != nullptr)
        cs_time_plot_vals_write(p, nt_cur, t_cur, n, val);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output time plot for internal mobile structures.
 *
 * \param[in, out]  ms   pointer to mobile structures, or nullptr
 * \param[in]       ts   pointer to time step structure
 */
/*----------------------------------------------------------------------------*/

static void
_output_time_plots(cs_mobile_structures_t  *ms,
                   const cs_time_step_t    *ts)
{
  if (cs_glob_rank_id > 0)
    return;

  /* Main processing */

  cs_real_t *vartmp;
  BFT_MALLOC(vartmp, ms->n_structs, cs_real_t);

  int pl_id = 0;

  cs_real_3_t *v_pointers[4] = {ms->xstr, ms->xpstr, ms->xpstr, ms->forstr};

  for (int i = 0; i < 4; i++) {

    cs_real_3_t *v = v_pointers[i];
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < ms->n_structs; k++)
        vartmp[k] = v[k][j];

      _time_plot_write(pl_id, ms, ms->n_structs, ts->nt_cur, ts->t_cur, vartmp);
      pl_id++;
    }

  }

  BFT_FREE(vartmp);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize mobile structures with ALE for internal coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_setup(void)
{
  cs_mobile_structures_t *ms = _mobile_structures;

  if (ms == nullptr)
    return;

  /* Internal structures */

  int monitor = 1;

  cs_gui_mobile_mesh_init_structures(ms->n_structs,
                                     &(ms->aexxst),
                                     &(ms->bexxst),
                                     &(ms->cfopre),
                                     &monitor,
                                     (cs_real_t *)ms->xstp,
                                     (cs_real_t *)ms->xstreq,
                                     (cs_real_t *)ms->xpstr);

  const cs_time_step_t *ts = cs_glob_time_step;
  int is_restart = (ts->nt_prev > 0) ? 1 : 0;

  cs_user_fsi_structure_define(is_restart,
                               ms->n_structs,
                               &(ms->plot),
                               &(ms->plot_time_control),
                               &(ms->aexxst),
                               &(ms->bexxst),
                               &(ms->cfopre),
                               ms->xstp,
                               ms->xpstr,
                               ms->xstreq);

  if (ms->aexxst < -0.5*cs_math_big_r)
    ms->aexxst = 0.5;
  if (ms->bexxst < -0.5*cs_math_big_r)
    ms->bexxst = 0.;
  if (cs_glob_mobile_structures_i_max == 1) {
    if (ms->cfopre < -0.5*cs_math_big_r)
      ms->cfopre = 2.0;
  }
  else
    ms->cfopre = 1.0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize mobile structures with ALE for internal coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_initialize(void)
{
  cs_mobile_structures_t *ms = _mobile_structures;
  int n_structs = (ms != nullptr) ? ms->n_structs : 0;
  int n_ast_structs = cs_glob_ast_coupling_n_couplings;

  if (n_structs + n_ast_structs == 0)
    return;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  /* Reserve idfstr array */

  BFT_REALLOC(ms->idfstr, n_b_faces, int);

  int *idfstr = ms->idfstr;
  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    idfstr[i] = 0;

  /* Definition of idfstr
     -------------------- */

  /* Associate internal (code_saturne) and external (code_aster) structures */

  cs_gui_mobile_mesh_bc_structures(idfstr);
  cs_user_fsi_structure_num(cs_glob_domain, idfstr);

  /* Count number of internal and external structures

     Remark: currently, all external structures are associated to the
     same code_aster domain, so this information not used for those
     structures. */

  int m_vals[2] = {0, 0};
  cs_lnum_t n_ast_faces = 0;

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    int str_num = idfstr[i];
    if (str_num > m_vals[0])
      m_vals[0] = str_num;
    else if (-str_num > m_vals[1]) {
      m_vals[1] = -str_num;
      n_ast_faces += 1;
    }
  }

  cs_parall_max(2, CS_INT_TYPE, m_vals);

  /* Compare number os structures to restart */

  if (m_vals[0] > ms->n_structs) {
    cs_parameters_error
      (CS_ABORT_IMMEDIATE,
       _("Internal mobile structures"),
       _("The number of referenced structures is greater than the\n"
         "number of defined structures:\n"
         "  Number of defined structures: %d\n"
         "  Number of referenced structures: %d\n"
         "\n"
         "Check the coupled boundary structure associations."),
       ms->n_structs, m_vals[0]);
  }

  if (ms->n_structs > 0) {
    for (int i = 0; i < ms->n_structs; i++)
      ms->dtstr[i] = cs_glob_time_step->dt[0];
  }

  /* Prepare and exchange mesh info with code_aster
     ---------------------------------------------- */

  if (n_ast_structs > 0) {

    cs_lnum_t *face_ids = nullptr;
    BFT_MALLOC(face_ids, n_ast_faces, cs_lnum_t);

    cs_lnum_t count = 0;

    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      if (idfstr[i] < 0) {
        face_ids[count] = i;
        count += 1;
      }
    }

    const cs_real_t almax = cs_glob_turb_ref_values->almax;

    /* Exchange code_aster coupling parameters */
    cs_ast_coupling_initialize(cs_glob_mobile_structures_i_max,
                               cs_glob_mobile_structures_i_eps);

    /* Send geometric information to code_aster */
    cs_ast_coupling_geometry(n_ast_faces, face_ids, almax);

    BFT_FREE(face_ids);

  }

  /* In no mobile structures are used, deallocate idfstr
     and indicate that no implicitation iterations for the structure
     displacement will be needed. */

  if (n_structs + n_ast_structs == 0) {
    cs_glob_mobile_structures_i_max = 1;

    BFT_FREE(ms->idfstr);
    idfstr = nullptr;
  }

  if (n_structs && ms->plot > 0)
    _init_time_plot(ms);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize mobile structures with ALE for internal coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_finalize(void)
{
  if (_mobile_structures != nullptr)
    _mobile_structures_destroy(&_mobile_structures);

  BFT_FREE(_bc_coeffs_save);
  BFT_FREE(_pr_save);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log structures and coupling information
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_log_setup(void)
{
  cs_mobile_structures_t *ms = _mobile_structures;
  int n_structs = (ms != nullptr) ? ms->n_structs : 0;
  int n_ast_structs = cs_glob_ast_coupling_n_couplings;

  cs_log_t log = CS_LOG_SETUP;

  if (n_structs + n_ast_structs == 0) {
    cs_log_printf(log,
                  _("\n"
                    "ALE: no coupled structures\n\n"));
    return;
  }

  cs_log_printf(log,
                _("\n"
                  "ALE displacement with coupled structures\n"
                  "-----------------------------------------\n\n"));

  char fmt_type[32];

  if (n_structs > 0) {
    if (ms->plot & (CS_TIME_PLOT_DAT + 1)) {
      if (ms->plot & (CS_TIME_PLOT_CSV + 1))
        strncpy(fmt_type, ".dat, .csv", 31);
      else
        strncpy(fmt_type, ".dat", 31);
    }
    else if (ms->plot & (CS_TIME_PLOT_CSV + 1))
      strncpy(fmt_type, ".dat", 31);
    else
      strncpy(fmt_type, "none", 31);
  }
  else
    strncpy(fmt_type, "none", 31);

  fmt_type[31] = '\0';

  /* Set Newmark coefficients if not defined by user */

  if (ms->betnmk < -0.5*cs_math_big_r)
    ms->betnmk = cs_math_pow2(1. - ms->alpnmk) / 4.;
  if (ms->gamnmk < -0.5*cs_math_big_r)
    ms->gamnmk = (1. - 2.*ms->alpnmk) / 2.;

  assert(ms != nullptr);
  cs_log_printf(log,
                ("  Number of internal structures: %d\n"
                 "\n"
                 "    Newmark coefficients:\n"
                 "      alpnmk: %12.4e\n"
                 "      betnmk: %12.4e\n"
                 "      gamnmk: %12.4e\n"
                 "\n"
                 "    Monitoring output interval for structures:\n"
                 "      format: %s\n"
                 "      nthist: %d\n"
                 "      frhist: %g\n"),
                ms->n_structs,
                ms->alpnmk, ms->betnmk, ms->gamnmk,
                fmt_type,
                ms->plot_time_control.interval_nt,
                ms->plot_time_control.interval_t);

  if (cs_glob_mobile_structures_i_max == 1) {
    cs_log_printf(log,
                  ("\n"
                   "  Explicit coupling scheme\n"
                   "    Coefficients:\n"
                   "      aexxst: %12.4e\n"
                   "      bexxst: %12.4e\n"
                   "      cfopre: %12.4e\n\n"),
                  ms->aexxst, ms->bexxst, ms->cfopre);
  }
  else {
    cs_log_printf(log,
                  ("\n"
                   "  Implicit coupling scheme\n"
                   "    maximum number of inner iterations: %d\n"
                   "    convergence threshold:              %g\n\n"),
                  cs_glob_mobile_structures_i_max,
                  cs_glob_mobile_structures_i_eps);
  }

  if (n_ast_structs > 0) {
    cs_log_printf(log,
                  ("  Number of coupled code_aster structures: %d\n\n"),
                  n_ast_structs);
  }

  cs_log_separator(log);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query number of internal mobile structures defined.
 *
 * \return  number of internal mobile structures
 */
/*----------------------------------------------------------------------------*/

int
cs_mobile_structures_get_n_structures(void)
{
  int retval = 0;

  if (_mobile_structures != nullptr)
    retval = _mobile_structures->n_structs;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add internal mobile structures.
 *
 * This function may be called multiple time to change the number of
 * mobile structures.
 *
 * \param[in]   n_structures  number of internal mobile structures
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_add_n_structures(int  n_structures)
{
  if (n_structures > 0) {
    cs_mobile_structures_t *ms = _mobile_structures;

    if (ms == nullptr) {
      ms = _mobile_structures_create();
      _mobile_structures = ms;
    }

    _init_internal_structures(ms, ms->n_structs + n_structures);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Newmark coefficients for internal mobile structures.
 *
 * \param[in]   alpha  alpha coefficient for Newmark algorithm
 * \param[in]   beta   beta coefficient for Newmark algorithm
 * \param[in]   gamma  gamma coefficient for Newmark algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_set_newmark_coefficients(cs_real_t  alpha,
                                              cs_real_t  beta,
                                              cs_real_t  gamma)
{
  if (   alpha < 0 || alpha > 1
      || beta < 0  || beta > 0.5
      || gamma < 0 || gamma > 1)
    cs_parameters_error
      (CS_ABORT_IMMEDIATE,
       _("Internal mobile structures"),
       _("%s: The Newmark coefficients should be in the following ranges:\n"
         "\n"
         "  alpha: [0, 1]\n"
         "  beta:  [0, 0.5]\n"
         "  gamma: [0, 1]\n"
         "\n"
         "Here, we have:\n"
         "  alpha: %g\n"
         "  beta:  %g\n"
         "  gamma: %g\n"
         "\n"
         "Check the provided parameters."),
       __func__, alpha, beta, gamma);

  cs_mobile_structures_t *ms = _mobile_structures;
  if (ms == nullptr)
    ms = _mobile_structures_create();

  ms->alpnmk = alpha;
  ms->betnmk = beta;
  ms->gamnmk = gamma;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predict displacement of mobile structures with ALE.
 *
 * \param[in]   itrale   ALE iteration number
 * \param[in]   italim   implicit coupling iteration number
 * \param[in]   ineefl   indicate whether fluxes should be saved
 * \param[out]  impale   imposed displacement indicator
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_prediction(int  itrale,
                                int  italim,
                                int  ineefl,
                                int  impale[])
{
  cs_mobile_structures_t *ms = _mobile_structures;
  int n_structs = (ms != nullptr) ? ms->n_structs : 0;
  int n_ast_structs = cs_glob_ast_coupling_n_couplings;

  if (n_structs + n_ast_structs == 0)
    return;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  cs_real_3_t *disale = nullptr;

  /* Predict structures displacement
     ------------------------------- */

  /* Internal structures:
   *
   * When initializing ALE (itrale=0), xstp contains:
   * - The value of the initial structures displacement if the user has touched
   *   them (restart or not).
   * - 0 if computation starts with structures
   * - The displacement used for the previous computation if restarted with no
   *   modification by user.
   *
   * Its value must be transferred to xstr (which is used by  Newmark).
   * In the following iterations (itrale>0) we use the standard computation
   * schemes for xstp. */

  if (n_structs > 0) {

    cs_real_t dt_ref = cs_glob_time_step->dt[0];
    for (int i = 0; i < ms->n_structs; i++) {
      ms->dtstr[i] = dt_ref;
    }

    if (itrale == 0) {
      for (int i = 0; i < n_structs; i++) {
        for (int j= 0; j < 3; j++)
          ms->xstr[i][j] = ms->xstp[i][j];
      }
    }

    else {

      /* Explicit coupling scheme */
      if (cs_glob_mobile_structures_i_max == 1) {
        cs_real_t aexxst = ms->aexxst;
        cs_real_t bexxst = ms->bexxst;

        for (int i = 0; i < n_structs; i++) {
          cs_real_t dt = ms->dtstr[i];
          for (int j= 0; j < 3; j++) {
            ms->xstp[i][j] =   ms->xstr[i][j]
                             + aexxst*dt*ms->xpstr[i][j]
                             + bexxst*dt*(ms->xpstr[i][j]-ms->xpsta[i][j]);
          }
        }
      }

      /* Implicit coupling scheme */

      else {
        for (int i = 0; i < n_structs; i++) {
          for (int j= 0; j < 3; j++) {
            ms->xstp[i][j] = ms->xstr[i][j];
          }
        }
      }

    }

    const cs_lnum_t *b_face_vtx_idx = m->b_face_vtx_idx;
    const cs_lnum_t *b_face_vtx = m->b_face_vtx_lst;

    cs_field_t  *f_displ = cs_field_by_name("mesh_displacement");
    disale = (cs_real_3_t *)(f_displ->val);

    const int *idfstr = ms->idfstr;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      int str_num = idfstr[face_id];

      /* Internal structures */

      if (str_num > 0) {
        int str_id = str_num -1;
        cs_lnum_t s_id = b_face_vtx_idx[face_id];
        cs_lnum_t e_id = b_face_vtx_idx[face_id+1];
        for (cs_lnum_t j = s_id; j < e_id; j++) {
          cs_lnum_t vtx_id = b_face_vtx[j];
          impale[vtx_id] = 1;
          for (cs_lnum_t k = 0; k < 3; k++)
            disale[vtx_id][k] = ms->xstp[str_id][k];
        }
      }

      /* External structures (code_aster) */

      else if (str_num < 0) {

        cs_lnum_t s_id = b_face_vtx_idx[face_id];
        cs_lnum_t e_id = b_face_vtx_idx[face_id+1];
        for (cs_lnum_t j = s_id; j < e_id; j++) {
          cs_lnum_t vtx_id = b_face_vtx[j];
          impale[vtx_id] = 1;
        }

      }

    }

  }

  if (n_ast_structs > 0) {

    /* If itrale = 0, we do nothing for now, but should receive the
       initial displacements coming from code_aster. */

    if (itrale > 0) {
      /* Receive predicted displacements and fill disale */

      cs_ast_coupling_compute_displacement(disale);
    }

  }

  /* Displacement at previous time step and save flux and pressure.
     -------------------------------------------------------------- */

  if (italim == 1) {

    for (int i = 0; i < n_structs; i++) {
      for (int j= 0; j < 3; j++) {
        ms->xsta[i][j] = ms->xstr[i][j];
        ms->xpsta[i][j] = ms->xpstr[i][j];
        ms->xppsta[i][j] = ms->xppstr[i][j];
      }
    }

    if (ineefl == 1) {

      /* Save BC coefficients.

         Using separate values for velocity and pressure could
         also make this more readable and safer. */

      cs_real_3_t *coefau = (cs_real_3_t *)CS_F_(vel)->bc_coeffs->a;
      cs_real_33_t *coefbu = (cs_real_33_t *)CS_F_(vel)->bc_coeffs->b;

      cs_real_t *coefap = CS_F_(p)->bc_coeffs->a;
      cs_real_t *coefbp = CS_F_(p)->bc_coeffs->b;

      BFT_REALLOC(_bc_coeffs_save, n_b_faces, _cs_real_11_t);
      _cs_real_11_t *cofale = _bc_coeffs_save;

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

        cofale[face_id][0] = coefap[face_id];

        cofale[face_id][1] = coefau[face_id][0];
        cofale[face_id][2] = coefau[face_id][1];
        cofale[face_id][3] = coefau[face_id][2];

        cofale[face_id][4] = coefbp[face_id];

        /* coefficient B is supposed to be symmetric */

        cofale[face_id][5]  = coefbu[face_id][0][0];
        cofale[face_id][6]  = coefbu[face_id][1][1];
        cofale[face_id][7]  = coefbu[face_id][2][2];
        cofale[face_id][8]  = coefbu[face_id][1][0];
        cofale[face_id][9]  = coefbu[face_id][2][1];
        cofale[face_id][10] = coefbu[face_id][2][0];

      }

      /* Backup of pressure */

      if (cs_glob_velocity_pressure_param->nterup > 1) {
        cs_lnum_t n_vals = m->n_cells_with_ghosts;
        const cs_real_t *cvara_pr = CS_F_(p)->val_pre;

        BFT_REALLOC(_pr_save, n_b_faces, cs_real_t);
        cs_real_t *xprale = _pr_save;

        for (cs_lnum_t i = 0; i < n_vals; i++)
          xprale[i] = cvara_pr[i];
      }

    } /* ineefl */

  } /* italim */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Displacement of mobile structures with ALE for internal coupling.
 *
 * \param[in]       itrale   ALE iteration number
 * \param[in]       italim   implicit coupling iteration number
 * \param[in, out]  itrfin   indicator for last iteration of implicit coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_displacement(int   itrale,
                                  int   italim,
                                  int  *itrfin)
{
  /* Number of internal and external couplings */

  cs_mobile_structures_t *ms = _mobile_structures;

  int n_structs = (ms != nullptr) ? ms->n_structs : 0;
  int n_ast_structs = cs_glob_ast_coupling_n_couplings;

  if (n_structs + n_ast_structs == 0)
    return;

  /* Initialization */

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  const cs_field_t *f_b_forces = cs_field_by_name("boundary_forces");
  cs_real_3_t *b_force = (cs_real_3_t *)f_b_forces->val;

  /* Divide by face surface if boundary forces field is extensive */
  bool  div_by_surf = (f_b_forces->type & CS_FIELD_EXTENSIVE) ? true : false;
  cs_real_t  *b_face_surf = cs_glob_mesh_quantities->b_face_surf;

  cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(mesh_u));

  const cs_time_step_t *ts = cs_glob_time_step;

  /* Compute forces on the structures
     -------------------------------- */

  for (int i = 0; i < n_structs; i++) {
    for (int j = 0; j < 3; j++) {
      ms->forsta[i][j] = ms->forstr[i][j];
      ms->forstr[i][j] = 0.;
    }
  }

  cs_real_3_t *forast = nullptr;
  if (n_ast_structs > 0)
    forast = cs_ast_coupling_get_fluid_forces_pointer();

  /* Allocate a temporary array */

  int *idfstr = ms->idfstr;

  cs_lnum_t indast = 0;
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    int str_num = idfstr[face_id];
    if (str_num > 0) {
      int i = str_num - 1;
      for (cs_lnum_t j = 0; j < 3; j++)
        ms->forstr[i][j] += b_force[face_id][j];
    }
    else if (str_num < 0) {
      if (div_by_surf) {
        for (cs_lnum_t j = 0; j < 3; j++)
          forast[indast][j] = b_force[face_id][j] / b_face_surf[face_id];
      }
      else {
        for (cs_lnum_t j = 0; j < 3; j++)
          forast[indast][j] = b_force[face_id][j];
      }
      indast += 1;
    }
  }

  cs_parall_sum(n_structs*3, CS_REAL_TYPE, (cs_real_t *)ms->forstr);

  /* Compute effort sent to internal structures */

  const cs_real_t cfopre = ms->cfopre;

  for (int i = 0; i < n_structs; i++) {
    for (int j = 0; j < 3; j++) {
      ms->forstp[i][j] =          cfopre * ms->forstr[i][j]
                         + (1.0 - cfopre) * ms->forsta[i][j];
    }
  }

  /* Send effort applied to external structures */

  if (n_ast_structs > 0) {
    cs_ast_coupling_exchange_fields();
  }

  /* Structure characteristics defined by the user
     --------------------------------------------- */

  if (n_structs > 0) {

    cs_gui_mobile_mesh_internal_structures(ms->xmstru,
                                           ms->xcstru,
                                           ms->xkstru,
                                           ms->forstp);

    cs_user_fsi_structure_values(n_structs,
                                 ts,
                                 ms->xstreq,
                                 ms->xstr,
                                 ms->xpstr,
                                 ms->xmstru,
                                 ms->xcstru,
                                 ms->xkstru,
                                 ms->forstp,
                                 ms->dtstr);

  }

  /* If the fluid is initializing, we do not read structures */
  if (itrale <= cs_glob_ale_n_ini_f) {
    *itrfin = -1;
    return;
  }

  /* Internal structures displacement
     -------------------------------- */

  for (int i = 0; i < n_structs; i++) {

    _newmark(i, ms->alpnmk, ms->betnmk, ms->gamnmk,
             ms->xmstru[i], ms->xcstru[i], ms->xkstru[i],
             ms->xstreq[i],
             ms->xstr[i], ms->xpstr[i], ms->xppstr[i],
             ms->xsta[i], ms->xpsta[i], ms->xppsta[i],
             ms->forstp[i], ms->forsta[i], ms->dtstr[i]);

  }

  /* Convergence test
     ---------------- */

  int icvext = 0, icvint = 0, icved  = 0;
  double delta = 0.;

  for (int i = 0; i < n_structs; i++) {
    delta += cs_math_3_square_distance(ms->xstr[i], ms->xstp[i]);
  }

  if (n_structs > 0) {
    const cs_real_t almax = cs_glob_turb_ref_values->almax;

    delta = sqrt(delta) / almax / n_structs;
    if (delta < cs_glob_mobile_structures_i_eps)
      icvint = 1;
  }

  if (n_ast_structs > 0)
    icvext = cs_ast_coupling_get_ext_cvg();

  if (n_structs > 0) {
    if (n_ast_structs > 0)
      icved = icvext*icvint;
    else
      icved = icvint;
  }
  else if (n_ast_structs > 0)
    icved = icvext;

  if (eqp->verbosity >= 2)
    bft_printf(_("            Implicit ALE: iter=%5d drift=%12.5e\n"),
               italim, delta);

  /* if converged */
  if (icved == 1) {
    if (*itrfin == 1) {
      /* if itrfin=1 we exit */
      if (eqp->verbosity >= 1)
        bft_printf(_("            Implicit ALE: iter=%5d drift=%12.5e\n"),
                   italim, delta);
      *itrfin = -1;
    }
    else {
      /* Otherwise one last iteration for SYRTHES/T1D/radiation
         and reset icved to 0 so code_aster also runs an iteration;
         FIXME: this can probably be simplified, as "last iteration" for
         Syrthes was only required long ago (for Syrthes 3). */
      *itrfin = 1;
      icved = 0;
    }
  }
  else if (*itrfin == 0  &&  italim == cs_glob_mobile_structures_i_max-1) {
    /* this will be the last iteration */
    *itrfin = 1;
  }
  else if (italim == cs_glob_mobile_structures_i_max) {
    /* we have itrfin=1 and are finished */
    if (cs_glob_mobile_structures_i_max > 1)
      bft_printf(_("@\n"
                   "@ @@ Warning: implicit ALE'\n"
                   "@    =======\n"
                   "@  Maximum number of iterations (%d) reached\n"
                   "@  Normed drift: %12.5e\n"
                   "@\n"),
                 italim, delta);
    *itrfin = -1;
    /* Set icved to 1 so code_aster also stops. */
    icved = 1;
  }

  /* Return the final convergence indicator to code_aster */
  if (n_ast_structs > 0)
    cs_ast_coupling_send_cvg(icved);

  /* Restore previous values if required
     ----------------------------------- */

  /* If nterup  >  1, values at previous time step have been modified
     after cs_solve_navier_stokes; We must then go back to a previous value. */

  if (*itrfin != -1) {

    const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
    const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;

    const int n_fields = cs_field_n_fields();

    for (int field_id = 0; field_id < n_fields; field_id++) {
      cs_field_t *f = cs_field_by_id(field_id);
      if (   f->type & CS_FIELD_VARIABLE
          && f->location_id == CS_MESH_LOCATION_CELLS
          && (f->type & CS_FIELD_CDO) == 0) {

        cs_real_t *cvar_var = f->val;
        cs_real_t *cvara_var = f->val_pre;
        cs_lnum_t n_vals = (cs_lnum_t)f->dim*n_cells_ext;

        if (   f == CS_F_(p)
            && cs_glob_velocity_pressure_param->nterup > 1) {
          cs_real_t *xprale = _pr_save;
          for (cs_lnum_t i = 0; i < n_vals; i++)
            cvara_var[i] = xprale[i];
        }

        for (cs_lnum_t i = 0; i < n_vals; i++)
          cvar_var[i] = cvara_var[i];

      }
    }

    /* Restore mass fluxes */

    int kimasf = cs_field_key_id("inner_mass_flux_id");
    int kbmasf = cs_field_key_id("boundary_mass_flux_id");

    cs_field_t *f_i = cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kimasf));
    cs_field_t *f_b = cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf));

    cs_real_t *imasfl = f_i->val;
    cs_real_t *bmasfl = f_b->val;
    cs_real_t *imasfl_pre = f_i->val_pre;
    cs_real_t *bmasfl_pre = f_b->val_pre;

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
      imasfl[face_id] = imasfl_pre[face_id];

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
      bmasfl[face_id] = bmasfl_pre[face_id];

    /* Restore BC coefficients.

       Using separate values for velocity and pressure could
       also make this more readable and safer. */

    cs_real_3_t *coefau = (cs_real_3_t *)CS_F_(vel)->bc_coeffs->a;
    cs_real_33_t *coefbu = (cs_real_33_t *)CS_F_(vel)->bc_coeffs->b;

    cs_real_t *coefap = CS_F_(p)->bc_coeffs->a;
    cs_real_t *coefbp = CS_F_(p)->bc_coeffs->b;

    _cs_real_11_t *cofale = _bc_coeffs_save;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      coefap[face_id]    = cofale[face_id][0];

      coefau[face_id][0] = cofale[face_id][1];
      coefau[face_id][1] = cofale[face_id][2];
      coefau[face_id][2] = cofale[face_id][3];

      coefbp[face_id]    = cofale[face_id][4];

      coefbu[face_id][0][0] = cofale[face_id][5];
      coefbu[face_id][1][1] = cofale[face_id][6];
      coefbu[face_id][2][2] = cofale[face_id][7];
      coefbu[face_id][1][0] = cofale[face_id][8];
      coefbu[face_id][2][1] = cofale[face_id][9];
      coefbu[face_id][2][0] = cofale[face_id][10];

      /* coefficient B is supposed to be symmetric */
      coefbu[face_id][0][1] = cofale[face_id][8];
      coefbu[face_id][1][2] = cofale[face_id][9];
      coefbu[face_id][0][2] = cofale[face_id][10];

    }

  }

  else if (*itrfin == -1) {
    BFT_FREE(_bc_coeffs_save);
    BFT_FREE(_pr_save);

    if (cs_time_control_is_active(&(ms->plot_time_control), ts)) {

      int t_top_id = cs_timer_stats_switch(_post_out_stat_id);

      _output_time_plots(ms, ts);

      cs_timer_stats_switch(t_top_id);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read mobile structures data to checkpoint.
 *
 * \param[in, out]  r   associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_restart_read(cs_restart_t  *r)
{
  int retcode = CS_RESTART_SUCCESS;
  int n_errors = 0;

  /* When this function is called, the current number of
     coupled structures should be defined */

  cs_mobile_structures_t *ms = _mobile_structures;

  int n_str[2] = {0, 0}, n_str_prev[2] = {0, 0};

  n_str[0] = (ms != nullptr) ? ms->n_structs : 0;
  n_str[1] = cs_glob_ast_coupling_n_couplings;

  if (n_str[0] + n_str[1] == 0)
    return;

  char sec_name[64];
  strcpy(sec_name, "nombre_structures");

  retcode = cs_restart_check_section(r, sec_name, 0, 2, CS_TYPE_int);

  if (retcode == CS_RESTART_ERR_EXISTS)  /* new name */
    strcpy(sec_name, "number_of_mobile_structures");

  retcode = cs_restart_read_section(r,
                                    sec_name,
                                    0, /* location_id */
                                    2,
                                    CS_TYPE_int,
                                    n_str_prev);

  if (retcode != CS_RESTART_SUCCESS) {
    n_str_prev[0] = 0;
    n_str_prev[1] = 0;
  }

  if (n_str_prev[0] > 0 &&
      n_str_prev[0] != n_str[0]) {
    cs_parameters_error
      (CS_ABORT_IMMEDIATE,
       _("Internal mobile structures"),
       _("The number of defined structures is different from the\n"
         "previous calculation:\n"
         "  Number of structures in previous calculation: %d\n"
         "  Number of structures in current calculation: %d\n"
         "\n"
         "Check the coupled boundary structure associations."),
       n_str_prev[0], n_str[0]);
  }
  if (n_str_prev[1] > 0 &&
      n_str_prev[1] != n_str[1]) {
    cs_parameters_error
      (CS_ABORT_IMMEDIATE,
       _("External (code_aster) mobile structures"),
       _("The number of defined structures is different from the\n"
         "previous calculation:\n"
         "  Number of structures in previous calculation: %d\n"
         "  Number of structures in current calculation: %d\n"
         "\n"
         "Check the coupled boundary structure associations."),
       n_str_prev[1], n_str[1]);
  }

  /* Read structure info if present; if we have more structures
     than in the previous run, we assume the first structures
     match and the next ones are added, so we read the available
     data. If we have less structures than previously, we only
     read the only the required data. */

  int n_struct_read = CS_MIN(n_str_prev[0], n_str[0]);

  for (int str_id = 0; str_id < n_struct_read; str_id++) {

    snprintf(sec_name, 63, "donnees_structure_%02d", str_id+1);
    sec_name[63] = '\0';

    retcode = cs_restart_check_section(r, sec_name, 0, 27, CS_TYPE_cs_real_t);
    if (retcode == CS_RESTART_ERR_EXISTS) { /* new name */
      snprintf(sec_name, 63, "mobile_structure_%02d", str_id+1);
      sec_name[63] = '\0';
    }

    cs_real_t tmpstr[27];

    retcode = cs_restart_read_section(r,
                                      sec_name,
                                      0, /* location_id */
                                      27,
                                      CS_TYPE_cs_real_t,
                                      tmpstr);

    if (retcode != CS_RESTART_SUCCESS) {
      n_errors += 1;
      continue;
    }

    for (int i = 0; i < 3; i++) {
      ms->xstr[str_id][i] = tmpstr[i];
      ms->xpstr[str_id][i] = tmpstr[3+i];
      ms->xppstr[str_id][i] = tmpstr[6+i];
      ms->xsta[str_id][i] = tmpstr[9+i];
      ms->xpsta[str_id][i] = tmpstr[12+i];
      ms->xppsta[str_id][i] = tmpstr[15+i];
      ms->xstp[str_id][i] = tmpstr[18+i];
      ms->forstr[str_id][i] = tmpstr[21+i];
      ms->forsta[str_id][i] = tmpstr[24+i];
    }

  }

  /* Abort in case of unhandled error types */
  if (n_errors > 0) {
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: %d error(s) reading mobile structures data\n"
                " in auxiliairy restart file."),
              __func__, n_errors);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write mobile structures data to checkpoint.
 *
 * \param[in, out]  r   associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_restart_write(cs_restart_t  *r)
{
  cs_mobile_structures_t *ms = _mobile_structures;

  int n_str[2] = {0, 0};

  n_str[0] = (ms != nullptr) ? ms->n_structs : 0;
  n_str[1] = cs_glob_ast_coupling_n_couplings;

  if (n_str[0] + n_str[1] == 0)
    return;

  cs_restart_write_section(r,
                           "number_of_mobile_structures",
                           0, /* location_id */
                           2,
                           CS_TYPE_int,
                           n_str);

  for (int str_id = 0; str_id < n_str[0]; str_id++) {

    char sec_name[64];
    snprintf(sec_name, 63, "mobile_structure_%02d", str_id+1);
    sec_name[63] = '\0';

    cs_real_t tmpstr[27];
    for (int i = 0; i < 3; i++) {
      tmpstr[i]   = ms->xstr[str_id][i];
      tmpstr[3+i] = ms->xpstr[str_id][i];
      tmpstr[6+i] = ms->xppstr[str_id][i];
      tmpstr[9+i] = ms->xsta[str_id][i];
      tmpstr[12+i] = ms->xpsta[str_id][i];
      tmpstr[15+i] = ms->xppsta[str_id][i];
      tmpstr[18+i] = ms->xstp[str_id][i];
      tmpstr[21+i] = ms->forstr[str_id][i];
      tmpstr[24+i] = ms->forsta[str_id][i];
    }

    cs_restart_write_section(r,
                             sec_name,
                             0, /* location_id */
                             27,
                             CS_TYPE_cs_real_t,
                             tmpstr);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
