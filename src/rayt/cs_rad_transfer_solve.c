/*============================================================================
 * Radiation solver main subroutine.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_blas.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_combustion_model.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_gui_util.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_parameters_check.h"
#include "cs_thermal_model.h"
#include "cs_equation_iterative_solve.h"
#include "cs_gradient.h"
#include "cs_face_viscosity.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_time_step.h"

#include "cs_gui_radiative_transfer.h"
#include "cs_rad_transfer.h"
#include "cs_rad_transfer_absorption.h"
#include "cs_rad_transfer_pun.h"
#include "cs_rad_transfer_bcs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_solve.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*  \file cs_rad_transfer_solve.c */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

static int ipadom = 0;

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions for fortran API
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Descend binary tree for the lexicographical ordering of axis coordinates.
 *
 * parameters:
 *   s     <-- axial coordinate
 *   level <-- level of the binary tree to descend
 *   n     <-- number of entities in the binary tree to descend
 *   order <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_axis_descend_tree(const cs_real_t   s[],
                         size_t            level,
                         const size_t      n,
                         cs_lnum_t         order[])
{
  const cs_real_t eps = 1e-24;

  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (n/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      if (s[i1] > s[i2])
        lv_cur++;
      else if (s[i1] + eps > s[i2] && i1 > i2)
        lv_cur++;

    }

    if (lv_cur >= n) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (s[i1] > s[i2])
      break;
    if (s[i1] + eps >= s[i2] && i1 >= i2)
      break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order cells based on an axis coordinate and original numbering.
 *
 * parameters:
 *   s     <-- axial coordinate
 *   order --> pointer to pre-allocated ordering table
 *   n     <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_axis(const cs_real_t  s[],
            cs_lnum_t        order[],
            const size_t     n)
{
  size_t i;
  cs_lnum_t o_save;

  /* Initialize ordering array */

  for (i = 0; i < n; i++)
    order[i] = i;

  if (n < 2)
    return;

  /* Create binary tree */

  i = (n / 2);
  do {
    i--;
    _order_axis_descend_tree(s, i, n, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = n - 1; i > 0; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_axis_descend_tree(s, 0, i, order);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Order linear solvers for DOM radiative model.
 */
/*----------------------------------------------------------------------------*/

static void
_order_by_direction(void)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;

  int kdir = 0;

  cs_real_t *s;
  BFT_MALLOC(s, n_cells, cs_real_t);

  for (int ii = -1; ii < 2; ii+=2) {
    for (int jj = -1; jj < 2; jj+=2) {
      for (int kk = -1; kk < 2; kk+=2) {
        for (int dir_id = 0; dir_id < cs_glob_rad_transfer_params->ndirs; dir_id++) {

          cs_real_t v[3] = {ii*cs_glob_rad_transfer_params->vect_s[dir_id][0],
                            jj*cs_glob_rad_transfer_params->vect_s[dir_id][1],
                            kk*cs_glob_rad_transfer_params->vect_s[dir_id][2]};

          /* Gloal direction id */
          kdir++;

          char name[32];
          sprintf(name, "radiation_%03d", kdir);

          cs_sles_t *sles = cs_sles_find(-1, name);

          if (sles == NULL) {

            if (cs_glob_rad_transfer_params->dispersion == false) {

              cs_sles_it_t *sc = cs_sles_it_define(-1,
                                                   name,
                                                   CS_SLES_P_GAUSS_SEIDEL,
                                                   0,      /* poly_degree */
                                                   1000);  /* n_max_iter */

              for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
                s[c_id] =   v[0]*cell_cen[c_id][0]
                          + v[1]*cell_cen[c_id][1]
                          + v[2]*cell_cen[c_id][2];

              cs_lnum_t *order;
              BFT_MALLOC(order, n_cells, cs_lnum_t);

              _order_axis(s, order, n_cells);

              cs_sles_it_assign_order(sc, &order); /* becomes owner of order */

            }
            else { /* In case of dispersion, Jacobi and Gauss-Seidel
                      usually exhibit quite bad convergence. */

              (void)cs_sles_it_define(-1,
                                      name,
                                      CS_SLES_BICGSTAB,
                                      0,       /* poly_degree */
                                      10000);  /* n_max_iter */

            }

          } /* If linear solver as not already associated */

        }
      }
    }
  }

  BFT_FREE(s);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Radiative flux and source term computation
 *
 * 1/ Luminance data at domain boundaries
 *       (BC: reflection and isotropic emission)
 *                              ->  ->           ->
 * 2/ Compute radiance L(x, s) at point x
 *                          d L
 *    by solving equation : --- = -TK.L +TS
 *                          d S
 *                              ->                o
 *    also expressed as: div (L.S ) = -TK.L + TK.L
 *                            ->   /   -> -> ->
 * 3/ Compute flux densities  Q = /  L(X, S).S domega
 *                               /4.pi
 *                                 /   -> ->
 *        and absorption      Sa= /  L(X, S). domega
 *                               /4.pi
 *    by integration of radiances over solid angles.
 *    Note: it is useful to compute the heating rate
 *    -----
 *                                        /   -> ->  ->  ->
 * 4/ Compute the incident flux Qincid = /  L(X ,S ).S . N domega
 *                                      /->->
 *       ->                            / S.N >0
 *       N fluid to wall normal
 *
 * \param[in]       gg_id     number of the i-th gray gas
 * \param[in]       w_gg      Weights of the i-th gray gas at boundaries
 * \param[in]       tempk     temperature in Kelvin
 * \param[in]       bc_type   boundary face types
 * \param[in, out]  coefap    boundary condition work array for the radiance
 *                             (explicit part)
 * \param[in, out]  coefbp    boundary condition work array for the radiance
 *                             (implicit part)
 * \param[in, out]  cofafp    boundary condition work array for the diffusion
 *                             of the radiance (explicit part)
 * \param[in, out]  cofbfp    boundary condition work array for the diffusion
 *                             of the radiance (implicit part)
 * \param[in, out]  flurds    pseudo mass flux work array (interior faces)
 * \param[in, out]  flurdb    pseudo mass flux work array (boundary faces)
 * \param[in, out]  viscf     visc*surface/dist work array at interior faces
 * \param[in, out]  viscb     visc*surface/dist work array at boundary faces
 * \param[in, out]  rhs       work array for RHS
 * \param[in, out]  rovsdt    work array for unsteady term
 * \param[out]      q         explicit flux density vector
 * \param[out]      int_rad_domega integral of I dOmega
 * \param[out]      int_abso  work array for absorption
 * \param[out]      int_emi   work array for emission
 * \param[out]      int_rad_ist work array for implicit source term
 */
/*----------------------------------------------------------------------------*/

static void
_cs_rad_transfer_sol(int                        gg_id,
                     cs_real_t                  w_gg[],
                     const cs_real_t            tempk[restrict],
                     cs_real_t        *restrict ckg,
                     int                        bc_type[],
                     cs_real_t        *restrict coefap,
                     cs_real_t        *restrict coefbp,
                     cs_real_t        *restrict cofafp,
                     cs_real_t        *restrict cofbfp,
                     cs_real_t        *restrict flurds,
                     cs_real_t        *restrict flurdb,
                     cs_real_t        *restrict viscf,
                     cs_real_t        *restrict viscb,
                     cs_real_t        *restrict rhs,
                     cs_real_t        *restrict rovsdt,
                     cs_real_3_t      *restrict q,
                     cs_real_t        *restrict int_rad_domega,
                     cs_real_t        *restrict int_abso,
                     cs_real_t        *restrict int_emi,
                     cs_real_t        *restrict int_rad_ist)
{
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  cs_lnum_t n_i_faces  = cs_glob_mesh->n_i_faces;
  cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells   = cs_glob_mesh->n_cells;

  cs_real_3_t *b_face_normal = (cs_real_3_t *)cs_glob_mesh_quantities->b_face_normal;
  cs_real_3_t *i_face_normal = (cs_real_3_t *)cs_glob_mesh_quantities->i_face_normal;
  cs_real_t   *b_face_surf = cs_glob_mesh_quantities->b_face_surf;

  const cs_real_t *cell_vol
    = (const cs_real_t *)cs_glob_mesh_quantities->cell_vol;

  const cs_real_t c_stefan = cs_physical_constants_stephan;
  const cs_real_t onedpi  = 1.0 / cs_math_pi;

  cs_field_t *f_qincid = cs_field_by_name("rad_incident_flux");
  cs_field_t *f_snplus = cs_field_by_name("rad_net_flux");

  /* Allocate work arrays */

  cs_real_t *rhs0, *dpvar, *radiance, *radiance_prev;
  cs_real_t *ck_u_d = NULL;
  BFT_MALLOC(rhs0,  n_cells_ext, cs_real_t);
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_t);
  BFT_MALLOC(radiance,    n_cells_ext, cs_real_t);
  BFT_MALLOC(radiance_prev,   n_cells_ext, cs_real_t);

  /* Specific heat capacity of the bulk phase */
  // CAUTION FOR NEPTUNE INTEGRATION HERE
  cs_field_t *f_cp = CS_F_(cp);

  cs_real_t *dcp;
  BFT_MALLOC(dcp, n_cells_ext, cs_real_t);
  if (cs_glob_fluid_properties->icp > 0) {
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      dcp[cell_id] = 1.0 / f_cp->val[cell_id];
  }
  else {
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      dcp[cell_id] = 1.0 / cs_glob_fluid_properties->cp0;
  }

  /* Upwards/Downwards atmospheric integration */
  /* Postprocessing atmospheric upwards and downwards flux */
  cs_field_t  *f_up = NULL, *f_down = NULL;
  cs_real_t *ck_u = NULL;
  cs_real_t *ck_d = NULL;
  if (cs_glob_rad_transfer_params->atmo_ir_absorption) {
    f_up = cs_field_by_name_try("rad_flux_up");//TODO distinguish IR with solar?
    f_down = cs_field_by_name_try("rad_flux_down");
    BFT_MALLOC(ck_u_d,  n_cells_ext, cs_real_t);
    ck_u = cs_field_by_name("rad_absorption_coeff_up")->val;
    ck_d = cs_field_by_name("rad_absorption_coeff_down")->val;
  }

  /* Initialization */

  cs_field_t *f_qinspe;
  if (cs_glob_rad_transfer_params->imoadf >= 1)
    /* Pointer to the spectral flux density field */
    f_qinspe = cs_field_by_name_try("spectral_rad_incident_flux");

  cs_var_cal_opt_t vcopt = cs_parameters_var_cal_opt_default();

  vcopt.iwarni =  cs_glob_rad_transfer_params->iimlum;
  vcopt.iconv  =  1; /* Pure convection */
  vcopt.istat  = -1;
  vcopt.ndircl = 1;/* There are Dirichlet BCs */
  vcopt.idiff  =  0; /* no face diffusion */
  vcopt.idifft = -1;
  vcopt.isstpc =  0;
  vcopt.nswrsm =  1; /* One sweep is sufficient because of the upwind scheme */
  vcopt.imrgra =  cs_glob_space_disc->imrgra;
  vcopt.blencv =  0; /* Pure upwind...*/
  vcopt.epsrsm =  1e-08;  /* TODO: try with default (1e-07) */

  if (cs_glob_rad_transfer_params->dispersion) {
    vcopt.idiff  =  1; /* Added face diffusion */
    vcopt.nswrgr = 20;
    vcopt.nswrsm =  2;
  }

  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_prev + 1)
    _order_by_direction();

  /*                              / -> ->
   * Correct BCs to ensure : pi= /  s. n domega
   *                            /2PI
   */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
    f_snplus->val[face_id] = 0.0;

  cs_real_3_t vect_s;
  cs_real_t domegat, aa;
  for (int ii = -1; ii <= 1; ii+=2) {
    for (int jj = -1; jj <= 1; jj+=2) {
      for (int kk = -1; kk <= 1; kk+=2) {

        for (int dir_id = 0; dir_id < cs_glob_rad_transfer_params->ndirs; dir_id++) {
          vect_s[0] = ii * cs_glob_rad_transfer_params->vect_s[dir_id][0];
          vect_s[1] = jj * cs_glob_rad_transfer_params->vect_s[dir_id][1];
          vect_s[2] = kk * cs_glob_rad_transfer_params->vect_s[dir_id][2];
          domegat = cs_glob_rad_transfer_params->angsol[dir_id];
          for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
            aa = cs_math_3_dot_product(vect_s, b_face_normal[face_id]);
            aa /= b_face_surf[face_id];
            f_snplus->val[face_id] += 0.5 * (-aa + CS_ABS(aa)) * domegat;
          }
        }

      }
    }
  }

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    coefap[face_id] *= cs_math_pi / f_snplus->val[face_id];
    cofafp[face_id] *= cs_math_pi / f_snplus->val[face_id];
  }

  /* initialization for integration in following loops */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    f_qincid->val[face_id] = 0.0;
    f_snplus->val[face_id] = 0.0;
    if (cs_glob_rad_transfer_params->imoadf >= 1)
      f_qinspe->val[gg_id + face_id * cs_glob_rad_transfer_params->nwsgg] = 0.0;

  }

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    int_rad_domega[cell_id] = 0.0;
    int_abso[cell_id] = 0.0;
    int_emi[cell_id] = 0.0;
    int_rad_ist[cell_id] = 0.0;
    q[cell_id][0] = 0.0;
    q[cell_id][1] = 0.0;
    q[cell_id][2] = 0.0;
    if (cs_glob_rad_transfer_params->atmo_ir_absorption) {
      f_up->val[cell_id] = 0.;
      f_down->val[cell_id] = 0.;
    }
  }

  /* Save rhs in buffer, reload at each change of direction */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    rhs0[cell_id] = rhs[cell_id];

  /* rovsdt loaded once only */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    rovsdt[cell_id] = CS_MAX(rovsdt[cell_id], 0.0);

  /* Angular discretization */

  int kdir = 0;

  for (int ii = -1; ii <= 1; ii+=2) {
    for (int jj = -1; jj <= 1; jj+=2) {
      for (int kk = -1; kk <= 1; kk+=2) {

        for (int dir_id = 0; dir_id < cs_glob_rad_transfer_params->ndirs; dir_id++) {
          vect_s[0] = ii * cs_glob_rad_transfer_params->vect_s[dir_id][0];
          vect_s[1] = jj * cs_glob_rad_transfer_params->vect_s[dir_id][1];
          vect_s[2] = kk * cs_glob_rad_transfer_params->vect_s[dir_id][2];
          domegat = cs_glob_rad_transfer_params->angsol[dir_id];
          /* Gloal direction id */
          kdir++;

          /* Update boundary condition coefficients */
          if (cs_glob_rad_transfer_params->atmo_ir_absorption)
            cs_rad_transfer_bc_coeffs(bc_type,
                                      vect_s,
                                      coefap, coefbp,
                                      cofafp, cofbfp,
                                      NULL, /* only usefull for P1 */
                                      w_gg,
                                      gg_id);


          char    cnom[80];
          snprintf(cnom, 80, "%s%03d", "radiation_", kdir);

          /* Spatial discretization */

          /* Explicit source term */


          /* Upwards/Downwards atmospheric integration */
          if (cs_glob_rad_transfer_params->atmo_ir_absorption) {

            for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
              if (cs_math_3_dot_product(cs_glob_physical_constants->gravity,
                                        vect_s) < 0.0)
                ck_u_d[cell_id] =  ck_u[cell_id] * 3./5.;
              else
                ck_u_d[cell_id] =  ck_d[cell_id] * 3./5.;

              rovsdt[cell_id] =  ck_u_d[cell_id] * cell_vol[cell_id];

              rhs[cell_id]  =  ck_u_d[cell_id] * cell_vol[cell_id]
                                 * c_stefan * cs_math_pow4(tempk[cell_id]) * onedpi;

            }
          }
          else {
            for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
              rhs[cell_id] = rhs0[cell_id];
          }

          /* Implicit source term (rovsdt seen above) */

          if (cs_glob_rad_transfer_params->dispersion) {
            const cs_real_t disp_coeff
              = cs_glob_rad_transfer_params->dispersion_coeff;
            const cs_real_t pi = 4.0 * atan(1.);
            const cs_real_t tan_alpha
              =    sqrt(domegat * (4.*pi - domegat))
                 / (2. * pi - domegat);
            const cs_real_t *i_face_surf = cs_glob_mesh_quantities->i_face_surf;

            for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
              viscf[face_id] = disp_coeff * tan_alpha * i_face_surf[face_id];
          }

          else {
            for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
              viscf[face_id] = 0.0;
          }

          for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
            viscb[face_id] = 0.0;

          for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
            radiance[cell_id] = 0.0;
            radiance_prev[cell_id] = 0.0;
          }

          for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
            flurds[face_id] = cs_math_3_dot_product(vect_s, i_face_normal[face_id]);


          for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
            flurdb[face_id] =  cs_math_3_dot_product(vect_s, b_face_normal[face_id]);

          /* Resolution
             ---------- */

          /* In case of a theta-scheme, set theta = 1;
             no relaxation in steady case either */

          cs_equation_iterative_solve_scalar(0,   /* idtvar */
                                             1,   /* external sub-iteration */
                                             -1,  /* f_id */
                                             cnom,
                                             0,   /* iescap */
                                             0,   /* imucpp */
                                             -1,  /* normp */
                                             &vcopt,
                                             radiance_prev,
                                             radiance,
                                             coefap,
                                             coefbp,
                                             cofafp,
                                             cofbfp,
                                             flurds,
                                             flurdb,
                                             viscf,
                                             viscb,
                                             viscf,
                                             viscb,
                                             NULL,
                                             NULL,
                                             NULL,
                                             0, /* icvflb (upwind) */
                                             NULL,
                                             rovsdt,
                                             rhs,
                                             radiance,
                                             dpvar,
                                             NULL,
                                             NULL);

          /* Integration of fluxes and source terms
           * Increment absorption and emission for Atmo on the fly
           * */
          if (cs_glob_rad_transfer_params->atmo_ir_absorption) {

            for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
              aa = radiance[cell_id] * domegat;
              int_rad_domega[cell_id]  += aa;
              /* Absorption */
              int_abso[cell_id] += ck_u_d[cell_id] * aa;
              /* Emmission */
              int_emi[cell_id] -= ck_u_d[cell_id]
                * c_stefan * domegat * onedpi * cs_math_pow4(tempk[cell_id]);

              int_rad_ist[cell_id] -= 4.0 * dcp[cell_id] * ck_u_d[cell_id]
                           * c_stefan * domegat * onedpi * cs_math_pow3(tempk[cell_id]);

              q[cell_id][0] += aa * vect_s[0];
              q[cell_id][1] += aa * vect_s[1];
              q[cell_id][2] += aa * vect_s[2];
            }

          }
          else {

            for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
              aa = radiance[cell_id] * domegat;
              int_rad_domega[cell_id]  += aa;
              q[cell_id][0] += aa * vect_s[0];
              q[cell_id][1] += aa * vect_s[1];
              q[cell_id][2] += aa * vect_s[2];
            }

          }

          /* Flux incident to wall */

          for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
            cs_lnum_t cell_id = cs_glob_mesh->b_face_cells[face_id];
            aa = cs_math_3_dot_product(vect_s, b_face_normal[face_id]);
            aa /= b_face_surf[face_id];
            aa = 0.5 * (aa + CS_ABS(aa)) * domegat;
            f_snplus->val[face_id] += aa;
            if (cs_glob_rad_transfer_params->imoadf >= 1)
              f_qinspe->val[gg_id + face_id * cs_glob_rad_transfer_params->nwsgg]
                += aa * radiance[cell_id];

            else
              f_qincid->val[face_id]
                += aa * radiance[cell_id];

          }

          if (cs_math_3_dot_product(cs_glob_physical_constants->gravity,
                                    vect_s) < 0.0 && f_up != NULL) {
            for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
              f_up->val[cell_id] += radiance[cell_id] * domegat * vect_s[2];//FIXME S.g/||g||
          }
          else if (cs_math_3_dot_product(cs_glob_physical_constants->gravity,
                                         vect_s) > 0.0 && f_down != NULL) {
            for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
              f_down->val[cell_id] += radiance[cell_id] * domegat * vect_s[2];
          }
        }
      }
    }
  } /* End of loop over directions */

  /* Absorption and emission if not atmo */
  if (!cs_glob_rad_transfer_params->atmo_ir_absorption) {

    /* Absorption */
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      int_abso[cell_id] = ckg[cell_id] * int_rad_domega[cell_id];

    /* Emission and implicit ST */
    if (   cs_glob_physical_model_flag[CS_COMBUSTION_3PT] == -1
        && cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == -1) {
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        int_emi[cell_id] = - ckg[cell_id]
                       * 4.0 * c_stefan
                       * cs_math_pow4(tempk[cell_id + n_cells * 0]);

        int_rad_ist[cell_id] = - 16.0 * dcp[cell_id] * ckg[cell_id]
                           * c_stefan * cs_math_pow3(tempk[cell_id + n_cells * 0]);
      }
    } else {
      cs_real_t *cpro_t4m = cs_field_by_name_try("temperature_4")->val;
      cs_real_t *cpro_t3m = cs_field_by_name_try("temperature_3")->val;

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        int_emi[cell_id] = - ckg[cell_id]
                           * 4.0 * c_stefan * cpro_t4m[cell_id];

        int_rad_ist[cell_id] = - 16.0 * dcp[cell_id] * ckg[cell_id]
                           * c_stefan * cpro_t3m[cell_id];
      }
    }

    BFT_FREE(dcp);
  }


#if 0
  /* TODO add clean generation and log of "per day source terms"
     for atmospheric radiative model */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    bft_printf("srad K/day[%d] = %f\n",
               cell_id, rad_estm[cell_id] *-86400.0
               / CS_F_(rho)->val[cell_id]
               / cs_glob_fluid_properties->cp0);
    if (f_up != NULL)
      bft_printf("f_up[%d] = %f\n",
                 cell_id, f_up->val[cell_id]);
    if (f_down != NULL)
      bft_printf("f_down[%d] = %f\n",
                 cell_id, f_down->val[cell_id]);
  }
#endif

  /* Free memory */

  if (ck_u_d != NULL)
    BFT_FREE(ck_u_d);
  BFT_FREE(rhs0);
  BFT_FREE(dpvar);
  BFT_FREE(radiance);
  BFT_FREE(radiance_prev);
}

/*-------------------------------------------------------------------------------*/
/*!
 * \brief Compute the net radiation flux.
 *
 * The density of net radiation flux must be calculated
 * consistently with the boundary conditions of the intensity.
 * The density of net flux is the balance between the radiative
 * emiting part of a boudary face (and not the reflecting one)
 * and the radiative absorbing part.
 *
 * \param[in]   bc_type   boundary face types
 * \param[in]   coefap    boundary condition work array for the radiance
 *                        (explicit part)
 * \param[in]   twall     inside current wall temperature (K)
 * \param[in]   qincid    radiative incident flux  (W/m2)
 * \param[in]   eps       emissivity (>0)
 * \param[out]  net_flux  net flux (W/m2)
 */
/*-------------------------------------------------------------------------------*/

static void
_compute_net_flux(const int        itypfb[],
                  const cs_real_t  coefap[],
                  const cs_real_t  twall[],
                  const cs_real_t  qincid[],
                  const cs_real_t  eps[],
                  cs_real_t        net_flux[])
{
  const cs_real_t c_stefan = cs_physical_constants_stephan;
  cs_real_t  xmissing = -cs_math_big_r * 0.2;

  /* Initializations */

  /* Net flux dendity for the boundary faces
   * The provided examples are sufficient in most of cases.*/

  /* If the boundary conditions given above have been modified
   *   it is necessary to change the way in which density is calculated from
   *   the net radiative flux consistently.*/

  /* The rule is:
   *   the density of net flux is a balance between the emitting energy from a
   *   boundary face (and not the reflecting energy) and the absorbing radiative
   *   energy. Therefore if a wall heats the fluid by radiative transfer, the
   *   net flux is negative */

  for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_b_faces; ifac++) {

    /* Wall faces */
    if (   itypfb[ifac] == CS_SMOOTHWALL
        || itypfb[ifac] == CS_ROUGHWALL)
      net_flux[ifac] = eps[ifac] * (qincid[ifac] - c_stefan * cs_math_pow4(twall[ifac]));

    /* Symmetry   */
    else if (itypfb[ifac] == CS_SYMMETRY)
      net_flux[ifac] = 0.0;

    /* Inlet/Outlet    */
    else if (   itypfb[ifac] == CS_INLET
             || itypfb[ifac] == CS_CONVECTIVE_INLET
             || itypfb[ifac] == CS_OUTLET
             || itypfb[ifac] == CS_FREE_INLET) {
      if (cs_glob_rad_transfer_params->type == CS_RAD_TRANSFER_DOM)
        net_flux[ifac] = qincid[ifac] - cs_math_pi * coefap[ifac];
      else if (cs_glob_rad_transfer_params->type == CS_RAD_TRANSFER_P1)
        net_flux[ifac] = 0.0;
    }

    /* Set "impossible" value for other faces, so that if this
       is not overwritten in \ref cs_user_rad_transfer_net_flux,
       top if there are forgotten faces   */
    else
      net_flux[ifac] = xmissing;

  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the radiative transfer equation.
 *
 *  Two types of method are available:
 *  - Discretes Ordinates Methods (DOM)
 *  - P-1 approximation (only recommended for pulverized coal)
 *
 *  \param[in, out]  bc_type       boundary face types
 *  \param[in]       dt            time step (per cell)
 *  \param[in]       cp2fol        fuel oil liquid CP
 *  \param[in]       cp2ch         pulverized coal CP's
 *  \param[in]       ichcor        pulverized coal indirection
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_solve(int               bc_type[],
                      const cs_real_t   dt[],
                      cs_real_t         cp2fol,
                      const cs_real_t   cp2ch[],
                      const int         ichcor[])
{
  /* Shorter notation */
  cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;

  int nwsgg = rt_params->nwsgg;

  /* Physical constants */
  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
  const cs_real_t c_stefan = cs_physical_constants_stephan;

  /* Mesh params */
  cs_lnum_t n_cells     = cs_glob_mesh->n_cells;
  cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t n_b_faces   = cs_glob_mesh->n_b_faces;
  cs_lnum_t n_i_faces   = cs_glob_mesh->n_i_faces;

  cs_real_3_t *b_face_normal = (cs_real_3_t *)cs_glob_mesh_quantities->b_face_normal;
  const cs_real_t *b_face_surf = cs_glob_mesh_quantities->b_face_surf;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  char fname[80];

  /* Number of passes */
  ipadom++;

  /* Allocate temporary arrays for the radiative equations resolution */
  cs_real_t *viscf, *viscb, *rhs, *rovsdt;
  BFT_MALLOC(viscf,  n_i_faces, cs_real_t);
  BFT_MALLOC(viscb,  n_b_faces, cs_real_t);
  BFT_MALLOC(rhs,  n_cells_ext, cs_real_t);
  BFT_MALLOC(rovsdt, n_cells_ext, cs_real_t);

  /* Allocate specific arrays for the radiative transfer module */
  cs_real_t *tempk, *coefap, *coefbp, *cofafp, *cofbfp, *flurds, *flurdb;
  BFT_MALLOC(tempk,  n_cells_ext * rt_params->nrphas, cs_real_t);
  BFT_MALLOC(coefap, n_b_faces, cs_real_t);
  BFT_MALLOC(coefbp, n_b_faces, cs_real_t);
  BFT_MALLOC(cofafp, n_b_faces, cs_real_t);
  BFT_MALLOC(cofbfp, n_b_faces, cs_real_t);
  BFT_MALLOC(flurds, n_i_faces, cs_real_t);
  BFT_MALLOC(flurdb, n_b_faces, cs_real_t);

  /* Allocate work arrays */
  /* Absorption coeffcient of the bulk phase  */
  cs_real_t *ckmel;
  BFT_MALLOC(ckmel, n_cells_ext, cs_real_t);

  cs_real_t *twall;
  BFT_MALLOC(twall, n_b_faces, cs_real_t);

  /* Map field arrays */
  cs_field_t *f_tempb = CS_F_(t_b);
  cs_field_t *f_qinci = CS_F_(qinci);
  cs_field_t *f_xlam  = CS_F_(xlam);
  cs_field_t *f_epa   = CS_F_(epa);
  cs_field_t *f_eps   = CS_F_(emissivity);
  cs_field_t *f_fnet  = CS_F_(fnet);

  /* ADF model parameters */
  /* Irradiating spectral flux density   */
  cs_field_t *f_qinsp = NULL;
  if (   rt_params->imoadf >= 1
      || rt_params->imfsck == 1)
    f_qinsp = cs_field_by_name("spectral_rad_incident_flux");

  /* Radiation coeffcient kgi and the corresponding weight
     agi of the i-th grey gas
     (the sum over the grey gases is  CS_FI_(rad_cak, 0)->val)*/
  cs_real_t *kgi, *agi;
  BFT_MALLOC(kgi, n_cells_ext * nwsgg, cs_real_t);
  BFT_MALLOC(agi, n_cells_ext * nwsgg, cs_real_t);

  cs_real_t *int_rad_domega;
  BFT_MALLOC(int_rad_domega,  n_cells_ext, cs_real_t);

  /* Flux density components   */
  cs_real_3_t *iqpar;
  BFT_MALLOC(iqpar, n_cells_ext, cs_real_3_t);

  /* Numer of classes for Coal or Fuel combustion */
  int n_classes = 0;
  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0)
    n_classes = cs_glob_combustion_model->coal.nclacp;
  else if (cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] >= 0)
    n_classes = cs_glob_combustion_model->fuel.nclafu;

  /* Irradiating flux density at walls.
     Careful: Should not be confused with qinci */
  cs_real_t *iqpato;
  BFT_MALLOC(iqpato, n_b_faces, cs_real_t);

  /* Weight of the i-th grey gas at walls     */
  cs_real_t *w_gg;
  BFT_MALLOC(w_gg, n_b_faces * nwsgg, cs_real_t);

  /* Wall temperature */
  cs_real_t xptk;
  if (cs_glob_thermal_model->itpscl == 2)
    xptk = tkelvi;
  else
    xptk = 0.0;

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {
    if (   bc_type[ifac] == CS_SMOOTHWALL
        || bc_type[ifac] == CS_ROUGHWALL)
      twall[ifac]  = f_tempb->val[ifac] + xptk;
    else
      twall[ifac]  = 0.0;
  }

  cs_real_t *wq = cs_glob_rad_transfer_params->wq;

  /* FSCK model parameters */
  if (wq == NULL) {
    /* Weight of the i-the gaussian quadrature  */
    BFT_MALLOC(wq, nwsgg, cs_real_t);
    cs_glob_rad_transfer_params->wq = wq;

    /* Must be set to 1 in case of using the standard as well as */
    /* the ADF radiation models  */
    for (int i = 0; i < nwsgg; i++)
      wq[i] = 1.0;
  }


  /* Initializations
     --------------- */

  if (   ipadom > 1
      && cs_glob_time_step->nt_cur%rt_params->nfreqr != 0)
    return;

  cs_log_printf(CS_LOG_DEFAULT,
                _("   ** Information on the radiative source term\n"
                  "      ----------------------------------------\n"));

  /* Constants initialization */
  cs_real_t onedpi  = 1.0 / cs_math_pi;

  /* Bulk absoption:
   * Radiation absorbed by the gasphase and the solid phase
     (all particles classes) */
  cs_real_t *absom = NULL;
  if (n_classes > 0)
    BFT_MALLOC(absom, n_cells_ext, cs_real_t);
  else
    absom = CS_FI_(rad_abs, 0)->val;

  /* Bulk emission by the gasphase and the solid phase
     (all particles classes)*/
  cs_real_t *emim = NULL;
  if (n_classes > 0)
    BFT_MALLOC(emim, n_cells_ext, cs_real_t);
  else
    emim = CS_FI_(rad_emi, 0)->val;

  /* Bulk implicit source term */
  cs_real_t *rad_istm = NULL;
  if (n_classes > 0)
    BFT_MALLOC(rad_istm, n_cells_ext, cs_real_t);
  else
    rad_istm = CS_FI_(rad_ist, 0)->val;

  /* Bulk explicit source term */
  cs_real_t *rad_estm = NULL;
  if (n_classes > 0)
    BFT_MALLOC(rad_estm, n_cells_ext, cs_real_t);
  else
    rad_estm = CS_FI_(rad_est, 0)->val;

  /* Medium (gas) Absorption coefficient */
  cs_real_t *ckg = NULL;
  if (n_classes > 0)
    BFT_MALLOC(ckg, n_cells_ext, cs_real_t);
  else
    ckg = CS_FI_(rad_cak, 0)->val;

  /* Work arays */
  cs_real_t *int_abso, *int_emi, *int_rad_ist;
  BFT_MALLOC(int_abso, n_cells_ext, cs_real_t);
  BFT_MALLOC(int_emi, n_cells_ext, cs_real_t);
  BFT_MALLOC(int_rad_ist, n_cells_ext, cs_real_t);

  cs_real_t *cpro_lumin = CS_F_(rad_lumin)->val;
  cs_real_3_t *cpro_q = (cs_real_3_t *)(CS_F_(rad_q)->val);

  /* Work arrays */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Radiation coefficient k of the gas phase */

    /* Bulk Implicit ST due to emission
     * for the gas phase and the solid/droplet phase (all classes) */
    rad_istm[cell_id] = 0.0;

    /* Explicit ST due to emission and absorption */
    rad_estm[cell_id] = 0.0;

    /* Absortion: Sum, i((kg, i+kp) * Integral(Ii)dOmega):
     * for the gas phase and the solid/droplet phase (all classes) */
    absom[cell_id] = 0.0;

    /* Emmitted radiation: Sum, i((kg, i+kp) * c_stefan * T^4 *agi):
     * for the gas phase and the solid/droplet phase (all classes) */
    emim[cell_id]  = 0.0;

    /* radiative flux vector */
    cpro_q[cell_id][0] = 0.0;
    cpro_q[cell_id][1] = 0.0;
    cpro_q[cell_id][2] = 0.0;

    /* Total emitted intensity   */
    cpro_lumin[cell_id] = 0.0;

    /* Radiation coeffcient of the bulk phase:
     * for the gas phase and the solid/droplet phase (all classes) */
    ckmel[cell_id] = 0.0;

    /* In case of grey gas radiation properties (kgi!=f(lambda))    */
    /* agi must be set to 1.     */
    for (int i = 0; i < nwsgg; i++) {
      kgi[cell_id + n_cells * i]  = 0.0;
      agi[cell_id + n_cells * i]  = 1.0;
    }
  }

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {
    iqpato[ifac] = 0.0;
    /* In case of grey gas radiation properties (kgi!=f(lambda))    */
    /* w_gg must be set to 1.    */
    for (int i = 0; i < nwsgg; i++) {
      w_gg[ifac + i * n_b_faces]     = 1.0;
    }
  }

  for (int class_id = 0; class_id < n_classes; class_id++) {
    int ipcla = 1 + class_id;
    /* Absorbed and emmitted radiation of a single size class
     * (needed to compute the source terms of the particle enthalpy equation) */
    cs_real_t *cpro_abso = CS_FI_(rad_abs, ipcla)->val;
    cs_real_t *cpro_emi  = CS_FI_(rad_emi, ipcla)->val;
    cs_real_t *cpro_tsri = CS_FI_(rad_ist, ipcla)->val;
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      cpro_abso[cell_id] = 0.0;
      cpro_emi[cell_id]  = 0.0;
      cpro_tsri[cell_id] = 0.0;
    }
  }

  /* Store temperature (in Kelvin) in tempk(cell_id, irphas) */

  /* --> Temperature transport */

  if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE) {

    /* val index to access, necessary for compatibility with neptune */
    cs_field_t *temp_field = cs_field_by_name_try("temperature");
    cs_real_t *cvara_scalt;
    if (temp_field != NULL)
      cvara_scalt = temp_field->vals[1];
    else
      cvara_scalt = CS_FI_(t,0)->val;

    /* temperature is in Celsius */
    if (cs_glob_thermal_model->itpscl == 2) {
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        tempk[cell_id] = cvara_scalt[cell_id] + tkelvi;

    }
    else {
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        tempk[cell_id] = cvara_scalt[cell_id];
    }

  }

  /* Enthalpy transport (flurdb is a temporary array) */

  else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

    const cs_real_t *cvara_scalt = CS_F_(h)->vals[1];

    CS_PROCF(c_h_to_t, C_H_TO_T)(cvara_scalt, tempk);

    /* Coal particles or fuel temperature */
    for (int class_id = 0; class_id < n_classes; class_id++) {
      int ipcla = 1 + class_id;
      if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0)
        snprintf(fname, 80, "t_p_%02d", class_id+1);
      else if (cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] >= 0)
        snprintf(fname, 80, "t_fuel_%02d", class_id+1);
      cs_field_t *f_temp2 = cs_field_by_name(fname);
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        tempk[cell_id + n_cells * ipcla] = f_temp2->val[cell_id];
      }
    }

  }
  else
    cs_parameters_error
      (CS_ABORT_IMMEDIATE,
       _("Radiative transfer module"),
       _("Compatible thermal model should be temperature or enthaly-based,\n"
         "but here, \"cs_glob_thermal_model->itherm\" = %d."),
       cs_glob_thermal_model->itherm);

  /* Absorption coefficient for different modules;
     Warning: for the P-1 approximation, the absorption coefficient
     is required for boundary conditions. */

  if (cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] >= 2)
    cs_rad_transfer_absorption(tempk, kgi, agi, w_gg);

  else {

    /* Absorption coefficient;
       initialization to a non-admissible value
       for testing after cs_user_rad_transfer_absorption(). */

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      ckg[cell_id] = -cs_math_big_r;

    /* Data from GUI */

    if (cs_gui_file_is_loaded()) {

      /* FIXME for ADF */

      cs_gui_rad_transfer_absorption(ckg);

      if (   rt_params->type == CS_RAD_TRANSFER_P1
          && cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] <= 1
          && ipadom <= 3)

        cs_rad_transfer_absorption_check_p1(ckg);

    }

    /* Only necessary when grey gas radiation properties are applied.
       In case of the ADF model this test doesnt make sense. */

    if (rt_params->imoadf == 0 && rt_params->imfsck == 0) {

      cs_user_rad_transfer_absorption(bc_type,
                                      dt,
                                      ckg);

      if (cs_glob_rad_transfer_params->type == CS_RAD_TRANSFER_P1)
        cs_rad_transfer_absorption_check_p1(ckg);

    }

  }

  /* -> Test if the radiation coeffcient has been assigned */
  if (rt_params->type > CS_RAD_TRANSFER_NONE) {

    cs_real_t ckmin = 0.0;

    if (   rt_params->imoadf == 0
        && rt_params->imfsck == 0) {

      ckmin = ckg[0];
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        ckmin = CS_MIN(ckmin, ckg[cell_id]);

    }
    else {
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int gg_id = 0; gg_id < nwsgg; gg_id++)
          ckmin = CS_MIN(ckmin, kgi[cell_id + n_cells * gg_id]);
      }

    }

    cs_parall_min(1, CS_DOUBLE, &ckmin);

    if (ckmin < 0.0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s (in %s)\n"
           "The absorption coefficient must be > 0, but here the\n"
           "minimal value encountered is %g."),
         _("Radiative transfer module:\n"
           "-------------------------\n"), __func__,
         ckmin);

  }

  /* Check for transparent case => no need to compute absoption or emission */
  int idiver = rt_params->idiver;

  /* Solving the ETR.
     Loop over all gray gases. In case of the basic radiation models
     of Code_Saturne, nwsgg=1 */

  cs_real_t *cpro_cak;

  for (int gg_id = 0; gg_id < nwsgg; gg_id++) {

    if (   rt_params->imoadf >= 1
        || rt_params->imfsck == 1) {

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        ckg[cell_id] = kgi[cell_id + n_cells * gg_id];

    }
    else if (!cs_glob_rad_transfer_params->atmo_ir_absorption) {
      cs_real_t ckmax = 0.0;

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        ckmax = CS_MAX(ckmax, ckg[cell_id]);

      cs_parall_max(1, CS_REAL_TYPE, &ckmax);

      if (ckmax <= cs_math_epzero) {
        cs_log_printf(CS_LOG_DEFAULT,
                      _("      Radiative transfer with transparent medium."));
        idiver = -1;
      }
    }
    /* Infra red atmospheric model */
    else {
      cs_real_t *ck_u = cs_field_by_name("rad_absorption_coeff_up")->val;
      cs_real_t *ck_d = cs_field_by_name("rad_absorption_coeff_down")->val;

      cs_real_t ckumax = 0.0;

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        ckumax = CS_MAX(ckumax, ck_u[cell_id]);

      cs_parall_max(1, CS_REAL_TYPE, &ckumax);

      if (ckumax <= cs_math_epzero)
        cs_log_printf(CS_LOG_DEFAULT,
            _("      Atmo radiative transfer with transparent medium for upward directions.\n"));

      cs_real_t ckdmax = 0.0;

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        ckdmax = CS_MAX(ckdmax, ck_d[cell_id]);

      cs_parall_max(1, CS_REAL_TYPE, &ckdmax);

      if (ckdmax <= cs_math_epzero)
        cs_log_printf(CS_LOG_DEFAULT,
            _("      Atmo radiative transfer with transparent medium for downward directions.\n"));

      if (ckumax <= cs_math_epzero && ckdmax <= cs_math_epzero)
        idiver = -1;

    }

    /* P-1 radiation model
       ------------------- */

    if (rt_params->type == CS_RAD_TRANSFER_P1) {

      /* Gas phase: Explicit source term in the transport of theta4 */

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        rhs[cell_id] =  3.0 * ckg[cell_id]
                            * cs_math_pow4(tempk[cell_id])
                            * agi[cell_id + n_cells * gg_id]
                            * cell_vol[cell_id];

      /* Solid phase/coal particles or Fuel droplets:
         Explicit source term in the transport eqn. of theta4 */
      for (int class_id = 0; class_id < n_classes; class_id++) {

        int ipcla = class_id + 1;
        cpro_cak = CS_FI_(rad_cak, ipcla)->val;
        snprintf(fname, 80, "x_p_%02d", class_id+1);
        cs_field_t *f_x2 = cs_field_by_name(fname);

        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
          rhs[cell_id] +=  3.0 * f_x2->val[cell_id] * cpro_cak[cell_id]
                               * cs_math_pow4(tempk[cell_id + n_cells * ipcla])
                               * agi[cell_id + n_cells * gg_id]
                               * cell_vol[cell_id];
      }

      /* -> Gas phase: Implicit source term in the transport eqn. of theta4 */
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        rovsdt[cell_id] =  3.0 * ckg[cell_id]
                               * cell_vol[cell_id];

      /* -> Coal solid phase or fuel droplets: */
      /* Implicit source term in the transport eqn. of theta4   */
      for (int class_id = 0; class_id < n_classes; class_id++) {

        cpro_cak = CS_FI_(rad_cak, class_id+1)->val;

        snprintf(fname, 80, "x_p_%02d", class_id + 1);
        cs_field_t *f_x2 = cs_field_by_name(fname);

        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
          rovsdt[cell_id] +=  3.0 * f_x2->val[cell_id] * cpro_cak[cell_id]
                                  * cell_vol[cell_id];
      }

      /* Radiation coeffcient of the bulk phase   */
      /* Gas phase: */
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        ckmel[cell_id] = ckg[cell_id];

      /* Coal solid phase or fuel droplets */
      for (int class_id = 0; class_id < n_classes; class_id++) {

        cpro_cak = CS_FI_(rad_cak, class_id+1)->val;

        snprintf(fname, 80, "x_p_%02d", class_id + 1);
        cs_field_t *f_x2 = cs_field_by_name(fname);

        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
          ckmel[cell_id] += f_x2->val[cell_id] * cpro_cak[cell_id];
      }

      /* Test if ckmel is gt zero  */
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        if (ckmel[cell_id] <= 0.0)
          bft_error
            (__FILE__, __LINE__, 0,
             _("%s (in %s)\n"
               "The local radiation coeffcient of the bulk phase ckmel\n"
               "takes the value 0 somewhere. This often occurs during\n"
               "the very first iterations of the simulation.\n"
               "To avoid this, ensure the coal and/or the char mass fraction\n"
               "have been initialized to values different from zero."),
             _("Radiative transfer module (P-1 radiation):\n"
               "-------------------------\n"),
             __func__);
      }

      /* Update Boundary condition coefficients */
      cs_rad_transfer_bc_coeffs(bc_type,
                                NULL, /* No specific direction */
                                coefap, coefbp,
                                cofafp, cofbfp,
                                ckmel,
                                w_gg,   gg_id);

      /* Solving */
      cs_rad_transfer_pun(bc_type,
                          coefap, coefbp,
                          cofafp, cofbfp,
                          flurds, flurdb,
                          viscf, viscb,
                          rhs, rovsdt,
                          twall, ckmel,
                          iqpar,
                          w_gg,
                          int_rad_domega,
                          gg_id);

      /* Precomputed absoption and emission */
      /* Absorption */
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        int_abso[cell_id] = ckg[cell_id] * int_rad_domega[cell_id];

      /* Specific heat capacity of the bulk phase */
      // CAUTION FOR NEPTUNE INTEGRATION HERE
      cs_field_t *f_cp = CS_F_(cp);

      cs_real_t *dcp;
      BFT_MALLOC(dcp, n_cells_ext, cs_real_t);
      if (cs_glob_fluid_properties->icp > 0) {
        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
          dcp[cell_id] = 1.0 / f_cp->val[cell_id];
      }
      else {
        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
          dcp[cell_id] = 1.0 / cs_glob_fluid_properties->cp0;
      }

      /* Emission and implicit ST */
      if (   cs_glob_physical_model_flag[CS_COMBUSTION_3PT] == -1
          && cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == -1) {
        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
          int_emi[cell_id] = - ckg[cell_id]
            * 4.0 * c_stefan
            * cs_math_pow4(tempk[cell_id + n_cells * 0]);

          int_rad_ist[cell_id] = - 16.0 * dcp[cell_id] * ckg[cell_id]
            * c_stefan * cs_math_pow3(tempk[cell_id + n_cells * 0]);
        }
      } else {
        cs_real_t *cpro_t4m = cs_field_by_name_try("temperature_4")->val;
        cs_real_t *cpro_t3m = cs_field_by_name_try("temperature_3")->val;

        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
          int_emi[cell_id] = - ckg[cell_id]
            * 4.0 * c_stefan * cpro_t4m[cell_id];

          int_rad_ist[cell_id] = - 16.0 * dcp[cell_id] * ckg[cell_id]
            * c_stefan * cpro_t3m[cell_id];
        }
      }

      BFT_FREE(dcp);

    }

    /* Solving of the radiative transfer equation (DOM)
       ------------------------------------------------ */

    else if (rt_params->type == CS_RAD_TRANSFER_DOM) {

      /* -> Gas phase: Explicit source term of the ETR */

      if (   cs_glob_physical_model_flag[CS_COMBUSTION_3PT] == -1
          && cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == -1) {
        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
          rhs[cell_id] =  c_stefan * ckg[cell_id]
                                     * cs_math_pow4(tempk[cell_id])
                                     * agi[cell_id + n_cells * gg_id]
                                     * cell_vol[cell_id]
                                     * onedpi;
      } else {
        cs_real_t *cpro_t4m = cs_field_by_name_try("temperature_4")->val;

        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
          rhs[cell_id] =  c_stefan * ckg[cell_id]
                                     * cpro_t4m[cell_id]
                                     * agi[cell_id + n_cells * gg_id]
                                     * cell_vol[cell_id]
                                     * onedpi;
      }

      /* -> Solid phase: */
      /* Coal particles or fuel droplets:
       * Explicit source term of the ETR    */
      for (int class_id = 0; class_id < n_classes; class_id++) {

        int ipcla = class_id + 1;
        cpro_cak = CS_FI_(rad_cak, ipcla)->val;
        snprintf(fname, 80, "x_p_%02d", class_id+1);
        cs_field_t *f_x2 = cs_field_by_name(fname);

        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
          rhs[cell_id] +=   f_x2->val[cell_id]
                            * agi[cell_id + n_cells * gg_id]
                            * c_stefan
                            * cpro_cak[cell_id]
                            * cs_math_pow4(tempk[cell_id + n_cells * ipcla])
                            * cell_vol[cell_id]
                            * onedpi;
      }

      /* -> Gas phase: Implicit source term of the ETR */
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        rovsdt[cell_id] = ckg[cell_id] * cell_vol[cell_id];

      /* -> Coal solid phase or fuel droplets:
       * Implicit source term of the ETR    */
      for (int class_id = 0; class_id < n_classes; class_id++) {

        cpro_cak = CS_FI_(rad_cak, class_id+1)->val;

        snprintf(fname, 80, "x_p_%02d", class_id + 1);
        cs_field_t *f_x2 = cs_field_by_name(fname);

        for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
          rovsdt[cell_id] +=   f_x2->val[cell_id]
                             * cpro_cak[cell_id]
                             * cell_vol[cell_id];
      }

      /* Update boundary condition coefficients:
       * default ones, identical for each directions, may be overwritten
       * afterwards */
      cs_rad_transfer_bc_coeffs(bc_type,
                                NULL, /*no specific direction */
                                coefap, coefbp,
                                cofafp, cofbfp,
                                ckmel,
                                w_gg  , gg_id);

      /* Solving */
      _cs_rad_transfer_sol(gg_id,
                           w_gg,
                           tempk,
                           ckg,
                           bc_type,
                           coefap, coefbp,
                           cofafp, cofbfp,
                           flurds, flurdb,
                           viscf, viscb,
                           rhs, rovsdt,
                           iqpar,
                           int_rad_domega,
                           int_abso,
                           int_emi,
                           int_rad_ist);

    }

    /* Summing up the quantities of each grey gas    */

    /* Absorption */
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      absom[cell_id] += ckg[cell_id] * int_rad_domega[cell_id] * wq[gg_id];

    /* Coal solid phase or fuel droplets */
    for (int class_id = 0; class_id < n_classes; class_id++) {

      int ipcla = class_id + 1;
      cpro_cak = CS_FI_(rad_cak, ipcla)->val;
      /* Absorbed and emmitted radiation of a single size class */
      cs_real_t *cpro_abso = CS_FI_(rad_abs, ipcla)->val;

      snprintf(fname, 80, "x_p_%02d", ipcla);
      cs_field_t *f_x2 = cs_field_by_name(fname);
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        /* Absoprtion of particles is added to absom */
        absom[cell_id] += f_x2->val[cell_id] * cpro_cak[cell_id]
                      * int_rad_domega[cell_id] * wq[gg_id];
        cpro_abso[cell_id]
          +=  cpro_cak[cell_id] * int_rad_domega[cell_id] * wq[gg_id];
      }
    }

    /* Emission (gas phase, precomputed)  */

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      emim[cell_id] += int_emi[cell_id]
        * agi[cell_id + n_cells * gg_id]
        * wq[gg_id];

      rad_istm[cell_id] += int_rad_ist[cell_id]
        * agi[cell_id + gg_id * n_cells]
        * wq[gg_id];
    }

    /* Coal solid phase or fuel droplets */
    for (int class_id = 0; class_id < n_classes; class_id++) {

      int ipcla = class_id + 1;
      cpro_cak = CS_FI_(rad_cak, ipcla)->val;

      /* emmitted radiation of a single size class */
      cs_real_t *cpro_emi  = CS_FI_(rad_emi, ipcla)->val;
      cs_real_t *cpro_tsri = CS_FI_(rad_ist, ipcla)->val;
      snprintf(fname, 80, "x_p_%02d", class_id+1);
      cs_field_t *f_x2 = cs_field_by_name(fname);

      cs_real_t cp2 = 0.;
      if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0)
        cp2 = cp2ch[ichcor[class_id]-1];
      else if (cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] >= 0)
       cp2 = cp2fol;

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

        cs_real_t sig_ck_t4 = 4. * c_stefan * cpro_cak[cell_id]
                            * cs_math_pow4(tempk[cell_id + n_cells * ipcla])
                            * agi[cell_id + n_cells * gg_id]
                            * wq[gg_id];
        cs_real_t sig_ck_t3dcp2 = 16. * c_stefan * cpro_cak[cell_id]
                            * cs_math_pow3(tempk[cell_id + n_cells * ipcla])
                            * agi[cell_id + n_cells * gg_id]
                            * wq[gg_id] / cp2;

        /* Add Emission of particles: kp * c_stefan * T^4 *agi
         * to emim */
        emim[cell_id] -= sig_ck_t4 * f_x2->val[cell_id];

        cpro_emi[cell_id] -= sig_ck_t4;

        /* Implicit ST of the solid phase is added to rad_istm */
        rad_istm[cell_id] -= sig_ck_t3dcp2 * f_x2->val[cell_id];

        cpro_tsri[cell_id] -= sig_ck_t3dcp2;

      }
    }

    /* Storing of the total emitted intensity:
   *      / -> ->
   * SA= / L( X, S ). DOMEGA
   *    /4.PI
   */

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      /* Emitted intensity    */
      cpro_lumin[cell_id]  += (int_rad_domega[cell_id] * wq[gg_id]);

      /* Flux vector components    */
      cpro_q[cell_id][0] += iqpar[cell_id][0] * wq[gg_id];
      cpro_q[cell_id][1] += iqpar[cell_id][1] * wq[gg_id];
      cpro_q[cell_id][2] += iqpar[cell_id][2] * wq[gg_id];
    }

    /* If the ADF model is activated we have to sum
       the spectral flux densities */

    if (rt_params->imoadf >= 1) {
      for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++)
        iqpato[ifac] += f_qinsp->val[gg_id + ifac * nwsgg] * wq[gg_id];
    }

  } /* end loop on grey gas */

  BFT_FREE(int_rad_domega);

  /* The total radiative flux is copied in bqinci   */
  /* a) for post-processing reasons and   */
  /* b) in order to calculate bfnet  */
  if (rt_params->imoadf >= 1) {
    for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++)
      f_qinci->val[ifac] = iqpato[ifac];
  }

  /* Net radiative flux at walls: computation and integration */

  /* -> Initialization to a non-admissible value for testing after
   *    cs_user_rad_transfer_net_flux   */
  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++)
    f_fnet->val[ifac] = -cs_math_big_r;

  /* Basic definition for net flux */

  _compute_net_flux(bc_type,
                    coefap,
                    twall,
                    f_qinci->val,
                    f_eps->val,
                    f_fnet->val);

  /*---> Reading of User data
   * CAREFUL: The user has access to the radiation coeffcient (field f_cak1)
   * in cs_user_rad_transfer_net_flux.
   * However, only when the standard radiation models of code_saturne are
   * applied, this table contains the true value given by the user.
   * Thus, the radiation coefficient in cs_user_rad_transfer_net_flux must
   * be used with caution.
   * In its present version cs_user_rad_transfer_net_flux does NOT use the
   * radiation coeffcient, and thus, cs_user_rad_transfer_net_flux can
   * still be called here, even if the ADF model is activated.
   */

  cs_user_rad_transfer_net_flux(bc_type,
                                dt,
                                coefap,
                                coefbp,
                                cofafp,
                                cofbfp,
                                twall,
                                f_qinci->val,
                                f_xlam->val,
                                f_epa->val,
                                f_eps->val,
                                ckg,
                                f_fnet->val);

  /* Check net flux */
  cs_real_t xlimit = -cs_math_big_r * 0.1;

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {
    if (f_fnet->val[ifac] <= xlimit)
      bc_type[ifac] = - CS_ABS(bc_type[ifac]);
  }

  cs_boundary_conditions_error(bc_type, "Net flux BC values");

  /* Integrate net flux net on different boundary zones */

  cs_boundary_zone_update_face_class_id();
  const int n_zones = cs_boundary_zone_max_class_or_zone_id() + 1;
  const int *b_face_class_id = cs_boundary_zone_face_class_or_zone_id();

  int *iflux;
  BFT_MALLOC(iflux, n_zones, int);

  cs_real_t *flux;
  BFT_MALLOC(flux, n_zones, cs_real_t);

  for (int izone = 0; izone < n_zones; izone++) {
    flux[izone]  = 0.0;
    iflux[izone] = 0;
  }

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {
    int izone = b_face_class_id[ifac];
    flux[izone]  = flux[izone] + f_fnet->val[ifac] * b_face_surf[ifac];
    iflux[izone] = 1;
  }

  if (cs_glob_rank_id >= 0) {
    cs_parall_sum(n_zones, CS_REAL_TYPE, flux);
    cs_parall_max(n_zones, CS_INT_TYPE, iflux);
  }

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n"
       "Zone         Radiative net flux (Watt) (outward-facing unit normal)\n"
       "--------------------------------------\n"));

  for (int izone = 0; izone < n_zones; izone++) {

    if (iflux[izone] == 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    _("%6d             %11.4e\n"),
                    izone,
                    flux[izone]);
  }
  cs_log_printf(CS_LOG_DEFAULT, "\n");

  /* -> Integration de la densite de flux net aux frontieres */
  cs_real_t aa = cs_dot(n_b_faces, f_fnet->val, b_face_surf);

  cs_parall_sum(1, CS_REAL_TYPE, &aa);

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("Net radiative flux on all boundaries:      Fnet = %11.4e Watt\n"),
     aa);

  /* Semi-analitical radiative source terms */
  if (idiver >= 0) {

    /* Emission + Absorption of gas and particles --> explicit ST */
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      rad_estm[cell_id] = absom[cell_id] + emim[cell_id];


    /* In order to determine the source terms of the particle enthalpy
       tranport equation, we have to copy the approriate determined above
       into the corressponding tables. */

    /* Coal solid phase or fuel droplets */
    for (int class_id = 0; class_id < n_classes; class_id++) {
      int ipcla = 1 + class_id;

      cs_real_t *cpro_tsre = CS_FI_(rad_est, ipcla)->val;
      cs_real_t *cpro_abso = CS_FI_(rad_abs, ipcla)->val;
      cs_real_t *cpro_emi  = CS_FI_(rad_emi, ipcla)->val;

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        cpro_tsre[cell_id] = cpro_abso[cell_id]
                             + cpro_emi[cell_id];
      }

    }

  }
  else {

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      absom[cell_id]  = 0.0;
      emim[cell_id]   = 0.0;
      rad_estm[cell_id] = 0.0;
      rad_istm[cell_id] = 0.0;
    }

  }

  /* Explicit conservative radiative source terms
   * ---------------------------------------------*/

  if (idiver == 1 || idiver == 2) {

    /* Allocate temporary arrays for gradient computation */
    cs_real_3_t *coefaq;
    BFT_MALLOC(coefaq, n_b_faces, cs_real_3_t);
    cs_real_33_t *coefbq;
    BFT_MALLOC(coefbq, n_b_faces, cs_real_33_t);

    for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {
      coefaq[ifac][0] = f_fnet->val[ifac] * b_face_normal[ifac][0] / b_face_surf[ifac];
      coefaq[ifac][1] = f_fnet->val[ifac] * b_face_normal[ifac][1] / b_face_surf[ifac];
      coefaq[ifac][2] = f_fnet->val[ifac] * b_face_normal[ifac][2] / b_face_surf[ifac];
    }

    for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
          coefbq[ifac][jsou][isou] = 0.;
      }
    }

    cs_real_33_t *grad;
    BFT_MALLOC(grad, n_cells_ext, cs_real_33_t);

    /* Data for computation of divergence */

    cs_halo_type_t halo_type = CS_HALO_STANDARD;
    cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

    cs_gradient_type_by_imrgra(cs_glob_space_disc->imrgra,
                               &gradient_type,
                               &halo_type);

    cs_gradient_vector("Work array",
                       gradient_type,
                       halo_type,
                       1,      /* inc */
                       100,    /* n_r_sweeps, */
                       rt_params->iimlum,  /* iwarnp */
                       -1,     /* imligp */
                       1e-8,   /* epsrgp */
                       1.5,    /* climgp */
                       (const cs_real_3_t *)coefaq,
                       (const cs_real_33_t *)coefbq,
                       cpro_q,
                       NULL, /* weighted gradient */
                       NULL, /* coupling */
                       grad);

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      rad_estm[cell_id] = - grad[cell_id][0][0]
                          - grad[cell_id][1][1]
                          - grad[cell_id][2][2];
    }

    /* Free memory */
    BFT_FREE(grad);
    BFT_FREE(coefbq);
    BFT_FREE(coefaq);

  } /* End of computation of divergence */

  /* Explicit radiative semi-analytical corrected source term */
  if (idiver == 2) {

    /* Comparison of the semi-analytical and conservative source terms */

    cs_real_t s[2] = {0, 0};

    s[0] = cs_dot(n_cells, rad_estm, cell_vol);

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      s[1] +=   (absom[cell_id] + emim[cell_id])
              * cell_vol[cell_id];

    cs_parall_sum(2, CS_REAL_TYPE, s);

    s[0] /= s[1];

    /* Correction of the semi-analytical source term
       by the conservative source term  */
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      rad_estm[cell_id] = (  absom[cell_id]
                           + emim[cell_id]) * s[0];

  }

  /* Log information */
  if (idiver >= 0) {

    /* -> Integration volumique du terme source explicite
     * Le resultat de cette integration DOIT etre le meme que l'integration
     * surfacique de la densite de flux net radiatif faite plus haut
     * si idiver = 1 ou 2 */

    cs_real_t balance = 0.0;
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      balance += rad_estm[cell_id] * cell_vol[cell_id];

    cs_parall_sum(1, CS_REAL_TYPE, &balance);

    cs_log_printf
      (CS_LOG_DEFAULT,
       _("Volume integral of radiative source term:  Srad = %11.4e Watt\n"
         "(If cs_glob_rad_transfer_params->idiver = 1 or 2, we must have Srad = -Fnet)\n"),
       balance);

    /* Correction of explicit source term in raysca to allow a correct
     * post-processing of that term when the transported variable is the
     * temperature (for combustion, it is always enthalpy) */

  }

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("-------------------------------------------------------------------\n"));

  /* Free memory */

  BFT_FREE(iflux);
  BFT_FREE(flux);
  BFT_FREE(iqpato);
  BFT_FREE(viscf);
  BFT_FREE(viscb);
  BFT_FREE(rhs);
  BFT_FREE(rovsdt);
  BFT_FREE(tempk);
  BFT_FREE(coefap);
  BFT_FREE(coefbp);
  BFT_FREE(cofafp);
  BFT_FREE(cofbfp);
  BFT_FREE(flurds);
  BFT_FREE(flurdb);
  BFT_FREE(int_abso);
  BFT_FREE(int_emi);
  BFT_FREE(int_rad_ist);
  BFT_FREE(ckmel);
  BFT_FREE(twall);
  BFT_FREE(kgi);
  BFT_FREE(agi);
  BFT_FREE(w_gg);
  BFT_FREE(iqpar);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
