/*============================================================================
 * User definitions for fluid-structure interaction using ALE.
 *============================================================================*/

/* VERS */

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_fluid_structure_interaction-code_aster.c
 *
 * \brief User-defined functions dedicated to Fluid-Structure interaction
 *        modeling. Examples for code_aster.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Definition of internal mobile structures and corresponding initial
 *        conditions (initial displacement and velocity ).
 *
 * \param[in]       is_restart         indicate if computation is restarted
 * \param[in]       n_structs          number of mobile structures
 * \param[in, out]  plot;              monitoring format mask
 *                                       0: no plot
 *                                       1: plot to text (.dat) format
 *                                       2: plot to .csv format
 *                                       3: plot to both formats
 * \param[in, out]  plot_time_control  plot time output frequency control
 * \param[in, out]  aexxst             coefficient for predicted displacement
 * \param[in, out]  bexxst             coefficient for predicted displacement
 * \param[in, out]  cfopre             coefficient for predicted force
 * \param[in, out]  xstr0              initial displacement per structure
 * \param[in, out]  vstr0              initial velocity per structure
 * \param[in, out]  xstreq             displacement of initial mesh relative to
 *                                     structures position at equilibrium
 *
 * \param[in, out]  plot_time_control  time control associated to plotting
 */
/*----------------------------------------------------------------------------*/

void
cs_user_fsi_structure_define(int                 is_restart,
                             int                 n_structs,
                             int                *plot,
                             cs_time_control_t  *plot_time_control,
                             cs_real_t          *aexxst,
                             cs_real_t          *bexxst,
                             cs_real_t          *cfopre,
                             cs_real_t           xstr0[][3],
                             cs_real_t           vstr0[][3],
                             cs_real_t           xstreq[][3])
{
  CS_NO_WARN_IF_UNUSED(n_structs);
  CS_NO_WARN_IF_UNUSED(plot);
  CS_NO_WARN_IF_UNUSED(plot_time_control);
  CS_NO_WARN_IF_UNUSED(aexxst);
  CS_NO_WARN_IF_UNUSED(bexxst);
  CS_NO_WARN_IF_UNUSED(cfopre);
  CS_NO_WARN_IF_UNUSED(xstr0);
  CS_NO_WARN_IF_UNUSED(vstr0);
  CS_NO_WARN_IF_UNUSED(xstreq);

  /*! [fsi_i_ex_a] */
  if (! is_restart) { /* Avoid resetting in case of restart */
    xstr0[0][1]  = 2.0;
    xstreq[0][1] = 1.0;
    vstr0[1][2]  =-0.5;
  }
  /*! [fsi_i_ex_a] */

  /*! [fsi_i_ex_b] */
  *aexxst = 0.5;
  *bexxst = 0.0;
  *cfopre = 2.0;
  /*! [fsi_i_ex_b] */

  /*! [fsi_i_ex_c] */
  *plot = 2;  /* .csv format */

  cs_time_control_init_by_time_step(plot_time_control,
                                    -1,
                                    -1,
                                    5,       /* nt_interval */
                                    true,    /* at_start */
                                    false);  /* at end */
  /*! [fsi_i_ex_c] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Time-based settings for internal mobile structures.
 *
 * \param[in]       n_structs  number of mobile structures
 * \param[in]       ts         time step structure
 * \param[in]       xstreq     displacement of initial mesh rel. to equilibrium
 * \param[in]       xstr       structural displacement
 * \param[in]       vstr       structural velocity
 * \param[in, out]  xmstru     matrix of structural mass
 * \param[in, out]  xcstru     matrix of structural friction
 * \param[in, out]  xkstru     matrix of structural stiffness
 * \param[in, out]  forstr     forces acting on structures (take forces)
 * \param[in, out]  dtstr      structural time step
 */
/*----------------------------------------------------------------------------*/

void
cs_user_fsi_structure_values(int                    n_structs,
                             const cs_time_step_t  *ts,
                             const cs_real_t        xstreq[][3],
                             const cs_real_t        xstr[][3],
                             const cs_real_t        vstr[][3],
                             cs_real_t              xmstru[][3][3],
                             cs_real_t              xcstru[][3][3],
                             cs_real_t              xkstru[][3][3],
                             cs_real_t              forstr[][3],
                             cs_real_t              dtstr[])
{
  CS_NO_WARN_IF_UNUSED(n_structs);
  CS_NO_WARN_IF_UNUSED(ts);
  CS_NO_WARN_IF_UNUSED(xstreq);
  CS_NO_WARN_IF_UNUSED(xstr);
  CS_NO_WARN_IF_UNUSED(vstr);
  CS_NO_WARN_IF_UNUSED(xmstru);
  CS_NO_WARN_IF_UNUSED(xcstru);
  CS_NO_WARN_IF_UNUSED(xkstru);
  CS_NO_WARN_IF_UNUSED(forstr);
  CS_NO_WARN_IF_UNUSED(dtstr);

  /* Definition of structure parameters: mass, friction, stiffness, stresses
     ----------------------------------------------------------------------- */

  /*! [fsi_i_ex_d] */
  /* In the following example structure 0 is defined as an isotropic
   * system (i.e. matrices M, C and K are diagonal): mass equals 5 kg,
   * damping coefficient equals 2 N.s/m and stiffness equals 3 N/m. */

  for (int i = 0; i < 3; i++) {
    xmstru[i][i][0] = 5.0;
    xcstru[i][i][0] = 2.0;
    xkstru[i][i][0] = 3.0;
  }
  /*! [fsi_i_ex_d] */

  /*! [fsi_i_ex_e] */
  /* In this example structure '1' is subjected to the following movement:
   * - In plane xOy the movement is locally defined along an axis (OX).
   *   Structural parameters in X direction are called xm, xc and xk.
   *   The angle of inclination between global (Ox) axis and local (OX) axis
   *   is called theta. Movement in local (OY) direction is imposed to be rigid.
   * - In the 'z' direction the movement is modeled to be oscillating and
   *   harmonic (meaning that there is no friction). Mass equals 1. kg and
   *   stiffness equals 1. N/m. Fluid stresses in that direction are taken
   *   into account. Moreover the structure is also subjected to an external
   *   oscillating force Fz_ext = 3 * cos(4*t).
   *
   *   This leads to the following local equations:
   *     xm.X'' + xc.X' + xk.X = FX
   *                         Y = 0
   *        Z''         +    Z = FZ + 3.cos(4.t)
   */

  double theta = cs_math_pi / 6.0;
  double cost = cos(theta);
  double sint = sin(theta);

  /* FX, FY, and FZ stand for the local fluid forces components.
   * They are defined as follows, using gobal components of
   * fluid forces Fx, Fy and Fz.
   *   fx =  cost*Fx + sint*Fy
   *   fy = -sint*Fx + cost*Fy
   *   fz = Fz
   *
   * After changing of basis, the problem can be described as follows,
   * using global coordinates: */

  double xm = 1.0;
  double xc = 3.e-1;
  double xk = 2.0;
  double fx = forstr[1][0];
  double fy = forstr[1][1];

  xmstru[1][0][0] = xm*cost*cost;
  xmstru[1][1][0] = xm*cost*sint;
  xmstru[1][0][1] = xm*cost*sint;
  xmstru[1][1][1] = xm*sint*sint;
  xmstru[1][2][2] = 1.0;

  xcstru[1][0][0] = xc*cost*cost;
  xcstru[1][1][0] = xc*cost*sint;
  xcstru[1][0][1] = xc*cost*sint;
  xcstru[1][1][1] = xc*sint*sint;

  xkstru[1][0][0] = (xk-1.0)*cost*cost + 1.0;
  xkstru[1][1][0] = (xk-1.0)*cost*sint;
  xkstru[1][0][1] = (xk-1.0)*cost*sint;
  xkstru[1][1][1] = (xk-1.0)*sint*sint + 1.0;
  xkstru[1][2][2] = 1.0;

  forstr[1][0] = fx*cost*cost + fy*sint*cost;
  forstr[1][1] = fx*sint*cost + fy*sint*sint;
  forstr[1][2] += 3.0*cos(4.0*ts->t_cur);

  for (int i = 0; i < n_structs; i++) {
    dtstr[i] = ts->dt[0];
  }
  /*! [fsi_i_ex_e] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define structure numbers for faces associated with internal
 *        or external (code_aster) structures.
 *
 * Structure numbers associated to a given face have the following values:
 * - -i where coupled to  i-th (1-to n) external (code_aster) structure.
 * - 0 where not coupled with an internal or external structure.
 * - i  where coupled to  i-th (1-to n) internal (mass-spring) structure.
 *
 * \param[in, out]  domain         pointer to a cs_domain_t structure
 * \param[in, out]  structure_num  structure id associated to each face
 */
/*----------------------------------------------------------------------------*/

void
cs_user_fsi_structure_num(cs_domain_t  *domain,
                          int           structure_num[])
{
  CS_UNUSED(domain);

  /*! [fsi_i_str_num] */
  int n_structs = 2;
  const char *name[] = {"wall_1", "wall_2"};

  for (int st_id = 0; st_id < n_structs; st_id++) {

    const cs_zone_t  *zn = cs_boundary_zone_by_name(name[st_id]);

    for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

      const cs_lnum_t face_id = zn->elt_ids[e_idx];
      structure_num[face_id] = (st_id + 1);

    }

  }
  /*! [fsi_i_str_num] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
