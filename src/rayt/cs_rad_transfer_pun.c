/*============================================================================
 * Radiation solver boundary conditions treatment.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_math.h"

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

#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_thermal_model.h"
#include "cs_prototypes.h"
#include "cs_face_viscosity.h"
#include "cs_equation_iterative_solve.h"
#include "cs_gradient.h"
#include "cs_face_viscosity.h"

#include "cs_gui_radiative_transfer.h"
#include "cs_rad_transfer.h"
#include "cs_rad_transfer_wall_flux.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_pun.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_rad_transfer_pun.c */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Radiative flux and source term computation
 *
 * \param[in]       bc_type   boundary face types
 * \param[in, out]  coefap    boundary condition work array for the luminance
 *                             (explicit part)
 * \param[in, out]  coefbp    boundary condition work array for the luminance
 *                             (implicit part)
 * \param[in, out]  cofafp    boundary condition work array for the diffusion
 *                             of the luminance (explicit part)
 * \param[in, out]  cofbfp    boundary condition work array for the diffusion
 *                             of the luminance (implicit part)
 * \param[in, out]  flurds    pseudo mass flux work array (interior faces)
 * \param[in, out]  flurdb    pseudo mass flux work array (boundary faces)
 * \param[in, out]  viscf     visc*surface/dist work array at interior faces
 * \param[in, out]  viscb     visc*surface/dist work array at boundary faces
 * \param[in, out]  smbrs     work array for RHS
 * \param[in, out]  rovsdt    work array for unsteady term
 * \param[in]       twall     wall temperature in Kelvin
 * \param[in, out]  ckmel     absorption coefficient for gas-particles mix
 * \param[out]      q         explicit flux density vector
 * \param[in]       abo       weights of the i-th gray gas at boundaries
 * \param[in]       iband     number of the i-th gray gas
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_pun(cs_int_t         bc_type[],
                    cs_real_t        coefap[],
                    cs_real_t        coefbp[],
                    cs_real_t        cofafp[],
                    cs_real_t        cofbfp[],
                    cs_real_t        flurds[],
                    cs_real_t        flurdb[],
                    cs_real_t        viscf[],
                    cs_real_t        viscb[],
                    cs_real_t        smbrs[],
                    cs_real_t        rovsdt[],
                    cs_real_t        twall[],
                    cs_real_t        ckmel[],
                    cs_real_3_t      q[],
                    const cs_real_t  abo[],
                    int              iband)
{
  cs_real_t stephn = cs_physical_constants_stephan;

  /* Pointer to the spectral flux density field */
  cs_field_t *f_qinspe = NULL;
  if (cs_glob_rad_transfer_params->imoadf >= 1)
    f_qinspe = cs_field_by_name_try("spectral_rad_incident_flux");

  cs_field_t *f_qinci = CS_F_(qinci);
  cs_field_t *f_theta4 = CS_FI_(rad_abs, 0);
  cs_field_t *f_thetaa = CS_FI_(rad_emi, 0);
  cs_field_t *f_eps = CS_F_(emissivity);

  cs_real_t *rad_st_expl = CS_FI_(rad_est, 0)->val;

  /* Allocate temporary array  */
  cs_real_t *dpvar;
  BFT_MALLOC(dpvar, cs_glob_mesh->n_cells_with_ghosts, cs_real_t);

  /* Solver settings and initialization */

  cs_var_cal_opt_t  vcopt = cs_parameters_var_cal_opt_default();

  vcopt.imrgra = cs_glob_space_disc->imrgra;
  vcopt.istat  = -1;
  vcopt.ndircl =  1; /* There are Dirichlet BCs  */
  vcopt.isstpc =  0;
  vcopt.iwarni =  cs_glob_rad_transfer_params->iimlum;
  vcopt.blencv =  0.0;
  vcopt.epsrsm =  1e-08;  /* TODO: try with default (1e-07) */
  vcopt.iconv  =  0;      /* No convection for P1 model */
  vcopt.idiff  =  1;      /* Diffusion equation */
  vcopt.idifft = -1;

  int iescap = 0;
  int imucpp = 0;

  /* all boundary convective flux with upwind */
  int icvflb = 0;


  /* Reset arrays before solve */
  for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++) {
    f_theta4->val[iel]    = 0.0;
    f_thetaa->val[iel]    = 0.0;
  }

  for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_i_faces; ifac++)
    flurds[ifac]   = 0.0;

  for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_b_faces; ifac++)
    flurdb[ifac] = 0.0;

  /* Diffusion coefficients at faces  */

  for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++)
    ckmel[iel] = 1.0 / ckmel[iel];

  cs_face_viscosity(cs_glob_mesh,
                    cs_glob_mesh_quantities,
                    cs_glob_space_disc->imvisf,
                    ckmel,
                    viscf,
                    viscb);

  /* Resolution */
  /* Parameter for time scheme and steady case ? */

  cs_equation_iterative_solve_scalar(0,  /* idtvar */
                                     1,  /* external sub-iteration */
                                     -1, /* f_id */
                                     "radiation_p1",
                                     iescap,
                                     imucpp,
                                     -1, /* normp */
                                     &vcopt,
                                     f_thetaa->val,
                                     f_theta4->val,
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
                                     icvflb,
                                     NULL,
                                     rovsdt,
                                     smbrs,
                                     f_theta4->val,
                                     dpvar,
                                     NULL,
                                     NULL);

  /* Radiative flux density Q */

  int inc = 1;
  int iccocg = 1;
  int imligp =  -1;
  int iwarnp = cs_glob_rad_transfer_params->iimlum;
  cs_real_t epsrgp = 1e-08;
  cs_real_t climgp = 1.5;
  cs_real_t extrap = 0.0;
  int nswrgp = 100;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(vcopt.imrgra,
                             &gradient_type,
                             &halo_type);

  int idimtr = 0;
  int hyd_p_flag = 0;

  cs_gradient_scalar("radiative_flux",
                     gradient_type,
                     halo_type,
                     inc,
                     iccocg,
                     nswrgp,
                     idimtr,
                     hyd_p_flag,
                     1,             /* w_stride */
                     iwarnp,
                     imligp,
                     epsrgp,
                     extrap,
                     climgp,
                     NULL,
                     coefap,
                     coefbp,
                     f_theta4->val,
                     NULL,
                     NULL, /* internal coupling */
                     q);

  cs_real_t aa = - stephn * 4.0 / 3.0;
  for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++) {
    cs_real_t aaa = aa * ckmel[iel];
    q[iel][0] = q[iel][0] * aaa;
    q[iel][1] = q[iel][1] * aaa;
    q[iel][2] = q[iel][2] * aaa;
  }

  /* Absorption radiative source temr and incident flux density */

  /* Compute part of absorption or radiative source term */
  aa = 4.0 * stephn;
  for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++)
    rad_st_expl[iel] = aa * f_theta4->val[iel];

  const cs_real_t *b_dist = cs_glob_mesh_quantities->b_dist;

  /*     Calcul du flux incident Qincid  */
  for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_b_faces; ifac++) {
    cs_lnum_t iel  = cs_glob_mesh->b_face_cells[ifac];

    if (   bc_type[ifac] == CS_SMOOTHWALL
        || bc_type[ifac] == CS_ROUGHWALL) {

      if (cs_glob_rad_transfer_params->imoadf >= 1) {
        f_qinspe->val[iband + ifac * f_qinspe->dim] =
            stephn * (  (2.0 * f_theta4->val[iel])
                      + (  abo[ifac + (iband) * cs_glob_mesh->n_b_faces]
                         * f_eps->val[ifac] * (pow (twall[ifac], 4))))
          / (2.0 - f_eps->val[ifac]);
      } else {
        cs_real_t tw4 = pow(twall[ifac], 4.);
        cs_real_t aaa = 1.5*b_dist[ifac]/ckmel[iel]
                        * ( 2. /(2.-f_eps->val[ifac])-1.);
        aa = (aaa*tw4+f_theta4->val[iel])/(1.+aaa);

        f_qinci->val[ifac]
          =  stephn * ( 2.0 * aa - f_eps->val[ifac] * tw4)
                    / (2.0 - f_eps->val[ifac]);
      }
    }

    else {
      if (cs_glob_rad_transfer_params->imoadf >= 1)
        f_qinspe->val[iband + ifac * f_qinspe->dim]
          =   stephn * f_theta4->val[iel]
            + (  q[0][iel] * cs_glob_mesh_quantities->b_face_normal[ifac * 3]
               + q[1][iel] * cs_glob_mesh_quantities->b_face_normal[ifac * 3 + 1]
               + q[2][iel] * cs_glob_mesh_quantities->b_face_normal[ifac * 3 + 2])
            / (0.5 * cs_glob_mesh_quantities->b_face_surf[ifac]);

      else
        f_qinci->val[ifac]
          =   stephn * f_theta4->val[iel]
            + (  q[0][iel] * cs_glob_mesh_quantities->b_face_normal[ifac * 3]
               + q[1][iel] * cs_glob_mesh_quantities->b_face_normal[ifac * 3 + 1 ]
               + q[2][iel] * cs_glob_mesh_quantities->b_face_normal[ifac * 3 + 2])
            / (0.5 *  cs_glob_mesh_quantities->b_face_surf[ifac]);

    }
  }

  /* Free memory     */
  BFT_FREE(dpvar);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
