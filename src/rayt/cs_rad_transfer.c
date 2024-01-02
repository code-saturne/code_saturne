/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include "cs_boundary_zone.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_rad_transfer.c
           Radiation solver operations.
*/

/*----------------------------------------------------------------------------*/

/*!

  \struct cs_rad_transfer_params_t

  \brief Structure containing the radiation module parameters.

 \var  cs_rad_transfer_params_t::nrphas
        Phase which radiates (bulk by default, but may be coal class or fuel
        droplets phase).
  \var  cs_rad_transfer_params_t::iimpar
        Verbosity level in the log concerning the calculation of
        the wall temperatures:
        - 0: no display
        - 1: standard
        - 2: complete
  \var  cs_rad_transfer_params_t::verbosity
        Verbosity level in the log concerning the calculation of
        the radiative transfer equation:
        - 0: no display
        - 1: standard
        - 2: complete
  \var  cs_rad_transfer_params_t::imodak
        When gas or coal combustion is activated, \ref imodak indicates whether
        the absorption coefficient shall be calculated "automatically" (=1)
        or read from the data file (=0).
  \var  cs_rad_transfer_params_t::imoadf
        ADF model:
        - 0 no ADF model
        - 1 ADF model with 8 intervals of wave length
        - 2 ADF model with 50 intervals of wave length
  \var  cs_rad_transfer_params_t::iwrp1t
        P1 model transparency warnings counter.
  \var  cs_rad_transfer_params_t::imfsck
        FSCK model:
        - 0 no FSCK model
        - 1 FSCK model activated
  \var  cs_rad_transfer_params_t::xnp1mx
        For the P-1 model, percentage of cells for which we allow the optical
        thickness to exceed unity, although this should be avoided.
        (more precisely, where \f$ KL \f$ is lower than 1, where \f$ K \f$ is
        the absorption coefficient of the medium and \f$ L \f$ is a
        characteristic length of the domain).
  \var  cs_rad_transfer_params_t::idiver
        Indicates the method used to calculate the radiative source term:
        - -1: no renormalization
        - 0: semi-analytic calculation (compulsory with transparent media)
        - 1: conservative calculation
        - 2: semi-analytic calculation corrected in order to be globally
        conservative
  \var  cs_rad_transfer_params_t::i_quadrature
        Index of the quadrature and number of directions for a single octant.\n
        Sn quadrature (n(n+2) directions)
        - 1: S4 (24 directions)
        - 2: S6 (48 directions)
        - 3: S8 (80 directions)\n
        Tn quadrature (8n^2 directions)
        - 4: T2 (32 directions)
        - 5: T4 (128 directions)
        - 6: Tn (8*ndirec^2 directions)
        - 7: 120 directions (LC11)
        - 8: 48 directions (DCT020-2468)
  \var  cs_rad_transfer_params_t::ndirec
        Number of directions for the angular discretisation of the radiation
        propagation with the DOM model.\n
        No other possible value, because of the way the directions are calculated.\n
        The calculation with 32 directions may break the symmetry of physically
        axi-symmetric cases (but the cost in CPU time is much lower than with
        128 directions).\n
        Useful if and only if the radiation module is activated with the DOM method.
  \var  cs_rad_transfer_params_t::ndirs
        For the Tn quadrature, \ref ndirec squared
  \var  cs_rad_transfer_params_t::vect_s
        Direction vectors of angular values of the quadrature sx, sy, sz.
  \var  cs_rad_transfer_params_t::angsol
        Weight of the solid angle.
  \var  cs_rad_transfer_params_t::restart
        Indicates whether the radiation variables should be initialized or
        read from a restart file.
  \var  cs_rad_transfer_params_t::nwsgg
        Spectral radiation models (ADF and FSCK).\n
        Number of ETRs to solve.
  \var  cs_rad_transfer_params_t::wq
        Weights of the Gaussian quadrature
  \var  cs_rad_transfer_params_t::itpimp
        \deprecated See CS_BOUNDARY_RAD_WALL_GRAY.
  \var  cs_rad_transfer_params_t::ipgrno
        \deprecated See CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T.
  \var  cs_rad_transfer_params_t::iprefl
        \deprecated See CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T.
  \var  cs_rad_transfer_params_t::ifgrno
        \deprecated See CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX.
  \var  cs_rad_transfer_params_t::ifrefl
        \deprecated see CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX.
  \var  cs_rad_transfer_params_t::itpt1d
        See CS_BOUNDARY_RAD_WALL_GRAY_1D_T
  \var  cs_rad_transfer_params_t::ifinfe
        Modeling of an infinite extrusion for open boundaries.
  \var  cs_rad_transfer_params_t::time_control
        Determines at which time steps the variables are updated
        Also, in order to have proper initialization of the variables,
        the radiation module should always be called at
        the first time step of a calculation (restart or not).
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Model name */
const char *cs_rad_transfer_model_name[] = {
  N_("<none>"),
  N_("DOM (Discrete Ordinates Method)"),
  N_("P-1")};

/*! Quadrature name */
const char *cs_rad_transfer_quadrature_name[] = {
  "",
  "S4 (24 directions)",
  "S6 (48 directions)",
  "S8 (80 directions)",
  "T2 (32 directions)",
  "T4 (128 directions)",
  "TN (8*ndirec^2 directions)",
  "LC11",
  "DCT020_2468"};

cs_rad_transfer_params_t _rt_params = {
  .type = CS_RAD_TRANSFER_NONE,
  .nrphas = 1,
  .iimpar = 1,
  .verbosity = 0,
  .imodak = 0,
  .imoadf = 0,
  .iwrp1t = 0,
  .imfsck = 0,
  .xnp1mx = 10.,
  .idiver = -1,
  .i_quadrature = 1,
  .ndirec = 3,
  .ndirs = 0,
  .vect_s = NULL,
  .angsol = NULL,
  .restart = 0,
  .nwsgg = 1,
  .wq = NULL,
  .itpimp = CS_BOUNDARY_RAD_WALL_GRAY,
  .ipgrno = CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T,
  .iprefl = CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T,
  .ifgrno = CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX,
  .ifrefl = CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX,
  .itpt1d = CS_BOUNDARY_RAD_WALL_GRAY_1D_T,
  .ifinfe = 5,
  .atmo_model = CS_RAD_ATMO_3D_NONE,
  .atmo_dr_id = -1,
  .atmo_dr_o3_id = -1,
  .atmo_df_id = -1,
  .atmo_df_o3_id = -1,
  .atmo_ir_id = -1,
  .dispersion = false,
  .dispersion_coeff = 1.,
  .time_control = {
    .type = CS_TIME_CONTROL_TIME_STEP,
    .at_start = true,
    .at_first = true,
    .at_end = false,
    .start_nt = -1,
    .end_nt = -1,
    .interval_nt = 1,
    .control_func = NULL,
    .control_input = NULL,
    .current_state = false,
    .current_time_step = -1,
    .last_nt = -2,
    .last_t = -HUGE_VAL
  }
};

cs_rad_transfer_params_t *cs_glob_rad_transfer_params = &_rt_params;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_rad_transfer_get_pointers(int  **p_iirayo,
                             int  **p_nwsgg,
                             int  **p_nt_rad_prp,
                             int  **p_rad_atmo_model);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

void
cs_rad_transfer_get_pointers(int  **p_iirayo,
                             int  **p_nwsgg,
                             int  **p_nt_rad_prp,
                             int  **p_rad_atmo_model)
{
  *p_iirayo = (int *) &_rt_params.type;
  *p_nwsgg  = &_rt_params.nwsgg;
  *p_nt_rad_prp = &_rt_params.time_control.interval_nt;
  *p_rad_atmo_model = &_rt_params.atmo_model;
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize radiative transfer module.
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_finalize(void)
{
  BFT_FREE(_rt_params.vect_s);
  BFT_FREE(_rt_params.angsol);
  BFT_FREE(_rt_params.wq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Is time step for radiative transfer active?
 *
 * \return  true if active, false otherwise.
 */
/*----------------------------------------------------------------------------*/

bool
cs_rad_time_is_active(void) {
  cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;
  const cs_time_step_t *ts = cs_glob_time_step;
  return cs_time_control_is_active(&(rt_params->time_control), ts);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
