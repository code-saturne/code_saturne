/*============================================================================
 * General parameters management.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_parameters.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_parameters.c
        General parameters and options management.
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_space_disc_t

  \brief Space discretisation options descriptor.

  Members of this space discretisation structure are publicly accessible, to
  allow for concise syntax, as they are expected to be used in many places.

  \var  cs_space_disc_t::imvisf
        face viscosity field interpolation
        - 1: harmonic
        - 0: arithmetic (default)
  \var  cs_space_disc_t::imrgra
        type of gradient reconstruction
        - 0: iterative process
        - 1: standard least square method
        - 2: least square method with extended neighbourhood
        - 3: least square method with reduced extended neighbourhood
        - 4: iterative process initialized by the least square method
  \var  cs_space_disc_t::anomax
        non orthogonality angle of the faces, in radians.

        For larger angle values, cells with one node on the wall are kept in the
        extended support of the neighbouring cells.
  \var  cs_space_disc_t::iflxmw
        method to compute interior mass flux due to ALE mesh velocity
        - 1: based on cell center mesh velocity
        - 0: based on nodes displacement
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_piso_t

  \brief PISO options descriptor.

  Members of this PISO structure are publicly accessible, to allow for
  concise  syntax, as they are expected to be used in many places.

  \var  cs_piso_t::nterup
        number of interations on the pressure-velocity coupling on Navier-Stokes
  \var  cs_piso_t::epsup
        relative precision for the convergence test of the iterative process on
        pressure-velocity coupling
  \var  cs_piso_t::xnrmu
        norm  of the increment \f$ \vect{u}^{k+1} - \vect{u}^k \f$ of the
        iterative process on pressure-velocity coupling
  \var  cs_piso_t::xnrmu0
        norm of \f$ \vect{u}^0 \f$
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_var_cal_opt_t _var_cal_opt =
{
  0,     /* iwarni */
  1,     /* iconv  */
  1,     /* istat  */
  1,     /* idiff  */
  1,     /* idifft */
  1,     /* idften */
  0,     /* iswdyn */
  1,     /* ischcv */
  1,     /* isstpc */
  100,   /* nswrgr */
  1,     /* nswrsm */
  0,     /* imrgra */
  -1,    /* imligr */
  1,     /* ircflu */
  1.,    /* thetav */
  1.,    /* blencv */
  1.e-8, /* epsilo */
  1.e-7, /* epsrsm */
  1.e-5, /* epsrgr */
  1.5,   /* climgr */
  0.,    /* extrag */
  1.     /* relaxv */
};

/* Space discretisation options structure and associated pointer */

static cs_space_disc_t  _space_disc = {0, 0, -1e12*10.0, 1};

const cs_space_disc_t  *cs_glob_space_disc = &_space_disc;

/* PISO structure and associated pointer */

static cs_piso_t  _piso = {1, 1e-5, 0, 0};

const cs_piso_t  *cs_glob_piso = &_piso;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_space_disc_get_pointers(int     **imvisf,
                             int     **imrgra,
                             double  **anomax,
                             int     **iflxmw);

void
cs_f_piso_get_pointers(int     **nterup,
                       double  **epsup,
                       double  **xnrmu,
                       double  **xnrmu0);

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Log default values of the structure */

static void
_log_func_var_opt_cal(const void *t)
{
  const cs_var_cal_opt_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "iwarni", _t->iwarni);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "iconv ", _t->iconv );
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "istat ", _t->istat );
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "idiff ", _t->idiff );
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "idifft", _t->idifft);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "idften", _t->idften);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "iswdyn", _t->iswdyn);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "ischcv", _t->ischcv);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "isstpc", _t->isstpc);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "nswrgr", _t->nswrgr);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "nswrsm", _t->nswrsm);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "imrgra", _t->imrgra);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "imligr", _t->imligr);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-4d\n"),    "ircflu", _t->ircflu);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "thetav", _t->thetav);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "blencv", _t->blencv);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "epsilo", _t->epsilo);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "epsrsm", _t->epsrsm);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "epsrgr", _t->epsrgr);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "climgr", _t->climgr);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "extrag", _t->extrag);
  cs_log_printf(CS_LOG_SETUP, _("      %-19s  %-12.3g\n"), "relaxv", _t->relaxv);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Get pointers to members of the global space disc structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   imvisf  --> pointer to cs_glob_space_disc->imvisf
 *   imrgra  --> pointer to cs_glob_space_disc->imrgra
 *   anomax  --> pointer to cs_glob_space_disc->anomax
 *   iflxmw  --> pointer to cs_glob_space_disc->iflxmw
 *----------------------------------------------------------------------------*/

void
cs_f_space_disc_get_pointers(int     **imvisf,
                             int     **imrgra,
                             double  **anomax,
                             int     **iflxmw)
{
  *imvisf = &(_space_disc.imvisf);
  *imrgra = &(_space_disc.imrgra);
  *anomax = &(_space_disc.anomax);
  *iflxmw = &(_space_disc.iflxmw);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global piso structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   nterup  --> pointer to cs_glob_piso->nterup
 *   epsup   --> pointer to cs_glob_piso->epsup
 *   xnrmu   --> pointer to cs_glob_piso->xnrmu
 *   xnrmu0  --> pointer to cs_glob_piso->xnrmu0
 *----------------------------------------------------------------------------*/

void
cs_f_piso_get_pointers(int     **nterup,
                       double  **epsup,
                       double  **xnrmu,
                       double  **xnrmu0)
{
  *nterup = &(_piso.nterup);
  *epsup  = &(_piso.epsup);
  *xnrmu  = &(_piso.xnrmu);
  *xnrmu0 = &(_piso.xnrmu0);
}

/*! \endcond (end ignore by Doxygen) */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define general field keys.
 *
 * A recommened practice for different submodules would be to use
 * "cs_<module>_key_init() functions to define keys specific to those modules.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_define_field_keys(void)
{
  cs_field_define_key_int("inner_mass_flux_id", -1, 0);
  cs_field_define_key_int("boundary_mass_flux_id", -1, 0);

  cs_field_define_key_int("variable_id", -1, 0); /* inverse of ivarfl(ivar) */
  cs_field_define_key_int("post_id", -1, 0);     /* inverse of the ipp array */

  cs_field_define_key_int("scalar_diffusivity_id", -1, 0);
  cs_field_define_key_int("diffusivity_tensor", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("scalar_id", -1, 0); /* inverse of isca(iscal) */
  cs_field_define_key_int("drift_scalar_model", 0, 0);
  cs_field_define_key_int("scalar_class", -1, 0);
  cs_field_define_key_int("first_moment_id", -1, 0); // old iscavr(iscal)
  cs_field_define_key_double("min_scalar_clipping", -1.e12, 0);
  cs_field_define_key_double("max_scalar_clipping", 1.e12, 0);

  cs_field_define_key_int("property_id", -1, 0); /* inverse of iprpfl(iprop) */

  cs_field_define_key_struct("var_cal_opt",
                             &_var_cal_opt,
                             _log_func_var_opt_cal,
                             sizeof(cs_var_cal_opt_t),
                             CS_FIELD_VARIABLE);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read general restart info.
 *
 * This updates the previous time step info.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_read_restart_info(void)
{
  if (cs_restart_present()) {
    cs_restart_t *r
      = cs_restart_create("main", "restart", CS_RESTART_MODE_READ);
    cs_restart_read_time_step_info(r);
    r = cs_restart_destroy(r);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
