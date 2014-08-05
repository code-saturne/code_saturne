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

/* Definition of user variable */

typedef struct {

  char     *name;               /* Variable name */
  char     *ref_name;           /* Name of variable referred to */
  int       dim;                /* Variable dimension */
  bool      is_variance;        /* True if the variable is a variance */

} cs_user_variable_def_t;

/* Definition of user property */

typedef struct {

  char     *name;               /* Property name */
  int       dim;                /* Property dimension */
  int       location_id;        /* Propert location id */

} cs_user_property_def_t;

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

/* Definition of user variables and properties */

int                      _n_user_variables = 0;
int                      _n_user_properties = 0;
cs_user_variable_def_t  *_user_variable_defs = NULL;
cs_user_property_def_t  *_user_property_defs = NULL;

static cs_severe_acc_species_prop_t _severe_acc_species_prop =
{
  -1.,         /* molar mass              */
  -1.,          /* specific heat           */
  -1.,          /* volume diffusion        */
  -1.,   /* dynamic viscosity a     */
  -1.,   /* dynamic viscosity b     */
  -1.,        /* thermal condusctivity a */
  -1.,        /* thermal condusctivity b */
};

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

static void
_log_func_severe_acc_species_prop(const void *t)
{
  const cs_severe_acc_species_prop_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, _("      %-23s  %-12.3g\n"), "molar mass             ", _t->mol_mas);
  cs_log_printf(CS_LOG_SETUP, _("      %-23s  %-12.3g\n"), "specific heat          ", _t->cp);
  cs_log_printf(CS_LOG_SETUP, _("      %-23s  %-12.3g\n"), "volume diffusion       ", _t->vol_dif);
  cs_log_printf(CS_LOG_SETUP, _("      %-23s  %-12.3g\n"), "dynamic viscosity a    ", _t->mu_a);
  cs_log_printf(CS_LOG_SETUP, _("      %-23s  %-12.3g\n"), "dynamic viscosity b    ", _t->mu_b);
  cs_log_printf(CS_LOG_SETUP, _("      %-23s  %-12.3g\n"), "thermal condusctivity a", _t->lambda_a);
  cs_log_printf(CS_LOG_SETUP, _("      %-23s  %-12.3g\n"), "thermal condusctivity b", _t->lambda_b);
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

  cs_field_define_key_int("scalar_diffusivity_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_double("scalar_diffusivity_ref",
                             -1.e12*10., CS_FIELD_VARIABLE); /* old visls0(iscal) */
  cs_field_define_key_int("diffusivity_tensor", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("scalar_id", -1, 0); /* inverse of isca(iscal) */
  cs_field_define_key_int("drift_scalar_model", 0, 0);
  cs_field_define_key_int("scalar_class", 0, 0);
  cs_field_define_key_int("first_moment_id", -1, 0); // old iscavr(iscal)

  cs_field_define_key_int("source_term_prev_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_double("min_scalar_clipping", -1.e12, 0);
  cs_field_define_key_double("max_scalar_clipping", 1.e12, 0);

  cs_field_define_key_int("property_id", -1, 0); /* inverse of iprpfl(iprop) */

  cs_field_define_key_struct("var_cal_opt",
                             &_var_cal_opt,
                             _log_func_var_opt_cal,
                             sizeof(cs_var_cal_opt_t),
                             CS_FIELD_VARIABLE);

  /* Structure containing physical properties relative to
     severe accident scalars */
  cs_field_define_key_struct("severe_acc_species_prop",
                             &_severe_acc_species_prop,
                             _log_func_severe_acc_species_prop,
                             sizeof(cs_severe_acc_species_prop_t),
                             0);

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
/*!
 * \brief Define a user variable.
 *
 * Solved variables are always defined on cells.
 *
 * \param[in]  name  name of variable and associated field
 * \param[in]  dim   variable dimension
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_add_variable(const char  *name,
                           int          dim)
{
  BFT_REALLOC(_user_variable_defs,
              _n_user_variables + 1,
              cs_user_variable_def_t);

  BFT_MALLOC((_user_variable_defs + _n_user_variables)->name,
             strlen(name) + 1,
             char);
  strcpy((_user_variable_defs + _n_user_variables)->name, name);

  (_user_variable_defs + _n_user_variables)->dim = dim;
  (_user_variable_defs + _n_user_variables)->is_variance = 0;

  if (dim != 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Only user variables of dimension 1 are currently handled,\n"
                "but %s is defined with dimension %d."),
              name, dim);

  _n_user_variables++;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a user variable which is a variance of another variable.
 *
 * Only variances of thermal or user-defined variables are currently handled.
 *
 * \param[in]  name           name of variance and associated field
 * \param[in]  variable_name  name of associated variable
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_add_variable_variance(const char  *name,
                                    const char  *variable_name)
{
  BFT_REALLOC(_user_variable_defs,
              _n_user_variables + 1,
              cs_user_variable_def_t);
  BFT_MALLOC((_user_variable_defs + _n_user_variables)->name,
             strlen(name) + 1,
             char);
  BFT_MALLOC((_user_variable_defs + _n_user_variables)->ref_name,
             strlen(variable_name) + 1,
             char);

  strcpy((_user_variable_defs + _n_user_variables)->name, name);
  strcpy((_user_variable_defs + _n_user_variables)->ref_name, variable_name);
  (_user_variable_defs + _n_user_variables)->dim = -1;
  (_user_variable_defs + _n_user_variables)->is_variance = true;

  _n_user_variables++;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a user property.
 *
 * \param[in]  name         name of property and associated field
 * \param[in]  dim          property dimension
 * \param[in]  location_id  id of associated mesh location
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_add_property(const char  *name,
                           int          dim,
                           int          location_id)
{
  BFT_REALLOC(_user_property_defs,
              _n_user_properties + 1,
              cs_user_property_def_t);
  BFT_MALLOC((_user_property_defs + _n_user_properties)->name,
             strlen(name) + 1,
             char);

  strcpy((_user_property_defs + _n_user_properties)->name, name);
  (_user_property_defs + _n_user_properties)->dim = dim;
  (_user_property_defs + _n_user_properties)->location_id = location_id;

  _n_user_properties++;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of defined user variables not added yet.
 *
 * This number is reset to 0 when \ref cs_parameters_create_added_variables
 * is called.
 *
 * \return number of defined user variables
 */
/*----------------------------------------------------------------------------*/

int
cs_parameters_n_added_variables(void)
{
  return _n_user_variables;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of defined user properties not added yet.
 *
 * \return number of defined user properties
 */
/*----------------------------------------------------------------------------*/

int
cs_parameters_n_added_properties(void)
{
  return _n_user_properties;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create previously added user variables.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_create_added_variables(void)
{
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_USER;

  for (int i = 0; i < _n_user_variables; i++) {

    cs_field_t *f;

    const char *name = (_user_variable_defs + i)->name;

    int cmp_id = cs_field_id_by_name(name);

    if (cmp_id > -1)
      bft_error(__FILE__, __LINE__, 0,
                _("Error defining user variable \"%s\";\n"
                  "this name is already reserved for field with id %d."),
                name, cmp_id);

    /* Case where we define a variance */

    if ((_user_variable_defs + i)->is_variance) {

      const char *ref_name = (_user_variable_defs + i)->ref_name;
      const cs_field_t *f_ref = cs_field_by_name_try(ref_name);

      if (f_ref == NULL)
        bft_error(__FILE__, __LINE__, 0,
                  _("Error defining user variance \"%s\";\n"
                    "which refers to yet undefined variable \"%s\"."),
                  name, ref_name);

      f = cs_field_create(name,
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          f_ref->dim,
                          true,
                          true);
      int k_var = cs_field_key_id("first_moment_id");
      cs_field_set_key_int(f, k_var, f_ref->id);
      cs_field_lock_key(f, k_var);
      BFT_FREE((_user_variable_defs + i)->ref_name);

    }

    /* General case */

    else {

      f = cs_field_create(name,
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          (_user_variable_defs + i)->dim,
                          true,
                          true);

    }

    BFT_FREE((_user_variable_defs + i)->name);

    cs_field_set_key_int(f, cs_field_key_id("log"), 1);
    cs_field_set_key_int(f, cs_field_key_id("post_vis"), 1);

  }

  BFT_FREE(_user_variable_defs);
  _n_user_variables = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create previously added user properties.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_create_added_properties(void)
{
  /* Define variable diffusivities for the temperature or
     user-defined variables */

  /* Define regular user properties */

  for (int i = 0; i < _n_user_properties; i++) {

    const char *name = (_user_property_defs + i)->name;

    int cmp_id = cs_field_id_by_name(name);

    if (cmp_id > -1)
      bft_error(__FILE__, __LINE__, 0,
                _("Error defining user property \"%s\";\n"
                  "this name is already reserved for field with id %d."),
                name, cmp_id);

    (void) cs_field_create(name,
                           CS_FIELD_PROPERTY | CS_FIELD_USER,
                           (_user_property_defs + i)->location_id,
                           (_user_property_defs + i)->dim,
                           false,
                           false);

    BFT_FREE((_user_property_defs + i)->name);

  }

  BFT_FREE(_user_property_defs);
  _n_user_properties = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
