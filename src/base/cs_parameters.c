/*============================================================================
 * General parameters management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include <stdarg.h>
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

/*
  \enum parameter_error_behavior_t

  \brief File acces modes

  \var CS_WARN
       Warn only
  \var CS_ABORT_DELAYED
       Abort when \ref cs_parameters_error_barrier is called.
  \var CS_FILE_MODE_APPEND
       Abort immediately

  \struct cs_space_disc_t

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
        - 2: least square method with extended neighborhood
        - 3: least square method with reduced extended neighborhood
        - 4: iterative process initialized by the least square method
  \var  cs_space_disc_t::anomax
        non orthogonality angle of the faces, in radians.

        For larger angle values, cells with one node on the wall are kept in the
        extended support of the neighboring cells.
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

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Counter for parameter checking errors */

static int  _param_check_errors = 0;

/* Default variable compute options */

static cs_var_cal_opt_t _var_cal_opt =
{
  .iwarni = 0,
  .iconv  = 1,
  .istat  = 1,
  .idiff  = 1,
  .idifft = 1,
  .idften = 1,
  .iswdyn = 0,
  .ischcv = 1,
  .ibdtso = 1,
  .isstpc = 1,
  .nswrgr = 100,
  .nswrsm = 1,
  .imrgra = 0,
  .imligr = -1,
  .ircflu = 1,
  .iwgrec = 0,
  .thetav = 1.,
  .blencv = 1.,
  .epsilo = 1.e-8,
  .epsrsm = 1.e-7,
  .epsrgr = 1.e-5,
  .climgr = 1.5,
  .extrag = 0.,
  .relaxv = 1.
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

static cs_solving_info_t _solving_info =
{
  0,     /* n_it: number of iterations for the linear solver */
  0.,    /* rhs_norm: right hand side norm                   */
  0.,    /* res_norm: normed residual                        */
  0.,    /* derive: norm of the time derivative              */
  0.,    /* l2residual: L2 time residual                     */
};

static cs_gas_mix_species_prop_t _gas_mix_species_prop =
{
  -1.,   /* molar mass                              */
  -1.,   /* specific heat                           */
  -1.,   /* volume diffusion                        */
  -1.,   /* dynamic viscosity a                     */
  -1.,   /* dynamic viscosity b                     */
  -1.,   /* thermal conductivity a                  */
  -1.,   /* thermal conductivity b                  */
  -1.,   /* reference viscosity (Sutherland)        */
  -1.,   /* reference conductivity (Sutherland)     */
  -1.,   /* reference temperature for viscosity     */
  -1.,   /* reference temperature for conductivity  */
  -1.,   /* Sutherland temperature for viscosity    */
  -1.,   /* Sutherland temperature for conductivity */
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

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Log default values of the structure */

static void
_log_func_var_opt_cal(const void *t)
{
  const char fmt_i[] = N_("      %-19s  %-4d\n");
  const char fmt_r[] = N_("      %-19s  %-12.3g\n");
  const cs_var_cal_opt_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "iwarni", _t->iwarni);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "iconv ", _t->iconv );
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "istat ", _t->istat );
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "idiff ", _t->idiff );
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "idifft", _t->idifft);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "idften", _t->idften);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "iswdyn", _t->iswdyn);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "ischcv", _t->ischcv);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "ibdtso", _t->ibdtso);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "isstpc", _t->isstpc);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "nswrgr", _t->nswrgr);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "nswrsm", _t->nswrsm);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "imrgra", _t->imrgra);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "imligr", _t->imligr);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "ircflu", _t->ircflu);
  cs_log_printf(CS_LOG_SETUP, _(fmt_i), "iwgrec", _t->iwgrec);
  cs_log_printf(CS_LOG_SETUP, _(fmt_r), "thetav", _t->thetav);
  cs_log_printf(CS_LOG_SETUP, _(fmt_r), "blencv", _t->blencv);
  cs_log_printf(CS_LOG_SETUP, _(fmt_r), "epsilo", _t->epsilo);
  cs_log_printf(CS_LOG_SETUP, _(fmt_r), "epsrsm", _t->epsrsm);
  cs_log_printf(CS_LOG_SETUP, _(fmt_r), "epsrgr", _t->epsrgr);
  cs_log_printf(CS_LOG_SETUP, _(fmt_r), "climgr", _t->climgr);
  cs_log_printf(CS_LOG_SETUP, _(fmt_r), "extrag", _t->extrag);
  cs_log_printf(CS_LOG_SETUP, _(fmt_r), "relaxv", _t->relaxv);
}

static void
_log_func_gas_mix_species_prop(const void *t)
{
  const char fmt[] = N_("      %-23s  %-12.3g\n");
  const cs_gas_mix_species_prop_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, _(fmt), "molar mass            ", _t->mol_mas);
  cs_log_printf(CS_LOG_SETUP, _(fmt), "specific heat         ", _t->cp);
  cs_log_printf(CS_LOG_SETUP, _(fmt), "volume diffusion      ", _t->vol_dif);
  cs_log_printf(CS_LOG_SETUP, _(fmt), "dynamic viscosity a   ", _t->mu_a);
  cs_log_printf(CS_LOG_SETUP, _(fmt), "dynamic viscosity b   ", _t->mu_b);
  cs_log_printf(CS_LOG_SETUP, _(fmt), "thermal conductivity a", _t->lambda_a);
  cs_log_printf(CS_LOG_SETUP, _(fmt), "thermal conductivity b", _t->lambda_b);
  cs_log_printf(CS_LOG_SETUP, _(fmt),
                "reference thermal viscosity (Sutherland)",
                _t->muref);
  cs_log_printf(CS_LOG_SETUP, _(fmt),
                "reference thermal conductivity (Sutherland)",
                _t->lamref);
  cs_log_printf(CS_LOG_SETUP, _(fmt),
                "reference temperature (Sutherland for viscosity)",
                _t->trefmu);
  cs_log_printf(CS_LOG_SETUP, _(fmt),
                "reference temperature (Sutherland conductivity)",
                _t->treflam);
  cs_log_printf(CS_LOG_SETUP, _(fmt),
                "Sutherland temperature for viscosity",
                _t->smu);
  cs_log_printf(CS_LOG_SETUP, _(fmt),
                "Sutherland tempertaure for conductivity",
                _t->slam);

}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Provide acces to cs_glob_piso
 *
 * needed to initialize structure with GUI
 *
 * \return  piso information structure
 */
/*----------------------------------------------------------------------------*/

cs_piso_t *
cs_get_glob_piso(void)
{
  return &_piso;
}

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
  cs_field_define_key_int("scalar_id", -1, 0);   /* inverse of isca(iscal) */
  cs_field_define_key_int("post_id", -1, 0);     /* inverse of the ipp array */

  cs_field_define_key_int("scalar_diffusivity_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_double("scalar_diffusivity_ref",
                             -1.e12*10., CS_FIELD_VARIABLE); /* visls0(iscal) */

  cs_field_define_key_int("turbulent_flux_model", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("turbulent_flux_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("gradient_weighting_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("diffusivity_tensor", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("drift_scalar_model", 0, 0);

  cs_field_define_key_int("scalar_class", 0, 0);
  cs_field_define_key_int("first_moment_id", -1, 0); /* iscavr(iscal) */

  cs_field_define_key_int("source_term_prev_id", -1, CS_FIELD_VARIABLE);
  /* TODO merge with previous key word */
  cs_field_define_key_int("source_term_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("slope_test_upwind_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("boundary_value_id", -1, 0);

  cs_field_define_key_int("convection_limiter_id", -1, 0);

  cs_field_define_key_double("min_scalar_clipping", -1.e12, 0);
  cs_field_define_key_double("max_scalar_clipping", 1.e12, 0);

  /* Bounds of a given scalar which won't be used in clipping */

  cs_field_define_key_double("max_scalar", 1., 0);
  cs_field_define_key_double("min_scalar", 0., 0);

 /* Integer corresponding to the type of Roe-Sweby Limiter:
  * 1->minmod
  * 2->van-leer
  * 3->van-albada
  * 4->superbee */

  cs_field_define_key_int("limiter_choice", -1, CS_FIELD_VARIABLE);

  /* Structure containing the calculation options of the field variables */
  cs_field_define_key_struct("var_cal_opt",
                             &_var_cal_opt,
                             _log_func_var_opt_cal,
                             sizeof(cs_var_cal_opt_t),
                             CS_FIELD_VARIABLE);

  /* Structure containing the solving info of the field variables
     (used for listing, not setup, so set NULL setup logging function) */
  cs_field_define_key_struct("solving_info",
                             &_solving_info,
                             NULL,
                             sizeof(cs_solving_info_t),
                             CS_FIELD_VARIABLE);
  cs_field_key_disable_setup_log(cs_field_key_id("solving_info"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define field key for condensation.
 *
 * Note: this should be moved in the future to a condensation-specific file.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_define_field_key_gas_mix(void)
{
  /* Structure containing physical properties relative to
     species scalars used by the gas mixture modelling */
  cs_field_define_key_struct("gas_mix_species_prop",
                             &_gas_mix_species_prop,
                             _log_func_gas_mix_species_prop,
                             sizeof(cs_gas_mix_species_prop_t),
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
    cs_restart_destroy(&r);
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
  (_user_variable_defs + _n_user_variables)->is_variance = false;

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

    cs_field_t *fld =
      cs_field_create(name,
                      CS_FIELD_PROPERTY | CS_FIELD_USER,
                      (_user_property_defs + i)->location_id,
                      (_user_property_defs + i)->dim,
                      false);

    cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
    cs_field_set_key_int(fld, cs_field_key_id("post_vis"), 1);

    BFT_FREE((_user_property_defs + i)->name);

  }

  BFT_FREE(_user_property_defs);
  _n_user_properties = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a boundary values field for a variable field.
 *
 * \param[in, out]  f  pointer to field structure
 *
 * \return  pointer to boundary values field, or NULL if not applicable
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_parameters_add_boundary_values(cs_field_t  *f)
{
  cs_field_t *bf = NULL;

  /* Check we are on cells and don't already have such value */

  if (f->location_id != CS_MESH_LOCATION_CELLS)
    return bf;

  int kbf = cs_field_key_id_try("boundary_value_id");

  int bf_id = cs_field_get_key_int(f, kbf);
  if (bf_id > -1) {
    bf = cs_field_by_id(bf_id);
    return bf;
  }

  /* Currently only managed for scalars or temperature property */

  int ks = cs_field_key_id_try("scalar_id");
  if (ks < 0)
    return bf;

  int scalar_id = (f->type & CS_FIELD_VARIABLE) ?
    cs_field_get_key_int(f, ks) : -1;

  if (scalar_id < 0&& strcmp(f->name, "temperature") != 0)
    return bf;

  /* Build new field */

  char *b_name;
  size_t l = strlen("boundary_") + strlen(f->name) + 1;
  BFT_MALLOC(b_name, l, char);
  snprintf(b_name, l, "boundary_%s", f->name);

  /* Field may already have been defined */

  bf = cs_field_by_name_try(b_name);

  if (bf == NULL) {

    int type_flag =   (f->type & (CS_FIELD_INTENSIVE | CS_FIELD_EXTENSIVE))
                    | CS_FIELD_POSTPROCESS;

    bf = cs_field_create(b_name,
                         type_flag,
                         CS_MESH_LOCATION_BOUNDARY_FACES,
                         f->dim,
                         false);

    /* Set same label as parent */

    cs_field_set_key_str(bf,
                         cs_field_key_id("label"),
                         cs_field_get_label(f));

    /* Set same postprocessing and logging defaults as parent */

    int k_log = cs_field_key_id("log");
    cs_field_set_key_int(bf,
                         k_log,
                         cs_field_get_key_int(f, k_log));

    int k_vis = cs_field_key_id("post_vis");
    int f_vis = cs_field_get_key_int(f, k_vis);
    f_vis = CS_MAX(f_vis, 1);
    cs_field_set_key_int(bf, k_vis, f_vis);

  }
  else {

    if (   f->dim != bf->dim
        || bf->location_id != CS_MESH_LOCATION_BOUNDARY_FACES)
      bft_error(__FILE__, __LINE__, 0,
                _("Error defining variable boundary field:\n"
                  "  parent name:   \"%s\"\n"
                  "  name:          \"%s\"\n"
                  "  dimension:     %d\n\n"
                  "An incompatible field with matching name already exists:\n"
                  "  id:          %d\n"
                  "  location_id: %d\n"
                  "  dimension:   %d"),
                f->name, bf->name, f->dim,
                bf->id, bf->location_id, bf->dim);

  }

  BFT_FREE(b_name);

  cs_field_set_key_int(f, kbf, bf->id);
  cs_field_lock_key(f, kbf);

  return bf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a boundary values field for temperature, if applicable.
 *
 * \return  pointer to boundary values field, or NULL if not applicable
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_parameters_add_boundary_temperature(void)
{
  cs_field_t *bf = NULL;

  /* Check if we already have a temperature variable field
    (temperature or enthalpy */

  cs_field_t *f = cs_field_by_name_try("temperature");

  if (f != NULL)
    bf = cs_parameters_add_boundary_values(f);

  else {

    f = cs_field_by_name_try("enthalpy");

    if (f != NULL) {
      if (   f->location_id != CS_MESH_LOCATION_CELLS
          || (f->type & CS_FIELD_VARIABLE) == 0)
        f = NULL;
    }

    /* If we have a compatible cell enthalpy field,
       use if to define default output options */

    if (f != NULL) {

      char b_name[] = "boundary_temperature";

      bf = cs_field_by_name_try(b_name);

      if (bf == NULL) {

        int type_flag =   (f->type & (CS_FIELD_INTENSIVE | CS_FIELD_EXTENSIVE))
                           | CS_FIELD_POSTPROCESS;

        bf = cs_field_create(b_name,
                             type_flag,
                             CS_MESH_LOCATION_BOUNDARY_FACES,
                             f->dim,
                             false);

        /* Set same postprocessing and logging defaults as enthalpy */

        int k_log = cs_field_key_id("log");
        cs_field_set_key_int(bf,
                             k_log,
                             cs_field_get_key_int(f, k_log));

        int k_vis = cs_field_key_id("post_vis");
        int f_vis = cs_field_get_key_int(f, k_vis);
        f_vis = CS_MAX(f_vis, 1);
        cs_field_set_key_int(bf, k_vis, f_vis);

      }
      else {

        if (   1 != bf->dim
            || bf->location_id != CS_MESH_LOCATION_BOUNDARY_FACES)
          bft_error(__FILE__, __LINE__, 0,
                    _("Error defining variable \"boundary_temperature\" field:\n"
                      "An incompatible field with matching name already exists:\n"
                      "  id:          %d\n"
                      "  location_id: %d\n"
                      "  dimension:   %d"),
                    bf->id, bf->location_id, bf->dim);

      }

    }

  }

  return bf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a local variable calculation options structure,
 *        with default options.
 *
 * \return  variable calculations options structure
 */
/*----------------------------------------------------------------------------*/

cs_var_cal_opt_t
cs_parameters_var_cal_opt_default(void)
{
  return _var_cal_opt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print general parameters error or warning info.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param [in] format        format string, as printf() and family.
 * \param [in] ...           variable arguments based on format string.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_error(cs_parameter_error_behavior_t   err_behavior,
                    const char                     *section_desc,
                    const char                     *format,
                    ...)
{
  cs_parameters_error_header(err_behavior, section_desc);

  int log_id = CS_LOG_DEFAULT;

  va_list  arg_ptr;
  va_start(arg_ptr, format);

  cs_log_vprintf(log_id, format, arg_ptr);

  va_end(arg_ptr);

  cs_parameters_error_footer(err_behavior);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print header for a given parameters error message type.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_error_header(cs_parameter_error_behavior_t   err_behavior,
                           const char                     *section_desc)
{
  const int err_type_id = (err_behavior <= CS_WARNING) ? 0 : 1;
  const char *error_type[] = {N_("Warning"),
                              N_("Error")};

  int log_id = CS_LOG_DEFAULT;

  if (section_desc != NULL)
    cs_log_printf(log_id, "%s %s\n", _(error_type[err_type_id]), section_desc);
  else
    cs_log_printf(log_id, "%s\n", _(error_type[err_type_id]));
  size_t l = cs_log_strlen(_(error_type[err_type_id]));
  char underline[81];

  for (size_t i = 0; i < 80 && i < l; i++)
    underline[i] = '-';
  underline[CS_MIN(l,80)] = '\0';
  cs_log_printf(log_id, "%s\n", underline);

  if (err_behavior > CS_WARNING)
    _param_check_errors++;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print footer for a given parameters error message type.
 *
 * \param[in]  err_behavior  warn or abort ?
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_error_footer(cs_parameter_error_behavior_t   err_behavior)
{
  if (err_behavior == CS_ABORT_IMMEDIATE)
    bft_error
      (__FILE__, __LINE__, 0,
       _("\nCheck your data and parameters (GUI and user subroutines)."));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given integer keyword has values in a specified range.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  range_l       range lower bound (included)
 * \param[in]  range_u       range upper bound (excluded)
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_in_range_int(cs_parameter_error_behavior_t   err_behavior,
                              const char                     *section_desc,
                              const char                     *param_name,
                              int                             param_value,
                              int                             range_l,
                              int                             range_u)
{
  if (param_value < range_l || param_value >= range_u) {

    cs_parameters_error_header(err_behavior, section_desc);

    int log_id = CS_LOG_DEFAULT;

    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be in range [%d, %d].\n"),
                  param_name, param_value, range_l, range_u-1);

    cs_parameters_error_footer(err_behavior);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check that a given integer keyword has values in a specified range.
 *
 * \param[in]  err_behavior  warn or abort ?
 * \param[in]  section_desc  optional description of code section
 *                           containing this parameter, or NULL
 * \param[in]  param_name    name of parameter whose value we are checking
 * \param[in]  param_value   parameter's current_value
 * \param[in]  enum_size     size of possible enumeration
 * \param[in]  enum_values   optional list of enumerated values, or NULL
 *                           (in which case {0, ... enum_sizes-1} assumed
 * \param[in]  enum_names    optional list of value names, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_is_in_list_int(cs_parameter_error_behavior_t   err_behavior,
                             const char                     *section_desc,
                             const char                     *param_name,
                             int                             param_value,
                             int                             enum_size,
                             const int                      *enum_values,
                             const char                     *enum_names[])
{
  /* Check if we are in the defined range */

  if (enum_values != NULL) {
    for (int i = 0; i < enum_size; i++) {
      if (param_value == enum_values[i])
        return;
    }
  }
  else if (param_value >= 0 && param_value < enum_size)
    return;

  /* If we are not, report error */

  cs_parameters_error_header(err_behavior, section_desc);

  int log_id = CS_LOG_DEFAULT;

  if (enum_names != NULL) {
    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be one of:\n"),
                  param_name, param_value);
    for (int i = 0; i < enum_size; i++)
      cs_log_printf(log_id, "  %s\n", enum_names[i]);
  }
  else if (enum_values != NULL) {
    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be one of:\n"),
                  param_name, param_value);
    for (int i = 0; i < enum_size; i++)
      cs_log_printf(log_id, "  %d\n", i);
  }
  else {
    cs_log_printf(log_id,
                  _("Parameter: %s = %d\n"
                    "while its value must be in range [%d, %d].\n"),
                  param_name, param_value, 0, enum_size-1);
  }

  cs_parameters_error_footer(err_behavior);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Abort if the the parameter errors count is nonzero.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_error_barrier(void)
{
  cs_lnum_t n_errors = _param_check_errors;
  cs_parall_counter_max(&n_errors, 1);

  if (n_errors > 0)
    bft_error
      (__FILE__, __LINE__, 0,
       _("%d parameter error(s) reported.\n"
         "\n"
         "Read error messages above for details, then\n"
         "check your data and parameters (GUI and user subroutines)."),
       n_errors);

  _param_check_errors = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
