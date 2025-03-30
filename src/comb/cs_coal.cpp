/*============================================================================
 * Coal combustion model.
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

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "cdo/cs_equation_param.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_parameters.h"
#include "base/cs_parameters_check.h"
#include "base/cs_physical_constants.h"
#include "base/cs_physical_properties.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_post.h"
#include "base/cs_thermal_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "comb/cs_coal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal.cpp

  \brief Coal combustion model.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Coal combustion model parameters structure */

cs_coal_model_t  *cs_glob_coal_model = NULL;

/*!>  molar volume under normal pressure and temperature conditions
   (1 atmosphere, 0 degres C) in m-3 */

/*! reference temperature for molar volume */
const double  cs_coal_trefth = 25. + 273.15;

/*! reference pressure for molar volume */
const double  cs_coal_prefth = 1.01325e5;

/*! molar volume under normal pressure and temperature conditions
  (1 atmosphere, 0 \f$\text{\degresC}\f$) in \f$m^{-3}\f$ */
const double  cs_coal_volmol = 22.41e-3;

/* ids for atom types in wmolat */
const int  cs_coal_atom_id_c = 0;  /*!< id for C in wmolat */
const int  cs_coal_atom_id_h = 1;  /*!< id for H in wmolat */
const int  cs_coal_atom_id_o = 2;  /*!< id for O in wmolat */
const int  cs_coal_atom_id_n = 3;  /*!< id for N in wmolat */
const int  cs_coal_atom_id_s = 4;  /*!< id for S in wmolat */

/* precision for tests */
const double cs_coal_epsilon = 1.e-8;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/* Additional prototypes for Fortran mappings */

int
cs_add_model_field_indexes(int  f_id);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize coal model.
 *
 * \pram[in, out]  cm  pointer to coal model pointer to destroy.
 */
/*----------------------------------------------------------------------------*/

static void
_coal_model_finalize(void)
{
  CS_FREE(cs_glob_coal_model);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a coal model field.
 *
 * \param[in]  base_name   field base name
 * \param[in]  base_label  field base label
 * \param[in]  cc_id       class or coal id, or -1 if ignored
 * \param[in]  class_id    class or coal id
 * \param[in]  drift_flag  drift flag
 *
 * \return  pointer to field;
 */
/*----------------------------------------------------------------------------*/

static cs_field_t *
_add_coal_variable(const char  *base_name,
                   const char  *base_label,
                   int          cc_id,
                   int          class_id,
                   int          drift_flag)
{
  char name[64], label[64];
  if (cc_id > -1) {
    snprintf(name, 64, "%s_%02d", base_name, cc_id+1); name[63] = '\0';
    snprintf(label, 64, "%s_%02d", base_label, cc_id+1); label[63] = '\0';
  }
  else {
    strncpy(name, base_name, 64); name[63] = '\0';
    strncpy(label, base_label, 64); label[63] = '\0';
  }

  int f_id = cs_variable_field_create(name, label, CS_MESH_LOCATION_CELLS, 1);
  cs_field_t *f = cs_field_by_id(f_id);

  cs_add_model_field_indexes(f->id);

  // Set the index of the scalar class in the field structure
  // TODO: we could probably leave the id at the default value
  // (0, could be -1) when -1 is passed.

  const int keyccl = cs_field_key_id("scalar_class");
  int class_num = (class_id > -1) ? class_id + 1 : -1;
  cs_field_set_key_int(f, keyccl, class_num);

  if (drift_flag != 0) {
    const int keydri = cs_field_key_id("drift_scalar_model");
    cs_field_set_key_int(f, keydri, drift_flag);
  }

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a property field at cells.
 *
 * \param[in]  base_name   field base name
 * \param[in]  base_label  field base label
 * \param[in]  cc_id       class or coal id, or -1 if ignored
 * \param[in]  dim         field dimension
 *
 * Note that for consistency with older behavior, no "_" is systematically
 * added between parts of a composite name, so the caller can choose whether
 * to use the character or not.
 *
 * \return  pointer to field
 */
/*----------------------------------------------------------------------------*/

static cs_field_t *
_add_coal_property(const char  *base_name,
                   const char  *base_label,
                   int          cc_id,
                   int          dim)
{
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");
  const int keylbl = cs_field_key_id("label");

  char name[64], label[64];
  if (cc_id > -1) {
    snprintf(name, 64, "%s%02d", base_name, cc_id+1); name[63] = '\0';
    snprintf(label, 64, "%s%02d", base_label, cc_id+1); label[63] = '\0';
  }
  else {
    strncpy(name, base_name, 64); name[63] = '\0';
    strncpy(label, base_label, 64); label[63] = '\0';
  }

  if (cs_field_by_name_try(name) != NULL)
    cs_parameters_error(CS_ABORT_IMMEDIATE,
                        _("initial data setup"),
                        _("Field %s has already been assigned.\n"),
                        name);

  cs_physical_property_define_from_field(name,
                                         CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                                         CS_MESH_LOCATION_CELLS,
                                         dim,
                                         false);  // has previous

  int f_id = cs_physical_property_field_id_by_name(name);
  cs_field_t *f = cs_field_by_id(f_id);

  cs_field_set_key_int(f, keyvis, 0);
  cs_field_set_key_int(f, keylog, 1);
  cs_field_set_key_str(f, keylbl, label);

  return f;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return coal combustion model type.
 *
 * \return type of active coal combustion model
 *         (CS_COMBUSTION_COAL_NONE if model not active)
 */
/*----------------------------------------------------------------------------*/

cs_coal_model_type_t
cs_coal_model_get_type(void)
{
  if (cs_glob_coal_model == NULL)
   return CS_COMBUSTION_COAL_NONE;
  else
    return cs_glob_coal_model->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate coal combustion model.
 *
 * \return  pointer to coal combustion model structure.
 *
 * \param[in]  type  coal combustion model type
 */
/*----------------------------------------------------------------------------*/

cs_coal_model_t *
cs_coal_model_set_model(cs_coal_model_type_t  type)
{
  cs_glob_physical_model_flag[CS_COMBUSTION_COAL] = type;

  if (type == CS_COMBUSTION_COAL_NONE) {
    CS_FREE(cs_glob_coal_model);
    return NULL;
  }
  else if (cs_glob_coal_model != NULL) {
    cs_glob_coal_model->type = type;
    return cs_glob_coal_model;
  }

  /* Create and initialize model structure */

  cs_coal_model_t *cm = NULL;

  CS_MALLOC(cm, 1, cs_coal_model_t);
  memset(cm, 0, sizeof(cs_coal_model_t));

  cs_glob_coal_model = cm;

  /* Members whose initial value is not 0 */

  cm->type = type;

  cm->ieqnox = 1;
  cm->ieqco2 = 1;

  cm->ico = -1;
  cm->io2 = -1;
  cm->in2 = -1;
  cm->ico2 = -1;
  cm->ih2o = -1;

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_CLASSES; i++)
    cm->ighh2o[i] = -1;

  cm->ihgas = -1;
  cm->iyco2 = -1;
  cm->iyhcn = -1;
  cm->iynh3 = -1;
  cm->iyno  = -1;

  cm->ihox  = -1;

  for (int i = 0; i < CS_COMBUSTION_MAX_COALS; i++) {
    cm->if1m[i] = -1;
    cm->if2m[i] = -1;
  }

  cm->if4m = -1;
  cm->if5m = -1;
  cm->if6m = -1;
  cm->if7m = -1;
  cm->if8m = -1;
  cm->if9m = -1;
  cm->ifvp2m = -1;

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_CLASSES; i++) {
    cm->ixck[i] = -1;
    cm->ixch[i] = -1;
    cm->inp[i] = -1;
    cm->ih2[i] = -1;
    cm->ixwt[i] = -1;
  }

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS; i++)
    cm->iym1[i] = -1;

  cm->irom1 = -1;
  cm->immel = -1;

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_CLASSES; i++) {
    cm->itemp2[i] = -1;
    cm->irom2[i] = -1;
    cm->idiam2[i] = -1;
    cm->ix2[i] = -1;
    cm->igmdch[i] = -1;
    cm->igmhet[i] = -1;
    cm->igmtr[i] = -1;
    cm->ighco2[i] = -1;
    cm->igmdv1[i] = -1;
    cm->igmdv2[i] = -1;
    cm->igmsec[i] = -1;
  }

  cm->ibcarbone = -1;
  cm->iboxygen = -1;
  cm->ibhydrogen = -1;

  cm->ighcn1 = -1;
  cm->ighcn2 = -1;
  cm->ignoth = -1;
  cm->ignh31 = -1;
  cm->ignh32 = -1;
  cm->ifhcnd = -1;
  cm->ifhcnc = -1;
  cm->ifnh3d = -1;
  cm->ifnh3c = -1;
  cm->ifnohc = -1;
  cm->ifnonh = -1;
  cm->ifnoch = -1;
  cm->ifnoth = -1;
  cm->icnohc = -1;
  cm->icnonh = -1;
  cm->ifhcnr = -1;
  cm->icnorb = -1;

  cm->srrom = 0.95;

  /* Set finalization callback */

  cs_base_at_finalize(_coal_model_finalize);

  return cm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Specific setup operations for pulverized coal combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_setup(void)
{
  cs_coal_model_t *cm = cs_glob_coal_model;

  /* Transported variables
     --------------------- */

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();

  const int keysca = cs_field_key_id("scalar_id");
  const int kscacp  = cs_field_key_id("is_temperature");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");

  const int thm_f_id = cs_thermal_model_field()->id;

  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE))
      continue;
    if (f->type & CS_FIELD_CDO || f->type & CS_FIELD_USER)
      continue;
    if (cs_field_get_key_int(f, keysca) <= 0)
      continue;

    /* Field is a model variable (transported) */

    cs_field_set_key_int(f, kscacp, 0);  /* Should already be set or default */

    /* Solve on enthalpy with constant CP */

    if (f_id != thm_f_id && cs_field_get_key_int(f, kscavr) < 0) {

      /* For combustion, we consider that the turbulent viscosity dominates.
         We avoid computing laminar flames with  Le =/= 1. */

      cs_field_set_key_double(f, kvisl0, fp->viscl0);

    }

    /* Turbulent Schmidt or Prandtl */

    const cs_real_t turb_schmidt = 0.7;
    cs_field_set_key_double(f, ksigmas, turb_schmidt);

    /* Equation parameter defaults */

    cs_equation_param_t *eqp = cs_field_get_equation_param(f);

    if (eqp->isstpc == -999) {
      eqp->blencv = 0.;
      eqp->ischcv = 1;
      eqp->isstpc = 0;
      eqp->ircflu = 0;
    }

  }

  /* Additional information
     ---------------------- */

  /* Compute Ro0 based on T0 and P0
     (perfect gas law applied to air).

     Initialize RO0 with oxydant 1 which should be the dominant oxydant. */

  const int ioxy = 0;
  const int ico2 = cm->ico2 - 1;
  const int ih2o = cm->ih2o - 1;
  const int in2  = cm->in2 - 1;
  const int io2  = cm->io2 - 1;

  const cs_real_t wmolme =   (  cm->wmole[io2]  * cm->oxyo2[ioxy]
                              + cm->wmole[in2]  * cm->oxyn2[ioxy]
                              + cm->wmole[ih2o] * cm->oxyh2o[ioxy]
                              + cm->wmole[ico2] * cm->oxyco2[ioxy])
                           / (  cm->oxyo2[ioxy] + cm->oxyn2[ioxy]
                              + cm->oxyh2o[ioxy] + cm->oxyco2[ioxy]);

  fp->ro0 = fp->p0 * wmolme / (cs_physical_constants_r * fp->t0);

  // Initialization for coke density

  for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
    cm->rhock[coal_id] = cm->rho0ch[coal_id];
  }

  // Variable density and constant viscosity
  fp->irovar = 1;
  fp->ivivar = 0;

  /* Verification of user settings
     ----------------------------- */

  cs_parameters_is_in_range_double(CS_ABORT_IMMEDIATE,
                                   _("Pulverized coal combustion model setup"),
                                   "srrom",
                                   cm->srrom,
                                   0., 1.);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print pulverized combustion module options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_coal_log_setup(void)
{
  cs_coal_model_t *cm = cs_glob_coal_model;

  if (cm == nullptr)
    return;

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Pulverized coal combustion module options\n"
                  "-----------------------------------------\n\n"));

  cs_log_printf(CS_LOG_SETUP,
                _("  Time stepping relaxation coefficient\n"
                  "    rho(n+1) = srrom*rho(n) + (1-srrom)*rho(n+1)\n"
                  "    srrom: %14.5e\n"),
                cm->srrom);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add variable fields for pulverized coal combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_add_variable_fields(void)
{
  cs_coal_model_t *cm = cs_glob_coal_model;

  // Key id for drift scalar
  const int keydri = cs_field_key_id("drift_scalar_model");

  // Key ids for clipping
  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");

  // Key id for the diffusivity
  const int kvisl0 = cs_field_key_id("diffusivity_ref");

  // Drift flag
  int iscdri = (cm->idrift >= 1) ? 1 : 0;

  /* Definition of fields
     -------------------- */

  {
    // FIXME enth_bulk?
    int f_id = cs_variable_field_create("enthalpy", "Enthalpy",
                                        CS_MESH_LOCATION_CELLS, 1);
    cs_field_t *f = cs_field_by_id(f_id);
    cs_field_pointer_map(CS_ENUMF_(h), f);
    cs_add_model_field_indexes(f->id);

    cs_field_set_key_double(f, kscmin, -cs_math_big_r);
    cs_field_set_key_double(f, kscmax, cs_math_big_r);

    cs_field_set_key_double(f, kvisl0, 4.24e-5);

    /* set thermal model */
    cs_thermal_model_t *thermal = cs_get_glob_thermal_model();
    thermal->thermal_variable = CS_THERMAL_MODEL_ENTHALPY;
  }

  // Dispersed phase variables
  //--------------------------

  // NB: 'c' stands for continuous <> 'p' stands for particles

  // Drift: create additional mass flux for first element of eah class.
  if (iscdri != 0)
    iscdri = iscdri | CS_DRIFT_SCALAR_ADD_DRIFT_FLUX;

  // Number of particles of the class per kg of air-coal mixture (bulk)
  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cs_field_t *f = _add_coal_variable("n_p", "Np",
                                       class_id, class_id, iscdri);
    cm->inp[class_id] = f->id;
    cs_field_set_key_double(f, kscmin, 0);
    cs_field_set_key_double(f, kscmax, cs_math_infinite_r);
  }

  // Drift: do not create additional convective flux for
  // following elements of eah class (but use the one of the class).
  if (iscdri != 0)
    iscdri = iscdri & (~(CS_DRIFT_SCALAR_ADD_DRIFT_FLUX));

  // Mass fraction of reactive coal of the class icla per kg of bulk
  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cs_field_t *f = _add_coal_variable("x_p_coal", "Xp_Ch",
                                       class_id, class_id, iscdri);
    cm->ixch[class_id] = f->id;
    cs_field_set_key_double(f, kscmin, 0);
    cs_field_set_key_double(f, kscmax, 1);
  }

  // Mass fraction of char (coke in French) of the class per kg of bulk
  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cs_field_t *f = _add_coal_variable("x_p_char", "Xp_Ck",
                                       class_id, class_id, iscdri);
    cm->ixck[class_id] = f->id;
    cs_field_set_key_double(f, kscmin, 0);
    cs_field_set_key_double(f, kscmax, 1);
  }

  // Mass fraction of water (within the particle) of the class per kg of bulk
  if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      cs_field_t *f = _add_coal_variable("x_p_wt", "Xp_wt",
                                         class_id, class_id, iscdri);
      cm->ixwt[class_id] = f->id;
      cs_field_set_key_double(f, kscmin, 0);
      cs_field_set_key_double(f, kscmax, 1);
    }
  }

  // Enthalpy of the class per kg of bulk
  // (product of the mass fraction of the class by massic enthalpy of the class).
  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cs_field_t *f = _add_coal_variable("x_p_h", "Xp_Ent",
                                       class_id, class_id, iscdri);
    cm->ih2[class_id] = f->id;
    cs_field_set_key_double(f, kscmin, -cs_math_big_r);
    cs_field_set_key_double(f, kscmax, cs_math_big_r);
  }

  // Age of the class icla time the Np (number of particle per kg of bulk)
  if (cm->idrift >= 1) {
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      cs_field_t *f = _add_coal_variable("n_p_age", "Np_Age",
                                         class_id, class_id, iscdri);
      // TODO: test below used to reproduce previous behavior; is it desired ?
      if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING)
        cs_field_set_key_double(f, kscmin, 0);
      cs_field_set_key_double(f, kscmax, cs_math_big_r);
    }
  }

  // Particle velocities (when they are transported)

  if (cm->idrift >= 1) {
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      _add_coal_variable("v_x_p", "Vp_X", class_id, class_id, iscdri);
      _add_coal_variable("v_y_p", "Vp_Y", class_id, class_id, iscdri);
      _add_coal_variable("v_z_p", "Vp_Z", class_id, class_id, iscdri);
    }
  }

  // Continuous phase variables
  //---------------------------

  // NB: 'c' stands for continuous <> 'p' stands for particles

  // Field gas enthalpy

  // Enthalpy of the gas phase per kg of bulk
  // (The gas phase is a class with a negative icla, Enthalpy of the class
  //  is the product of the mass fraction of the class
  //  by massic enthalpy of the class).

  {
    cs_field_t *f = _add_coal_variable("x_c_h", "Xc_Ent", -1, -1, iscdri);
    cm->ihgas = f->id;

    // The first gas scalar contains the drift flux
    if (iscdri != 0) {
      iscdri = iscdri | CS_DRIFT_SCALAR_ADD_DRIFT_FLUX;
      cs_field_set_key_int(f, keydri, iscdri);
    }
  }

  // Drift: do not create additional convective fluxes
  if (iscdri != 0)
    iscdri = iscdri & (~(CS_DRIFT_SCALAR_ADD_DRIFT_FLUX));

  // Light (F8) and heavy (F9) volatile matter

  // Mass of the mean value of the tracer 1 (representing the light
  // volatiles released by the coal icha) divided by the mass of bulk
  // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)

  for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
    cs_field_t *f = _add_coal_variable("fr_mv1", "Fr_mv1",
                                       coal_id, -1, iscdri);
    cm->if1m[coal_id] = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);
  }

  // Mass of the mean value of the tracer 2 (representing the heavy
  // volatiles released by the coal icha) divided by the mass of bulk
  // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)

  for (int coal_id = 0; coal_id < cm->n_coals; coal_id++) {
    cs_field_t *f = _add_coal_variable("fr_mv2", "Fr_mv2",
                                       coal_id, -1, iscdri);
    cm->if2m[coal_id] = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);
  }

  // Mass of the Oxydant 2 divided by the mass of bulk
  // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  if (cm->noxyd >= 2) {
    cs_field_t *f = _add_coal_variable("fr_oxyd2", "FR_OXYD2", -1, -1, iscdri);
    cm->if4m = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);
  }

  // Mass of the Oxydant 3 divided by the mass of bulk
  // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  if (cm->noxyd >= 3) {
    cs_field_t *f = _add_coal_variable("fr_oxyd3", "FR_OXYD3", -1, -1, iscdri);
    cm->if5m = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);
  }

  // Mass of the water from coal drying divided by the mass of bulk
  // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {
    cs_field_t *f = _add_coal_variable("fr_h2o", "FR_H2O", -1, -1, iscdri);
    cm->if6m = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);
  }

  // Mass of the Carbon from coal oxydized by O2 divided by the mass of bulk
  // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  {
    cs_field_t *f = _add_coal_variable("fr_het_o2", "FR_HET_O2", -1, -1, iscdri);
    cm->if7m = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);
  }

  // Mass of the Carbon from coal gasified by CO2 divided by the mass of bulk
  // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  if (cm->ihtco2 == 1) {
    cs_field_t *f = _add_coal_variable("fr_het_co2", "FR_HET_CO2", -1, -1, iscdri);
    cm->if8m = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);
  }

  // Mass of the Carbon from coal gasified by H2O divided by the mass of bulk
  // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  if (cm->ihth2o == 1) {
    cs_field_t *f = _add_coal_variable("fr_het_h2o", "FR_HET_H2O", -1, -1, iscdri);
    cm->if9m = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);
  }

  // Variance
  {
    cs_field_t *f = _add_coal_variable("f1f2_variance", "Var_F1F2", -1, -1, iscdri);
    cm->ifvp2m = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 0.25);
  }

  // Mass of the Carbon dioxyde (CO or CO2) divided by the mass of bulk
  // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  //FIXME check for the oxycombustion, it would be more relevant to track CO
  if (cm->ieqco2 >= 1) {
    cs_field_t *f = _add_coal_variable("x_c_co2", "Xc_CO2", -1, -1, iscdri);
    cm->iyco2 = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);
  }

  if (cm->ieqnox == 1) {
    cs_field_t *f;

    // Mass of the HCN divided by the mass of bulk
    // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
    f = _add_coal_variable("x_c_hcn", "Xc_HCN", -1, -1, iscdri);
    cm->iyhcn = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);

    // Mass of the NH3 divided by the mass of bulk
    // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
    f = _add_coal_variable("x_c_nh3", "Xc_NH3", -1, -1, iscdri);
    cm->iynh3 = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);

    // Mass of the NO divided by the mass of bulk
    // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
    f = _add_coal_variable("x_c_no", "Xc_NO", -1, -1, iscdri);
    cm->iyno = f->id;
    cs_field_set_key_double(f, kscmin, 0.);
    cs_field_set_key_double(f, kscmax, 1.);

    // Enthalpy of the oxydizer times the fraction of gas
    // divided by the mass of bulk
    // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
    f = _add_coal_variable("x_c_h_ox", "Xc_Ent_Ox", -1, -1, iscdri);
    cm->ihox = f->id;
    cs_field_set_key_double(f, kscmin, -cs_math_big_r);
    cs_field_set_key_double(f, kscmax, cs_math_big_r);
  }

  // Age of bulk
  // FIXME give the possibility of having the age separately
  if (iscdri != 0) {
    int f_id = cs_variable_field_create("age", "Age",
                                        CS_MESH_LOCATION_CELLS, 1);
    cs_field_t *f = cs_field_by_id(f_id);
    cs_add_model_field_indexes(f->id);

    cs_field_set_key_double(f, kscmin, 0);
    cs_field_set_key_double(f, kscmax, cs_math_big_r);
  }

  /* Physical properties
     ------------------- */

  // Although we are in enthalpy formulation, we keep Cp constant.

  cs_fluid_properties_t *fluid_props = cs_get_glob_fluid_properties();
  fluid_props->icp = -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add property fields for pulverized coal combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_add_property_fields(void)
{
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");
  const int keylbl = cs_field_key_id("label");
  int iopchr = CS_POST_ON_LOCATION | CS_POST_MONITOR;
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;

  cs_coal_model_t *cm = cs_glob_coal_model;

  cs_field_t *f;

  // Definition of pointers related to state variables

  // Continuous phase (gaseous mix)
  f = _add_coal_property("temperature", "T_Gas", -1, 1);
  cs_field_pointer_map(CS_ENUMF_(t), f);
  cm->irom1 = _add_coal_property("rho_gas", "Rho_Gas", -1, 1)->id;

  // Gas mixture fractions
  cm->iym1[0]  = _add_coal_property("ym_chx1m", "Ym_CHx1m", -1, 1)->id;
  cm->iym1[1]  = _add_coal_property("ym_chx2m", "Ym_CHx2m", -1, 1)->id;
  cm->iym1[2]  = _add_coal_property("ym_co",    "Ym_CO",    -1, 1)->id;
  cm->iym1[3]  = _add_coal_property("ym_h2s",   "Ym_H2S",   -1, 1)->id;
  cm->iym1[4]  = _add_coal_property("ym_h2",    "Ym_H2",    -1, 1)->id;
  cm->iym1[5]  = _add_coal_property("ym_hcn",   "Ym_HCN",   -1, 1)->id;
  cm->iym1[6]  = _add_coal_property("ym_nh3",   "Ym_NH3",   -1, 1)->id;
  cm->iym1[7]  = _add_coal_property("ym_o2",    "Ym_O2",    -1, 1)->id;
  cm->iym1[8]  = _add_coal_property("ym_co2",   "Ym_CO2",   -1, 1)->id;
  cm->iym1[9]  = _add_coal_property("ym_h2o",   "Ym_H2O",   -1, 1)->id;
  cm->iym1[10] = _add_coal_property("ym_so2",   "Ym_SO2",   -1, 1)->id;
  cm->iym1[11] = _add_coal_property("ym_n2",    "Ym_N2",    -1, 1)->id;

  // Algebraic variables specific to gas - particles suspension
  f = _add_coal_property("xm", "Xm", -1, 1);
  cm->immel = f->id;
  cs_field_set_key_int(f, keylog, 0);
  cs_field_set_key_int(f, keyvis, 0);

  // Algebraic variables specific to continuous phase
  if (cm->ieqnox == 1) {
    cm->ighcn1 = _add_coal_property("exp1",      "EXP1",      -1, 1)->id;
    cm->ighcn2 = _add_coal_property("exp2",      "EXP1",      -1, 1)->id;
    cm->ignoth = _add_coal_property("exp3",      "EXP3",      -1, 1)->id;
    cm->ignh31 = _add_coal_property("exp4",      "EXP4",      -1, 1)->id;
    cm->ignh32 = _add_coal_property("exp5",      "EXP5",      -1, 1)->id;
    cm->ifhcnd = _add_coal_property("f_hcn_dev", "F_HCN_DEV", -1, 1)->id;
    cm->ifhcnc = _add_coal_property("f_hcn_het", "F_HCN_HET", -1, 1)->id;
    cm->ifnh3d = _add_coal_property("f_nh3_dev", "F_NH3_DEV", -1, 1)->id;
    cm->ifnh3c = _add_coal_property("f_nh3_het", "F_NH3_HET", -1, 1)->id;
    cm->ifnohc = _add_coal_property("f_no_hcn",  "F_NO_HCN", -1, 1)->id;
    cm->ifnonh = _add_coal_property("f_no_nh3",  "F_NO_NH3", -1, 1)->id;
    cm->ifnoch = _add_coal_property("f_no_het",  "F_NO_HET", -1, 1)->id;
    cm->ifnoth = _add_coal_property("f_no_the",  "F_NO_THE", -1, 1)->id;
    cm->icnohc = _add_coal_property("c_no_hcn",  "C_NO_HCN", -1, 1)->id;
    cm->icnonh = _add_coal_property("c_no_nh3",  "C_NO_NH3", -1, 1)->id;
    cm->ifhcnr = _add_coal_property("f_hcn_rb",  "F_HCN_RB", -1, 1)->id;
    cm->icnorb = _add_coal_property("c_no_rb",   "C_NO_RB",  -1, 1)->id;
    cm->igrb   = _add_coal_property("exp_rb",    "Exp_RB",   -1, 1)->id;
  }

  // Dispersed phase (particle classes)

  // NB: 'c' stands for continuous <> 'p' stands for particles

  // Temperature of particle class icla
  // NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  for (int class_id = 0; class_id < cm->nclacp; class_id++)
    cm->itemp2[class_id] = _add_coal_property("t_p_", "Tp_", class_id, 1)->id;

  for (int class_id = 0; class_id < cm->nclacp; class_id++)
    cm->ix2[class_id] = _add_coal_property("x_p_", "Xp_", class_id, 1)->id;

  for (int class_id = 0; class_id < cm->nclacp; class_id++)
    cm->irom2[class_id] = _add_coal_property("rho_p_", "Rhop_", class_id, 1)->id;

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cm->idiam2[class_id]
      = _add_coal_property("diam_p_", "Diamp_", class_id, 1)->id;
  }

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cm->igmdch[class_id]
      = _add_coal_property("dissapear_rate_p_", "D_Rate_Coal", class_id, 1)->id;
  }

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cm->igmdv1[class_id]
      = _add_coal_property("m_transfer_v1_p_", "D_V1_Coal", class_id, 1)->id;
  }

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cm->igmdv2[class_id]
      = _add_coal_property("m_transfer_v2_p_", "D_V2_Coal", class_id, 1)->id;
  }

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cm->igmhet[class_id]
      = _add_coal_property("het_ts_o2_p_", "Het_TS_O2_Coal", class_id, 1)->id;
  }

  for (int class_id = 0; class_id < cm->nclacp; class_id++) {
    cm->igmtr[class_id]
      = _add_coal_property("imp_m_transfer_to_g_p_", "Implicit_Mass_transfer",
                           class_id, 1)->id;
  }

  if (cm->idrift >= 1) {
    const int keyccl = cs_field_key_id("scalar_class");

    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      char name[64], label[64];
      int class_num = class_id+1;

      // Age of the particle class
      snprintf(name, 64, "age_p_%02d", class_num); name[63] = '\0';
      snprintf(label, 64, "Agep_%02d", class_num); label[63] = '\0';
      f = cs_field_create(name, field_type, CS_MESH_LOCATION_CELLS, 1, false);
      cs_field_set_key_str(f, keylbl, label);
      cs_field_set_key_int(f, keyccl, class_num);
      cs_field_set_key_int(f, keyvis, iopchr);

      // Limit velocity
      snprintf(name, 64, "vg_lim_p_%02d", class_num); name[63] = '\0';
      f = cs_field_create(name, field_type, CS_MESH_LOCATION_CELLS, 3, false);
      cs_field_set_key_int(f, keyccl, class_num);
      cs_field_set_key_int(f, keyvis, iopchr);
      cs_field_set_key_int(f, keylog, 1);

      snprintf(name, 64, "vg_p_%02d", class_num); name[63] = '\0';
      f = cs_field_create(name, field_type, CS_MESH_LOCATION_CELLS, 3, false);
      cs_field_set_key_int(f, keyccl, class_num);
      cs_field_set_key_int(f, keyvis, iopchr);
      cs_field_set_key_int(f, keylog, 1);

      // Additional drift velocity for the particle class
      snprintf(name, 64, "vd_p_%02d", class_num); name[63] = '\0';
      f = cs_field_create(name, field_type, CS_MESH_LOCATION_CELLS, 3, false);
      cs_field_set_key_int(f, keyccl, class_num);
      cs_field_set_key_int(f, keyvis, iopchr);
      cs_field_set_key_int(f, keylog, 1);
    }
  }

  if (cm->ihtco2 == 1) {
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      cm->ighco2[class_id]
        = _add_coal_property("het_ts_co2_p", "Het_TS_CO2_p", class_id, 1)->id;
    }
  }

  if (cm->ihth2o == 1) {
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      cm->ighh2o[class_id]
        = _add_coal_property("het_ts_h2o_p", "Het_TS_H2O_p", class_id, 1)->id;
    }
  }

  if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {
    for (int class_id = 0; class_id < cm->nclacp; class_id++) {
      cm->igmsec[class_id]
        = _add_coal_property("dry_ts_p", "Dry_TS_p", class_id, 1)->id;
      // FIXME is it a Source term?
    }
  }

  if (cm->idrift >= 1) {
    // Additional fields for drift velocity for the gas
    f = cs_field_create("vd_c",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,
                        false);
    cs_field_set_key_int(f, keyvis, iopchr);
    cs_field_set_key_int(f, keylog, 1);
  }

  // Mass fraction of the continuous phase (X1)
  cs_field_create("x_c", field_type, CS_MESH_LOCATION_CELLS, 1, false);

  // Mass fraction of the continuous phase (X1) BOUNDARY VALUE
  cs_field_create("b_x_c", field_type, CS_MESH_LOCATION_BOUNDARY_FACES,
                  1, false);

  // Explicit interfacial source terms for x1 h1 (deduced from thoses of x2 h2)
  cs_field_create("x_h_c_exp_st", field_type, CS_MESH_LOCATION_CELLS, 1, false);

  // Implicit interfacial source terms for x1 h1 (deduced from thoses of x2 h2)
  cs_field_create("x_h_c_imp_st", field_type, CS_MESH_LOCATION_CELLS, 1, false);

  /* Bulk
     ---- */

  cm->ibcarbone = _add_coal_property("x_carbon", "Z_Carbon", -1, 1)->id;
  cm->iboxygen  = _add_coal_property("x_oxygen", "Z_Oxygen", -1, 1)->id;
  cm->ibhydrogen = _add_coal_property("x_hydrogen", "Z_Hydrogen", -1, 1)->id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Take in account the radiative source terms in the particle equation
 *        of a given class for pulverized coal flame.
 *
 * \param[in]      f       pointer to scalar field
 * \param[in, out] smbrs   right and side (explicit ST contribution)
 * \param[in, out] rovsdt  system diagonal (implicit ST contribution)
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_rad_transfer_st(const cs_field_t  *f,
                        cs_real_t         *smbrs,
                        cs_real_t         *rovsdt)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_vol;

  /* Initialization
   * -------------- */

  const int keyccl = cs_field_key_id("scalar_class");
  const int numcla = cs_field_get_key_int(f, keyccl);
  const int ipcl   = 1 + numcla;

  /* Radiative source terms
   * ---------------------- */

  char f_name[64];

  snprintf(f_name, 63, "rad_st_implicit_%02d", ipcl);
  cs_real_t *cpro_tsri = cs_field_by_name(f_name)->val;

  snprintf(f_name, 63, "rad_st_%02d", ipcl);
  cs_real_t *cpro_tsre = cs_field_by_name(f_name)->val;

  snprintf(f_name, 63, "x_p_%02d", numcla);
  const cs_real_t *cval_x_p = cs_field_by_name(f_name)->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cpro_tsri[c_id] = cs::max(-cpro_tsri[c_id], 0.);

    if (cval_x_p[c_id] > cs_math_epzero) {
      /* Explicit part */
      smbrs[c_id] += cpro_tsre[c_id]*cell_f_vol[c_id]*cval_x_p[c_id];

      /* Implicit part */
      rovsdt[c_id] += cpro_tsri[c_id]*cell_f_vol[c_id];
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
