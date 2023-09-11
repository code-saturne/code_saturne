/*============================================================================
 * Base electrical model data.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

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
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_file.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_math.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_location.h"
#include "cs_time_step.h"
#include "cs_parameters.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_field_operator.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_thermal_model.h"
#include "cs_turbulence_model.h"
#include "cs_gui_specific_physics.h"
#include "cs_gui_util.h"
#include "cs_post.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_elec_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_elec_model.c
 *
 * \brief Base electrical model data.
 *
 * Please refer to the
 * <a href="../../theory.pdf#electric"><b>electric arcs</b></a>
 * section of the theory guide for more informations.
*/

/*----------------------------------------------------------------------------*/

/*! \struct cs_elec_option_t

  \brief option for electric model

  \var  cs_elec_option_t::ixkabe
        Model for radiative properties
        - 0: last column read but not use
        - 1: last column : absorption coefficient
        - 2: last column : radiative ST
  \var  cs_elec_option_t::ntdcla
        First iteration to take into account restrike model
  \var  cs_elec_option_t::irestrike
        Indicate if restrike or not
  \var  cs_elec_option_t::restrike_point
        Coordinates for restrike point
  \var  cs_elec_option_t::crit_reca
        Defines plane coordinates component used to calculate
        current in a plane.\n
        Useful if \ref modrec = 2.
  \var  cs_elec_option_t::ielcor
        Indicate if scaling or not.\n
        When \ref ielcor = 1, the boundary conditions for the potential
        will be tuned at each time step in order to reach a user-specified
        target dissipated power \ref puisim (Joule effect) or a user-specified
        target current intensity \ref couimp (electric arcs).\n The boundary
        condition tuning is controlled by subroutines \ref elreca or
        \ref cs_user_electric_scaling.
  \var  cs_elec_option_t::modrec
        Model for scaling
        - 1: volumic power for boundary conditions tuning,
        - 2: by plane for boundary conditions tuning,
        - 3: user function for boundary conditions tuning.
        In this case, we need to define the plane and the
        current density component used.
  \var  cs_elec_option_t::idreca
        Defines the current density component used to calculate
        current in a plane.\n
        Useful if \ref modrec = 2.
  \var  cs_elec_option_t::izreca
        Indicator for faces for scaling
  \var  cs_elec_option_t::couimp
        Imposed current.\n
        With the electric arcs module, \ref couimp is the target current
        intensity (\f$A\f$) for the calculations with boundary condition
        tuning for the potential.\n The target intensity will be reached
        if the boundary conditions are expressed using the variable
        \ref pot_diff or if the initial boundary conditions are multiplied by
        the variable \ref coejou.\n
        Useful with the electric arcs module if \ref ielcor = 1.
  \var  cs_elec_option_t::pot_diff
        Potential difference.\n
        \ref pot_diff is the potential difference (\f$V\f$) which generates
        the current (and the Joule effect) for the calculations with boundary
        conditions tuning for the potential. This value is initialised set by
        the user (\ref cs_user_parameters). It is then automatically tuned
        depending on the value of dissipated power (Joule effect module) or the
        intensity of current (electric arcs module). In order for the correct
        power or intensity to be reached, the boundary conditions for the
        potential must be expressed with \ref pot_diff . The tuning can be
        controlled in \ref cs_user_electric_scaling.\n
        Useful if \ref ielcor = 1.
  \var  cs_elec_option_t::puisim
        Imposed power.\n
        With the Joule effect module, \ref puisim is the target dissipated power ($W$)
        for the calculations with boundary condition tuning for the potential.\n
        The target power will be reached if the boundary conditions are expressed
        using the variable \ref pot_diff or if the initial boundary conditions are
        multiplied by the variable \ref coejou .
        Useful with the Joule effect module if \ref ielcor = 1.
  \var  cs_elec_option_t::coejou
        coefficient for scaling
  \var  cs_elec_option_t::elcou
        current in scaling plane
  \var  cs_elec_option_t::srrom
        Sub-relaxation coefficient for the density, following the formula:
        \f$\rho^{n+1}$\,=\,srrom\,$\rho^n$+(1-srrom)\,$\rho^{n+1}\f$
        Hence, with a zero value, there is no sub-relaxation.
*/

/*! \struct cs_data_joule_effect_t
²  \brief Structure to read transformer parameters in dp_ELE

  \var  cs_data_joule_effect_t::nbelec
        transformer number
  \var  cs_data_joule_effect_t::ielecc
  \var  cs_data_joule_effect_t::ielect
  \var  cs_data_joule_effect_t::ielecb
  \var  cs_data_joule_effect_t::nbtrf
  \var  cs_data_joule_effect_t::ntfref
  \var  cs_data_joule_effect_t::ibrpr
  \var  cs_data_joule_effect_t::ibrsec
  \var  cs_data_joule_effect_t::tenspr
  \var  cs_data_joule_effect_t::rnbs
  \var  cs_data_joule_effect_t::zr
  \var  cs_data_joule_effect_t::zi
  \var  cs_data_joule_effect_t::uroff
  \var  cs_data_joule_effect_t::uioff
*/

/*! \struct cs_data_elec_t

  \brief physical properties for electric model descriptor.

  \var  cs_data_elec_t::ngaz
        number of gas in electrical data file
  \var  cs_data_elec_t::npoint
        number of point in electrical data file for each gas
  \var  cs_data_elec_t::th
        temperature values
  \var  cs_data_elec_t::ehgaz
        enthalpy values
  \var  cs_data_elec_t::rhoel
        density values
  \var  cs_data_elec_t::cpel
        specific heat values
  \var  cs_data_elec_t::sigel
        electric conductivity values
  \var  cs_data_elec_t::visel
        dynamic viscosity
  \var  cs_data_elec_t::xlabel
        thermal conductivity
  \var  cs_data_elec_t::xkabel
        absorption coefficent
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define LG_MAX 1000

/*============================================================================
 * Type definitions
 *============================================================================*/

static cs_elec_option_t  _elec_option = {.ixkabe = -1,
                                         .ntdcla = -1,
                                         .irestrike = -1,
                                         .restrike_point = {0., 0., 0.},
                                         .crit_reca = {0, 0, 0, 0, 0},
                                         .ielcor = -1,
                                         .modrec = -1,
                                         .idreca = -1,
                                         .izreca = NULL,
                                         .couimp = 0.,
                                         .pot_diff = 0.,
                                         .puisim = 0.,
                                         .coejou = 0.,
                                         .elcou = 0.,
                                         .srrom = 0.};

static cs_data_elec_t  _elec_properties = {.ngaz = 0,
                                           .npoint = 0,
                                           .th = NULL,
                                           .ehgaz = NULL,
                                           .rhoel = NULL,
                                           .cpel = NULL,
                                           .sigel = NULL,
                                           .visel = NULL,
                                           .xlabel = NULL,
                                           .xkabel = NULL};

static cs_data_joule_effect_t  *_transformer     = NULL;

const cs_elec_option_t        *cs_glob_elec_option = NULL;
const cs_data_elec_t          *cs_glob_elec_properties = NULL;
const cs_data_joule_effect_t  *cs_glob_transformer     = NULL;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*!
 * vacuum magnetic permeability constant (H/m). (= 1.2566e-6)
 *
 */
const cs_real_t cs_elec_permvi = 1.2566e-6;

/*!
 * vacuum permittivity constant (F/m). (= 8.854e-12)
 *
 */
const cs_real_t cs_elec_epszer = 8.854e-12;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

void
cs_f_elec_model_get_pointers(int     **ielcor,
                             double  **pot_diff,
                             double  **coejou,
                             double  **elcou);

/*----------------------------------------------------------------------------
 * Get pointers to members of the global electric model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ielcor         --> pointer to cs_glob_elec_option->ielcor
 *   pot_diff       --> pointer to cs_glob_elec_option->pot_diff
 *   coejou         --> pointer to cs_glob_elec_option->coejou
 *   elcou          --> pointer to cs_glob_elec_option->elcou
 *----------------------------------------------------------------------------*/

void
cs_f_elec_model_get_pointers(int     **ielcor,
                             double  **pot_diff,
                             double  **coejou,
                             double  **elcou)
{
  *ielcor           = &(_elec_option.ielcor);
  *pot_diff         = &(_elec_option.pot_diff);
  *coejou           = &(_elec_option.coejou);
  *elcou            = &(_elec_option.elcou);
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

static void
_cs_electrical_model_verify(void)
{
  bool verif = true;

  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];
  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  if (ielarc != -1 && ielarc !=  2)
    bft_error(__FILE__, __LINE__, 0,
              _("Error for electric arc model\n"
                "only choice -1 or 2 are permitted yet\n"
                "model selected : \"%i\";\n"),
              ielarc);

  if (   ieljou != -1 && ieljou !=  1 && ieljou !=  2 && ieljou !=  3
      && ieljou !=  4)
    bft_error(__FILE__, __LINE__, 0,
              _("Error for joule model\n"
                "only choice -1, 1, 2, 3 or 4 are permitted yet\n"
                "model selected : \"%i\";\n"),
              ieljou);

  /* options */
  if (cs_glob_elec_option->ielcor != 0 && cs_glob_elec_option->ielcor != 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Error for scaling model\n"
                "only choice -1 or 2 are permitted yet\n"
                "model selected : \"%i\";\n"),
              cs_glob_elec_option->ielcor);

  if (cs_glob_elec_option->ielcor == 1) {
    if (ielarc > 0) {
      if (cs_glob_elec_option->couimp < 0.) {
        bft_printf("value for COUIMP must be strictly positive\n");
        verif = false;
      }
      if (cs_glob_elec_option->pot_diff < 0.) {
        bft_printf("value for DPOT must be strictly positive\n");
        verif = false;
      }
    }
    if (ieljou > 0) {
      if (cs_glob_elec_option->puisim < 0.) {
        bft_printf("value for PUISIM must be strictly positive\n");
        verif = false;
      }
      if (cs_glob_elec_option->coejou < 0.) {
        bft_printf("value for COEJOU must be strictly positive\n");
        verif = false;
      }
      if (cs_glob_elec_option->pot_diff < 0.) {
        bft_printf("value for DPOT must be strictly positive\n");
        verif = false;
      }
    }
  }

  if (!verif) {
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid or incomplete calculation parameter\n"
                "Verify parameters\n"));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers for electric arcs
 *
 * \param[in]  n_gasses    number of gasses
 */
/*----------------------------------------------------------------------------*/

static void
_field_pointer_map_electric_arcs(int  n_gasses)
{
  char s[64];

  cs_field_pointer_map(CS_ENUMF_(h),
                       cs_field_by_name_try("enthalpy"));

  cs_field_pointer_map(CS_ENUMF_(potr), cs_field_by_name_try("elec_pot_r"));
  cs_field_pointer_map(CS_ENUMF_(poti), cs_field_by_name_try("elec_pot_i"));

  cs_field_pointer_map(CS_ENUMF_(potva), cs_field_by_name_try("vec_potential"));

  for (int i = 0; i < n_gasses - 1; i++) {
    snprintf(s, 63, "esl_fraction_%02d", i+1); s[63] = '\0';
    cs_field_pointer_map_indexed(CS_ENUMF_(ycoel), i, cs_field_by_name_try(s));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers properties for electric arcs
 */
/*----------------------------------------------------------------------------*/

static void
_field_pointer_properties_map_electric_arcs(void)
{
  cs_field_pointer_map(CS_ENUMF_(t),
                       cs_field_by_name_try("temperature"));

  cs_field_pointer_map(CS_ENUMF_(joulp),
                       cs_field_by_name_try("joule_power"));
  cs_field_pointer_map(CS_ENUMF_(radsc),
                       cs_field_by_name_try("radiation_source"));
  cs_field_pointer_map(CS_ENUMF_(elech),
                       cs_field_by_name_try("elec_charge"));

  cs_field_pointer_map(CS_ENUMF_(curre),
                       cs_field_by_name_try("current_re"));
  cs_field_pointer_map(CS_ENUMF_(curim),
                       cs_field_by_name_try("current_im"));
  cs_field_pointer_map(CS_ENUMF_(laplf),
                       cs_field_by_name_try("laplace_force"));
  cs_field_pointer_map(CS_ENUMF_(magfl),
                       cs_field_by_name_try("magnetic_field"));
  cs_field_pointer_map(CS_ENUMF_(elefl),
                       cs_field_by_name_try("electric_field"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the imaginary potential gradient at specified cells.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        pointer to associated mesh structure
 *                               (to be cast as cs_mesh_t *) for interior
 *                               faces or vertices, unused otherwise
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_pot_gradient_im_f(int               location_id,
                   cs_lnum_t         n_elts,
                   const cs_lnum_t  *elt_ids,
                   void             *input,
                   void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_CELLS);

  cs_real_3_t *v = vals;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_field_t *f = cs_field_by_name("elec_pot_i");

  cs_real_3_t *grad;
  BFT_MALLOC(grad, m->n_cells_with_ghosts, cs_real_3_t);

  cs_field_gradient_scalar(f, false, 1, grad);

  if (elt_ids != NULL) {
    for (cs_lnum_t idx = 0; idx <  n_elts; idx++) {
      cs_lnum_t i = elt_ids[idx];
      for (cs_lnum_t j = 0; j < 3; j++)
        v[idx][j] = grad[i][j];
    }
  }

  else {
    for (cs_lnum_t i = 0; i <  n_elts; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        v[i][j] = grad[i][j];
    }
  }

  BFT_FREE(grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the imaginary current at specified cells.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        pointer to associated mesh structure
 *                               (to be cast as cs_mesh_t *) for interior
 *                               faces or vertices, unused otherwise
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_current_im_f(int               location_id,
              cs_lnum_t         n_elts,
              const cs_lnum_t  *elt_ids,
              void             *input,
              void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_CELLS);

  cs_real_3_t *v = vals;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_field_t *f = cs_field_by_name("elec_pot_i");

  cs_real_3_t *grad;
  BFT_MALLOC(grad, m->n_cells_with_ghosts, cs_real_3_t);

  cs_field_gradient_scalar(f, false, 1, grad);

  const int kivisl = cs_field_key_id("diffusivity_id");
  const int diff_id = cs_field_get_key_int(f, kivisl);

  if (diff_id > -1) {
    const cs_real_t *cvisii = cs_field_by_id(diff_id)->val;

    if (elt_ids != NULL) {
      for (cs_lnum_t idx = 0; idx <  n_elts; idx++) {
        cs_lnum_t i = elt_ids[idx];
        for (cs_lnum_t j = 0; j < 3; j++)
          v[idx][j] = -cvisii[i] * grad[i][j];
      }
    }

    else {
      for (cs_lnum_t i = 0; i <  n_elts; i++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          v[i][j] = -cvisii[i] * grad[i][j];
      }
    }
  }

  else {
    const int kvisls0 = cs_field_key_id("diffusivity_ref");
    const double visls_0 = cs_field_get_key_double(f, kvisls0);

    if (elt_ids != NULL) {
      for (cs_lnum_t idx = 0; idx <  n_elts; idx++) {
        cs_lnum_t i = elt_ids[idx];
        for (cs_lnum_t j = 0; j < 3; j++)
          v[idx][j] = -visls_0 * grad[i][j];
      }
    }

    else {
      for (cs_lnum_t i = 0; i <  n_elts; i++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          v[i][j] = -visls_0 * grad[i][j];
      }
    }
  }

  BFT_FREE(grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the module of the complex potential at specified cells.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        pointer to associated mesh structure
 *                               (to be cast as cs_mesh_t *) for interior
 *                               faces or vertices, unused otherwise
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_pot_module_f(int               location_id,
              cs_lnum_t         n_elts,
              const cs_lnum_t  *elt_ids,
              void             *input,
              void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_CELLS);

  cs_real_t *v = vals;

  const cs_real_t *cpotr = cs_field_by_name("elec_pot_r")->val;
  const cs_real_t *cpoti = cs_field_by_name("elec_pot_i")->val;

  if (elt_ids != NULL) {
    for (cs_lnum_t idx = 0; idx <  n_elts; idx++) {
      cs_lnum_t i = elt_ids[idx];
      v[idx] = sqrt(cpotr[i]*cpotr[i] + cpoti[i]*cpoti[i]);
    }
  }

  else {
    for (cs_lnum_t i = 0; i <  n_elts; i++) {
      v[i] = sqrt(cpotr[i]*cpotr[i] + cpoti[i]*cpoti[i]);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the argument of the complex potential at specified cells.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        pointer to associated mesh structure
 *                               (to be cast as cs_mesh_t *) for interior
 *                               faces or vertices, unused otherwise
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_pot_arg_f(int               location_id,
           cs_lnum_t         n_elts,
           const cs_lnum_t  *elt_ids,
           void             *input,
           void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_CELLS);

  cs_real_t *v = vals;

  const cs_real_t *cpotr = cs_field_by_name("elec_pot_r")->val;
  const cs_real_t *cpoti = cs_field_by_name("elec_pot_i")->val;

  cs_real_t pi_ov_4 = atan(1.);

  for (cs_lnum_t idx = 0; idx <  n_elts; idx++) {
    cs_lnum_t i = (elt_ids != NULL) ? elt_ids[idx] : idx;

    if (cpotr[i] > 0.)
      v[idx] = atan(cpoti[i] / cpotr[i]);
    else if (cpotr[i] < 0.){
      if (cpoti[i] > 0.)
        v[idx] = 4.*atan(1.) + atan(cpoti[i] / cpotr[i]);
      else
        v[idx] = -4.*atan(1.) + atan(cpoti[i] / cpotr[i]);
    }
    else {
      v[idx] = 2.*atan(1.);
    }

    if (v[idx] < 0)
      v[idx] += pow(8., pi_ov_4);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------*/

void
CS_PROCF (elini1, ELINI1) (void)
{
  cs_electrical_model_specific_initialization();
}

void
CS_PROCF (elflux, ELFLUX) (int      *iappel)
{
  cs_elec_compute_fields(cs_glob_mesh,
                         *iappel);
}

void
CS_PROCF (elthht, ELTHHT) (int       *mode,
                           cs_real_t *ym,
                           cs_real_t *enthal,
                           cs_real_t *temp)
{
  if (*mode == -1)
    *enthal = cs_elec_convert_t_to_h(ym, *temp);
  else if (*mode == 1)
    *temp = cs_elec_convert_h_to_t(ym, *enthal);
  else
    bft_error(__FILE__, __LINE__, 0,
              _("electric module:\n"
                "bad value for mode (integer equal to -1 or 1 : %i here.\n"),
              *mode);
}

void
CS_PROCF (ellecd, ELLECD) (void)
{
  cs_electrical_model_initialize();
  cs_electrical_properties_read();
}

void
CS_PROCF (elphyv, ELPHYV) (void)
{
  cs_elec_physical_properties(cs_glob_domain);
}

void
CS_PROCF (eliniv, ELINIV) (int       *isuite)
{
  cs_elec_fields_initialize(cs_glob_mesh,  *isuite);
}

void
CS_PROCF (elreca, ELRECA) (cs_real_t *dt)
{
  cs_elec_scaling_function(cs_glob_mesh, cs_glob_mesh_quantities, dt);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide access to cs_elec_option
 *----------------------------------------------------------------------------*/

cs_elec_option_t *
cs_get_glob_elec_option(void)
{
  return &_elec_option;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_transformer
 *----------------------------------------------------------------------------*/

cs_data_joule_effect_t *
cs_get_glob_transformer(void)
{
  return _transformer;
}

/*----------------------------------------------------------------------------
 * Initialize structures for electrical model
 *----------------------------------------------------------------------------*/

void
cs_electrical_model_initialize(void)
{
  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];

  if (ieljou >= 3)
    BFT_MALLOC(_transformer, 1, cs_data_joule_effect_t);

  _elec_option.ixkabe    = 0;
  _elec_option.ntdcla    = 1;
  _elec_option.irestrike = 0;
  for (int i = 0; i < 3; i++)
    _elec_option.restrike_point[i] = 0.;
  _elec_option.izreca    = NULL;
  _elec_option.elcou     = 0.;
  _elec_option.ielcor    = 0;
  _elec_option.couimp    = 0.;
  _elec_option.puisim    = 0.;
  _elec_option.pot_diff  = 0.;
  _elec_option.coejou    = 1.;
  _elec_option.modrec    = 1;    /* standard model */
  _elec_option.idreca    = 3;
  _elec_option.srrom     = 0.;

  for (int i = 0; i < 3; i++)
    _elec_option.crit_reca[i] = 0.;
  _elec_option.crit_reca[4] = 0.0002;

  cs_glob_elec_option     = &_elec_option;
  cs_glob_elec_properties = &_elec_properties;
  cs_glob_transformer     = _transformer;

  cs_fluid_properties_t *fluid_properties = cs_get_glob_fluid_properties();
  fluid_properties->icp = 0;
  fluid_properties->irovar = 1;
  fluid_properties->ivivar = 1;

  cs_base_at_finalize(cs_electrical_model_finalize);
}

/*----------------------------------------------------------------------------
 * Destroy structures for electrical model
 *----------------------------------------------------------------------------*/

void
cs_electrical_model_finalize(void)
{
  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];
  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  if (ielarc > 0) {
    BFT_FREE(_elec_properties.th);
    BFT_FREE(_elec_properties.ehgaz);
    BFT_FREE(_elec_properties.rhoel);
    BFT_FREE(_elec_properties.cpel);
    BFT_FREE(_elec_properties.sigel);
    BFT_FREE(_elec_properties.visel);
    BFT_FREE(_elec_properties.xlabel);
    BFT_FREE(_elec_properties.xkabel);
  }

  if (ieljou >= 3) {
    BFT_FREE(_transformer->tenspr);
    BFT_FREE(_transformer->rnbs);
    BFT_FREE(_transformer->zr);
    BFT_FREE(_transformer->zi);
    BFT_FREE(_transformer->ibrpr);
    BFT_FREE(_transformer->ibrsec);
    BFT_FREE(_transformer->tenspr);
    BFT_FREE(_transformer->uroff);
    BFT_FREE(_transformer->uioff);
    BFT_FREE(_transformer);
  }

  BFT_FREE(_elec_option.izreca);
}

/*----------------------------------------------------------------------------
 * Specific initialization for electric arc
 *----------------------------------------------------------------------------*/

void
cs_electrical_model_specific_initialization(void)
{
  cs_field_t *f = NULL;
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  const int kvisls0 = cs_field_key_id("diffusivity_ref");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  cs_var_cal_opt_t var_cal_opt;

  /* specific initialization for field */
  f = CS_F_(potr);
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
  var_cal_opt.iconv  = 0;
  var_cal_opt.istat  = 0;
  var_cal_opt.idiff  = 1;
  var_cal_opt.idifft = 0;

  cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);

  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];
  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  if (ieljou == 2 || ieljou == 4) {
    f = CS_F_(poti);
    cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
    var_cal_opt.iconv  = 0;
    var_cal_opt.istat  = 0;
    var_cal_opt.idiff  = 1;
    var_cal_opt.idifft = 0;
    cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);
  }

  if (ielarc > 1) {
    cs_field_t  *fp = cs_field_by_name_try("vec_potential");
    cs_field_get_key_struct(fp, key_cal_opt_id, &var_cal_opt);
    var_cal_opt.iconv  = 0;
    var_cal_opt.istat  = 0;
    var_cal_opt.idiff  = 1;
    var_cal_opt.idifft = 0;
    cs_field_set_key_struct(fp, key_cal_opt_id, &var_cal_opt);
    cs_field_set_key_double(fp, kvisls0, 1.0);
  }

  /* for all specific field */
  f = CS_F_(h);
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
  var_cal_opt.blencv = 1.;
  cs_field_set_key_double(f, ksigmas, 0.7);
  cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);

  f = CS_F_(potr);
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
  var_cal_opt.blencv = 1.;
  cs_field_set_key_double(f, ksigmas, 0.7);
  cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);

  if (ieljou == 2 || ieljou == 4) {
    f = CS_F_(poti);
    cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
    var_cal_opt.blencv = 1.;
    cs_field_set_key_double(f, ksigmas, 0.7);
    cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);
  }

  if (ielarc > 1) {
    cs_field_t  *fp = cs_field_by_name_try("vec_potential");
    cs_field_get_key_struct(fp, key_cal_opt_id, &var_cal_opt);
    var_cal_opt.blencv = 1.;
    cs_field_set_key_double(f, ksigmas, 0.7);
    cs_field_set_key_struct(fp, key_cal_opt_id, &var_cal_opt);
  }

  if (cs_glob_elec_properties->ngaz > 1) {
    for (int igaz = 0; igaz < cs_glob_elec_properties->ngaz - 1; igaz++) {
      f = CS_FI_(ycoel, igaz);
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      var_cal_opt.blencv = 1.;
      cs_field_set_key_double(f, ksigmas, 0.7);
      cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);
    }
  }

  cs_gui_elec_model();
  _elec_option.pot_diff = 1000.; //FIXME

  _cs_electrical_model_verify();
}

/*----------------------------------------------------------------------------
 * Read properties file
 *----------------------------------------------------------------------------*/

void
cs_electrical_properties_read(void)
{
  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];
  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  if (ielarc <= 0 && ieljou < 3)
    return;

  char str[LG_MAX];

  if (ielarc > 0) {

    /* read local file for electric properties if present,
       default otherwise */

    FILE *file = cs_base_open_properties_data_file("dp_ELE");

    /* Position at the beginning of the file */
    fseek(file, 0, SEEK_SET);

    int nb_line_tot = 0;

    int iesp = 0;
    int it = 0;

    while (fgets(str, LG_MAX, file) != NULL) {
      nb_line_tot++;
      if (nb_line_tot < 8)
        continue;

      /* read number of fluids and number of points */
      if (nb_line_tot == 8)
        sscanf(str, "%d %d",
               &(_elec_properties.ngaz),
               &(_elec_properties.npoint));

      if (_elec_properties.ngaz <= 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("incorrect number of species \"%i\";\n"),
                  _elec_properties.ngaz);

      cs_lnum_t size =   cs_glob_elec_properties->ngaz
                       * cs_glob_elec_properties->npoint;

      if (nb_line_tot == 8) {
        BFT_MALLOC(_elec_properties.th,
                   cs_glob_elec_properties->npoint,
                   cs_real_t);
        BFT_MALLOC(_elec_properties.ehgaz,  size, cs_real_t);
        BFT_MALLOC(_elec_properties.rhoel,  size, cs_real_t);
        BFT_MALLOC(_elec_properties.cpel,   size, cs_real_t);
        BFT_MALLOC(_elec_properties.sigel,  size, cs_real_t);
        BFT_MALLOC(_elec_properties.visel,  size, cs_real_t);
        BFT_MALLOC(_elec_properties.xlabel, size, cs_real_t);
        BFT_MALLOC(_elec_properties.xkabel, size, cs_real_t);
      }

      if (nb_line_tot < 14)
        continue;

      if (nb_line_tot == 14)
        sscanf(str, "%i", &(_elec_option.ixkabe));

      if (cs_glob_elec_option->ixkabe < 0 || cs_glob_elec_option->ixkabe >= 3)
        bft_error(__FILE__, __LINE__, 0,
                  _("incorrect choice for radiative model \"%i\";\n"),
                  cs_glob_elec_option->ixkabe < 0);

      if (nb_line_tot < 22)
        continue;

      if (nb_line_tot >= 22) {
        sscanf(str, "%lf %lf %lf %lf %lf %lf %lf %lf",
               &(_elec_properties.th[it]),
               &(_elec_properties.ehgaz[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties.rhoel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties.cpel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties.sigel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties.visel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties.xlabel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties.xkabel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]));
        it++;
        if (it == cs_glob_elec_properties->npoint) {
          iesp++;
          it = 0;
        }
      }
    }

    fclose(file);
  }

#if 0
  for (int it = 0; it < cs_glob_elec_properties->npoint; it++)
    bft_printf("read dp_ELE "
               "%15.8E %15.8E %15.8E %15.8E "
               "%15.8E %15.8E %15.8E %15.8E\n",
               _elec_properties.th[it],
               _elec_properties.ehgaz[it],
               _elec_properties.rhoel[it],
               _elec_properties.cpel[it],
               _elec_properties.sigel[it],
               _elec_properties.visel[it],
               _elec_properties.xlabel[it],
               _elec_properties.xkabel[it]);
#endif

  if (ieljou >= 3) {

    /* read local file for Joule effect if present,
       default otherwise */

    FILE *file = cs_base_open_properties_data_file("dp_transformers");

    /* Position at the beginning of the file */
    fseek(file, 0, SEEK_SET);

    int nb_line_tot = 0;

    int iesp = 0;
    int it = 0;
    while (fgets(str, LG_MAX, file) != NULL) {
      nb_line_tot++;
      if (nb_line_tot == 1)
        sscanf(str, "%i", &(_transformer->ntfref));

      if (nb_line_tot < 4)
        continue;

      if (nb_line_tot == 4) {
        sscanf(str, "%i", &(_transformer->nbtrf));

        BFT_MALLOC(_transformer->tenspr,  cs_glob_transformer->nbtrf, cs_real_t);
        BFT_MALLOC(_transformer->rnbs,    cs_glob_transformer->nbtrf, cs_real_t);
        BFT_MALLOC(_transformer->zr,      cs_glob_transformer->nbtrf, cs_real_t);
        BFT_MALLOC(_transformer->zi,      cs_glob_transformer->nbtrf, cs_real_t);
        BFT_MALLOC(_transformer->ibrpr,   cs_glob_transformer->nbtrf, int);
        BFT_MALLOC(_transformer->ibrsec,  cs_glob_transformer->nbtrf, int);

        // alloc for boundary conditions
        BFT_MALLOC(_transformer->uroff,  cs_glob_transformer->nbtrf, cs_real_t);
        BFT_MALLOC(_transformer->uioff,  cs_glob_transformer->nbtrf, cs_real_t);
      }

      if (nb_line_tot > 4 && nb_line_tot <= 4 + cs_glob_transformer->nbtrf * 6) {
        it++;
        if (it == 1)
          continue;
        if (it == 2)
          sscanf(str, "%lf", &(_transformer->tenspr[iesp]));
        if (it == 3)
          sscanf(str, "%lf", &(_transformer->rnbs[iesp]));
        if (it == 4)
          sscanf(str, "%lf %lf", &(_transformer->zr[iesp]), &(_transformer->zi[iesp]));
        if (it == 5)
          sscanf(str, "%i", &(_transformer->ibrpr[iesp]));
        if (it == 6) {
          sscanf(str, "%i", &(_transformer->ibrsec[iesp]));
          it = 0;
          iesp++;
        }
      }

      if (nb_line_tot < 7 + cs_glob_transformer->nbtrf * 6)
        continue;

      if (nb_line_tot == 7 + cs_glob_transformer->nbtrf * 6) {
        sscanf(str, "%i", &(_transformer->nbelec));
        BFT_MALLOC(_transformer->ielecc,  cs_glob_transformer->nbelec, int);
        BFT_MALLOC(_transformer->ielect,  cs_glob_transformer->nbelec, int);
        BFT_MALLOC(_transformer->ielecb,  cs_glob_transformer->nbelec, int);
        iesp = 0;
      }

      if (nb_line_tot > 7 + cs_glob_transformer->nbelec * 6) {
        sscanf(str, "%i %i %i",
               &(_transformer->ielecc[iesp]),
               &(_transformer->ielect[iesp]),
               &(_transformer->ielecb[iesp]));
        iesp++;
      }
    }

    fclose(file);
  }
}

/*----------------------------------------------------------------------------
 * compute physical properties
 *----------------------------------------------------------------------------*/

void
cs_elec_physical_properties(cs_domain_t  *domain)
{
  static long ipass = 0;
  int nt_cur = cs_glob_time_step->nt_cur;
  int isrrom = 0;
  const cs_lnum_t  n_cells = domain->mesh->n_cells;
  const int kivisl = cs_field_key_id("diffusivity_id");
  int diff_id = cs_field_get_key_int(CS_F_(potr), kivisl);
  cs_field_t *c_prop = NULL;
  if (diff_id > -1)
    c_prop = cs_field_by_id(diff_id);
  ipass++;

  const cs_data_elec_t  *e_props = cs_glob_elec_properties; /* local name */

  if (nt_cur > 1 && cs_glob_elec_option->srrom > 0.)
    isrrom = 1;

  /* Joule effect (law must be specified by user) */

  int ifcvsl = cs_field_get_key_int(CS_F_(h), kivisl);
  cs_field_t *diff_th = NULL;
  if (ifcvsl >= 0)
    diff_th = cs_field_by_id(ifcvsl);

  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  /* Electric arc */

  if (ielarc > 0) {
    if (ipass == 1)
      bft_printf("electric arc module: properties read on file.\n");

    /* compute temperature from enthalpy */
    int ngaz = e_props->ngaz;
    int npt  = e_props->npoint;

    cs_real_t *ym, *yvol, *roesp, *visesp, *cpesp;
    cs_real_t *sigesp, *xlabes, *xkabes, *coef;
    BFT_MALLOC(ym,     ngaz, cs_real_t);
    BFT_MALLOC(yvol,   ngaz, cs_real_t);
    BFT_MALLOC(roesp,  ngaz, cs_real_t);
    BFT_MALLOC(visesp, ngaz, cs_real_t);
    BFT_MALLOC(cpesp,  ngaz, cs_real_t);
    BFT_MALLOC(sigesp, ngaz, cs_real_t);
    BFT_MALLOC(xlabes, ngaz, cs_real_t);
    BFT_MALLOC(xkabes, ngaz, cs_real_t);
    BFT_MALLOC(coef,   ngaz * ngaz, cs_real_t);

    int ifcsig = cs_field_get_key_int(CS_F_(potr), kivisl);

    if (ngaz == 1) {
      ym[0] = 1.;

      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        CS_F_(t)->val[iel] = cs_elec_convert_h_to_t(ym, CS_F_(h)->val[iel]);
    }
    else {

      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        ym[ngaz - 1] = 1.;

        for (int ii = 0; ii < ngaz - 1; ii++) {
          ym[ii] = CS_FI_(ycoel, ii)->val[iel];
          ym[ngaz - 1] -= ym[ii];
        }

        CS_F_(t)->val[iel] = cs_elec_convert_h_to_t(ym, CS_F_(h)->val[iel]);
      }
    }

    /* Map some fields */

    cs_real_t *cpro_absco = NULL;

    if (cs_glob_elec_option->ixkabe == 1) {
      if (CS_FI_(rad_cak, 0) != NULL)
        cpro_absco = CS_FI_(rad_cak, 0)->val;
    }

    /* Interpolate properties */

#   pragma omp parallel for
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      // temperature
      cs_real_t tp = CS_F_(t)->val[iel];

      // determine point
      int it = -1;
      if (tp <= e_props->th[0])
        it = 0;
      else if (tp >= e_props->th[npt - 1])
        it = npt - 1;
      else {
        for (int iiii = 0; iiii < npt - 1; iiii++)
          if (tp > e_props->th[iiii] &&
              tp <= e_props->th[iiii + 1]) {
            it = iiii;
            break;
          }
      }
      if (it == -1)
        bft_error(__FILE__, __LINE__, 0,
                  _("electric module: properties read on file\n"
                    "Warning: error in cs_elec_physical_properties\n"
                    "Invalid reading with temperature : %f.\n"),
                  tp);

      /* mass fraction */
      ym[ngaz - 1] = 1.;

      for (int ii = 0; ii < ngaz - 1; ii++) {
        ym[ii] = CS_FI_(ycoel, ii)->val[iel];
        ym[ngaz - 1] -= ym[ii];
      }

      /* density, viscosity, ... for each species */
      if (it == 0) {
        for (int ii = 0; ii < ngaz; ii++) {
          roesp[ii]  = e_props->rhoel[ii * (npt - 1)];
          visesp[ii] = e_props->visel[ii * (npt - 1)];
          cpesp[ii]  = e_props->cpel[ii * (npt - 1)];
          sigesp[ii] = e_props->sigel[ii * (npt - 1)];
          xlabes[ii] = e_props->xlabel[ii * (npt - 1)];

          if (cs_glob_elec_option->ixkabe > 0)
            xkabes[ii] = e_props->xkabel[ii * (npt - 1)];
        }
      }
      else if (it == npt - 1) {
        for (int ii = 0; ii < ngaz; ii++) {
          roesp[ii]  = e_props->rhoel[ii * (npt - 1) + npt - 1];
          visesp[ii] = e_props->visel[ii * (npt - 1) + npt - 1];
          cpesp[ii]  = e_props->cpel[ii * (npt - 1) + npt - 1];
          sigesp[ii] = e_props->sigel[ii * (npt - 1) + npt - 1];
          xlabes[ii] = e_props->xlabel[ii * (npt - 1) + npt - 1];

          if (cs_glob_elec_option->ixkabe > 0)
            xkabes[ii] = e_props->xkabel[ii * (npt - 1) + npt - 1];
        }
      }
      else {
        double delt = e_props->th[it + 1] - e_props->th[it];

        for (int ii = 0; ii < ngaz; ii++) {
          double alpro = (e_props->rhoel[ii * (npt - 1) + it + 1] -
                          e_props->rhoel[ii * (npt - 1) + it]) / delt;
          roesp[ii]  =   e_props->rhoel[ii * (npt - 1) + it]
                       + alpro * (tp -e_props->th[it]);

          double alpvis =  (e_props->visel[ii * (npt - 1) + it + 1]
                          - e_props->visel[ii * (npt - 1) + it]) / delt;
          visesp[ii] =   e_props->visel[ii * (npt - 1) + it]
                       + alpvis * (tp -e_props->th[it]);

          double alpcp = (e_props->cpel[ii * (npt - 1) + it + 1] -
                         e_props->cpel[ii * (npt - 1) + it]) / delt;
          cpesp[ii]  = e_props->cpel[ii * (npt - 1) + it] +
                       alpcp * (tp -e_props->th[it]);

          double alpsig = (e_props->sigel[ii * (npt - 1) + it + 1] -
                          e_props->sigel[ii * (npt - 1) + it]) / delt;
          sigesp[ii] = e_props->sigel[ii * (npt - 1) + it] +
                       alpsig * (tp -e_props->th[it]);

          double alplab = (e_props->xlabel[ii * (npt - 1) + it + 1] -
                          e_props->xlabel[ii * (npt - 1) + it]) / delt;
          xlabes[ii] = e_props->xlabel[ii * (npt - 1) + it] +
                       alplab * (tp -e_props->th[it]);

          if (cs_glob_elec_option->ixkabe > 0) {
            double alpkab = (e_props->xkabel[ii * (npt - 1) + it + 1] -
                            e_props->xkabel[ii * (npt - 1) + it]) / delt;
            xkabes[ii] = e_props->xkabel[ii * (npt - 1) + it] +
                         alpkab * (tp -e_props->th[it]);
          }
        }
      }

      /* compute density */
      double rhonp1 = 0.;

      for (int ii = 0; ii < ngaz; ii++)
        rhonp1 += ym[ii] / roesp[ii];

      rhonp1 = 1. / rhonp1;

      if (isrrom == 1)
        CS_F_(rho)->val[iel] = CS_F_(rho)->val[iel] * cs_glob_elec_option->srrom
                               + (1. - cs_glob_elec_option->srrom) * rhonp1;
      else
        CS_F_(rho)->val[iel] = rhonp1;

      for (int ii = 0; ii < ngaz; ii++) {
        yvol[ii] = ym[ii] * roesp[ii] / CS_F_(rho)->val[iel];
        if (yvol[ii] <= 0.)
          yvol[ii] = cs_math_epzero * cs_math_epzero;
      }

      /* compute molecular viscosity : kg/(m s) */
      for (int iesp1 = 0; iesp1 < ngaz; iesp1++) {
        for (int iesp2 = 0; iesp2 < ngaz; iesp2++) {
          coef[iesp1 * (ngaz - 1) + iesp2]
            = 1. +   sqrt(visesp[iesp1] / visesp[iesp2])
                   * sqrt(sqrt(roesp[iesp2] / roesp[iesp1]));
          coef[iesp1 * (ngaz - 1) + iesp2] *= coef[iesp1 * (ngaz - 1) + iesp2];
          coef[iesp1 * (ngaz - 1) + iesp2] /=    (sqrt(1. + roesp[iesp1]
                                               / roesp[iesp2]) * sqrt(8.));
        }
      }

      CS_F_(mu)->val[iel] = 0.;

      for (int iesp1 = 0; iesp1 < ngaz; iesp1++) {
        if (yvol[iesp1] > 1e-30) {
          double somphi = 0.;
          for (int iesp2 = 0; iesp2 < ngaz; iesp2++) {
            if (iesp1 != iesp2)
              somphi +=   coef[iesp1 * (ngaz - 1) + iesp2]
                        * yvol[iesp2] / yvol[iesp1];
          }

          CS_F_(mu)->val[iel] += visesp[iesp1] / (1. + somphi);
        }
      }

      /* compute specific heat : J/(kg degres) */
      if (cs_glob_fluid_properties->icp > 0) {
        CS_F_(cp)->val[iel] = 0.;
        for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
          CS_F_(cp)->val[iel] += ym[iesp1] * cpesp[iesp1];
      }

      /* compute Lambda/Cp : kg/(m s) */
      if (diff_th != NULL) {

        for (int iesp1 = 0; iesp1 < ngaz; iesp1++) {
          for (int iesp2 = 0; iesp2 < ngaz; iesp2++) {
            coef[iesp1 * (ngaz - 1) + iesp2]
              = 1. +   sqrt(xlabes[iesp1] / xlabes[iesp2])
                     * sqrt(sqrt(roesp[iesp2] / roesp[iesp1]));
            coef[iesp1 * (ngaz - 1) + iesp2]
              *= coef[iesp1 * (ngaz - 1) + iesp2];
            coef[iesp1 * (ngaz - 1) + iesp2]
              /= (sqrt(1. + roesp[iesp1] / roesp[iesp2]) * sqrt(8.));
          }
        }
        /* Lambda */
        diff_th->val[iel] = 0.;

        for (int iesp1 = 0; iesp1 < ngaz; iesp1++) {
          if (yvol[iesp1] > 1e-30) {
            double somphi = 0.;
            for (int iesp2 = 0; iesp2 < ngaz; iesp2++) {
              if (iesp1 != iesp2)
                somphi +=   coef[iesp1 * (ngaz - 1) + iesp2]
                          * yvol[iesp2] / yvol[iesp1];
            }

            diff_th->val[iel] += xlabes[iesp1] / (1. + 1.065 * somphi);
          }
        }

        /* Lambda/Cp */
        if (cs_glob_fluid_properties->icp <= 0)
          diff_th->val[iel] /= cs_glob_fluid_properties->cp0;
        else
          diff_th->val[iel] /= CS_F_(cp)->val[iel];
      }

      /* compute electric conductivity: S/m */
      if (ifcsig >= 0) {
        c_prop->val[iel] = 0.;
        double val = 0.;

        for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
          val += yvol[iesp1] / sigesp[iesp1];

        c_prop->val[iel] = 1. / val;
      }

      /* compute radiative transfer : W/m3 */
      if (cs_glob_elec_option->ixkabe == 1) {
        if (cpro_absco != NULL) { /* May be NULL if no active radiation model */
          double val = 0.;
          for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
            val += yvol[iesp1] * xkabes[iesp1];
          cpro_absco[iel] = val;
        }
      }
      else if (cs_glob_elec_option->ixkabe == 2) {
        CS_F_(radsc)->val[iel] = 0.;
        double val = 0.;

        for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
          val += yvol[iesp1] * xkabes[iesp1];

        CS_F_(radsc)->val[iel] = val;
      }

      /* diffusivity for other properties
       * nothing to do
       * no other properties in this case */

    } /* End of loop on cells */

    BFT_FREE(ym);
    BFT_FREE(yvol);
    BFT_FREE(roesp);
    BFT_FREE(visesp);
    BFT_FREE(cpesp);
    BFT_FREE(sigesp);
    BFT_FREE(xlabes);
    BFT_FREE(xkabes);
    BFT_FREE(coef);
  }

  /* now user properties (for joule effect particulary) */
  cs_user_physical_properties(domain);
}

/*----------------------------------------------------------------------------
 * compute specific electric arc fields
 *----------------------------------------------------------------------------*/

void
cs_elec_compute_fields(const cs_mesh_t  *mesh,
                       int               call_id)
{
  cs_lnum_t  n_cells   = mesh->n_cells;
  cs_lnum_t  n_cells_ext = mesh->n_cells_with_ghosts;
  const int kivisl  = cs_field_key_id("diffusivity_id");

  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];
  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  bool log_active = cs_log_default_is_active();

  /* Reconstructed value */
  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  /* ----------------------------------------------------- */
  /* first call : J, E => J.E                              */
  /* ----------------------------------------------------- */

  if (call_id == 1) {
    /* compute grad(potR) */

    /* Get the calculation option from the field */
    cs_real_3_t *cpro_elefl = (cs_real_3_t *)(CS_F_(elefl)->val);

    cs_field_gradient_scalar(CS_F_(potr),
                             false, /* use_previous_t */
                             1,    /* inc */
                             grad);

    /* compute electric field E = - grad (potR) */
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      cpro_elefl[iel][0] = grad[iel][0];
      cpro_elefl[iel][1] = grad[iel][1];
      cpro_elefl[iel][2] = grad[iel][2];
    }

    /* compute current density j = sig E */
    int diff_id = cs_field_get_key_int(CS_F_(potr), kivisl);
    cs_field_t *c_prop = NULL;
    if (diff_id > -1)
      c_prop = cs_field_by_id(diff_id);

    if (ieljou > 0 || ielarc > 0) {
      cs_real_3_t *cpro_curre = (cs_real_3_t *)(CS_F_(curre)->val);
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        cpro_curre[iel][0] = -c_prop->val[iel] * grad[iel][0];
        cpro_curre[iel][1] = -c_prop->val[iel] * grad[iel][1];
        cpro_curre[iel][2] = -c_prop->val[iel] * grad[iel][2];
      }
    }

    /* compute joule effect : j . E */
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      CS_F_(joulp)->val[iel] =  c_prop->val[iel] *
                               (grad[iel][0] * grad[iel][0] +
                                grad[iel][1] * grad[iel][1] +
                                grad[iel][2] * grad[iel][2]);
    }

    /* compute min max for E and J */
    if (log_active) {
      bft_printf("-----------------------------------------\n"
                 "   Variable         Minimum       Maximum\n"
                 "-----------------------------------------\n");

      /* Grad PotR = -E */
      double vrmin[3], vrmax[3];

      for (int i = 0; i < 3; i++) {
        vrmin[i] = HUGE_VAL;
        vrmax[i] = -HUGE_VAL;
      }

      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (int i = 0; i < 3; i++) {
          vrmin[i] = CS_MIN(vrmin[i], grad[iel][i]);
          vrmax[i] = CS_MAX(vrmax[i], grad[iel][i]);
        }
      }

      cs_parall_min(3, CS_DOUBLE, vrmin);
      cs_parall_max(3, CS_DOUBLE, vrmax);

      for (int i = 0; i < 3; i++) {
        bft_printf("v  Gr_PotR%s    %12.5e  %12.5e\n",
                   cs_glob_field_comp_name_3[i],
                   vrmin[i], vrmax[i]);
      }

      /* current real */
      for (int i = 0; i < 3; i++) {
        vrmin[i] = HUGE_VAL;
        vrmax[i] = -HUGE_VAL;
      }

      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (int i = 0; i < 3; i++) {
          vrmin[i] = CS_MIN(vrmin[i], -c_prop->val[iel] * grad[iel][i]);
          vrmax[i] = CS_MAX(vrmax[i], -c_prop->val[iel] * grad[iel][i]);
        }
      }

      cs_parall_min(3, CS_DOUBLE, vrmin);
      cs_parall_max(3, CS_DOUBLE, vrmax);

      for (int i = 0; i < 3; i++) {
        bft_printf("v  Cour_Re%s    %12.5E  %12.5E\n",
                   cs_glob_field_comp_name_3[i],
                   vrmin[i], vrmax[i]);
      }
      bft_printf("-----------------------------------------\n");
    }

    if (ieljou == 2 || ieljou == 4) {
      /* compute grad(potI) */

      cs_field_gradient_scalar(CS_F_(poti),
                               false, /* use_previous_t */
                               1,    /* inc */
                               grad);

      /* compute electric field E = - grad (potI) */

      /* compute current density j = sig E */
      int diff_id_i = cs_field_get_key_int(CS_F_(poti), kivisl);
      cs_field_t *c_propi = NULL;
      if (diff_id_i > -1)
        c_propi = cs_field_by_id(diff_id_i);

      if (ieljou == 4) {
        cs_real_3_t *cpro_curim = (cs_real_3_t *)(CS_F_(curim)->val);
        for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
          cpro_curim[iel][0] = -c_propi->val[iel] * grad[iel][0];
          cpro_curim[iel][1] = -c_propi->val[iel] * grad[iel][1];
          cpro_curim[iel][2] = -c_propi->val[iel] * grad[iel][2];
        }
      }

      /* compute joule effect : j . E */
      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        CS_F_(joulp)->val[iel] +=   c_propi->val[iel]
                                  * cs_math_3_square_norm(grad[iel]);
      }

      /* compute min max for E and J */
      if (log_active) {

        double vrmin[3], vrmax[3];

        /* Grad PotR = -Ei */

        for (int i = 0; i < 3; i++) {
          vrmin[i] = HUGE_VAL;
          vrmax[i] = -HUGE_VAL;
        }

        for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
          for (int i = 0; i < 3; i++) {
            vrmin[i] = CS_MIN(vrmin[i], grad[iel][0]);
            vrmax[i] = CS_MAX(vrmax[i], grad[iel][0]);
          }
        }

        cs_parall_min(3, CS_DOUBLE, vrmin);
        cs_parall_max(3, CS_DOUBLE, vrmax);

        for (int i = 0; i < 3; i++) {
          bft_printf("v  Gr_PotI%s    %12.5E  %12.5E\n",
                     cs_glob_field_comp_name_3[i],
                     vrmin[i], vrmax[i]);
        }

        /* Imaginary current */

        for (int i = 0; i < 3; i++) {
          vrmin[i] = HUGE_VAL;
          vrmax[i] = -HUGE_VAL;
        }

        for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
          for (int i = 0; i < 3; i++) {
            vrmin[i] = CS_MIN(vrmin[i], -c_propi->val[iel] * grad[iel][i]);
            vrmax[i] = CS_MAX(vrmax[i], -c_propi->val[iel] * grad[iel][i]);
          }
        }

        cs_parall_min(3, CS_DOUBLE, vrmin);
        cs_parall_max(3, CS_DOUBLE, vrmax);

        for (int i = 0; i < 3; i++) {
          bft_printf("v  Cour_Im%s    %12.5E  %12.5E\n",
                     cs_glob_field_comp_name_3[i],
                     vrmin[i], vrmax[i]);
        }
      }
    }
  }

  /* ----------------------------------------------------- */
  /* second call : A, B, JXB                               */
  /* ----------------------------------------------------- */

  else if (call_id == 2) {

    cs_real_3_t *cpro_magfl = (cs_real_3_t *)(CS_F_(magfl)->val);

    if (ielarc == 2) {
      /* compute magnetic field component B */
      cs_field_t  *fp = cs_field_by_name_try("vec_potential");

      cs_real_33_t *gradv = NULL;
      BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

      cs_field_gradient_vector(fp,
                               false, /* use_previous_t */
                               1,    /* inc */
                               gradv);

      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        cpro_magfl[iel][0] = -gradv[iel][1][2]+gradv[iel][2][1];
        cpro_magfl[iel][1] =  gradv[iel][0][2]-gradv[iel][2][0];
        cpro_magfl[iel][2] = -gradv[iel][0][1]+gradv[iel][1][0];
      }

      BFT_FREE(gradv);
    }
    else if (ielarc == 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Error electric arc with ampere theorem not available\n"));

    /* compute laplace effect j x B */
    cs_real_3_t *cpro_laplf = (cs_real_3_t *)(CS_F_(laplf)->val);
    cs_real_3_t *cpro_curre = (cs_real_3_t *)(CS_F_(curre)->val);
    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      cpro_laplf[iel][0] =    cpro_curre[iel][1] * cpro_magfl[iel][2]
                            - cpro_curre[iel][2] * cpro_magfl[iel][1];
      cpro_laplf[iel][1] =    cpro_curre[iel][2] * cpro_magfl[iel][0]
                            - cpro_curre[iel][0] * cpro_magfl[iel][2];
      cpro_laplf[iel][2] =    cpro_curre[iel][0] * cpro_magfl[iel][1]
                            - cpro_curre[iel][1] * cpro_magfl[iel][0];
    }

    /* compute min max for B */
    if (ielarc > 1 && log_active) {
      /* Grad PotR = -E */
      double vrmin[3], vrmax[3];

      for (int i = 0; i < 3; i++) {
        vrmin[i] = HUGE_VAL;
        vrmax[i] = -HUGE_VAL;
      }

      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        for (int i = 0; i < 3; i++) {
          vrmin[i] = CS_MIN(vrmin[i], cpro_magfl[iel][i]);
          vrmax[i] = CS_MAX(vrmax[i], cpro_magfl[iel][i]);
        }
      }

      cs_parall_min(3, CS_DOUBLE, vrmin);
      cs_parall_max(3, CS_DOUBLE, vrmax);

      for (int i = 0; i < 3; i++) {
        bft_printf("v  Magnetic_field%s    %12.5E  %12.5E\n",
                   cs_glob_field_comp_name_3[i], vrmin[i], vrmax[i]);
      }
    }
  }

  /* Free memory */
  BFT_FREE(grad);
}

/*----------------------------------------------------------------------------
 * compute source terms for energy
 *----------------------------------------------------------------------------*/

void
cs_elec_source_terms(const cs_mesh_t             *mesh,
                     const cs_mesh_quantities_t  *mesh_quantities,
                     int                          f_id,
                     cs_real_t                   *smbrs)
{
  const cs_field_t  *f    = cs_field_by_id(f_id);
  const char        *name = f->name;
  cs_lnum_t          n_cells     = mesh->n_cells;
  cs_lnum_t          n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_real_t   *volume = mesh_quantities->cell_vol;

  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  cs_real_t *w1;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);

  /* enthalpy source term */
  if (strcmp(name, "enthalpy") == 0) {
    if (var_cal_opt.verbosity > 0)
      bft_printf("compute source terms for variable : %s\n", name);

    if (cs_glob_time_step->nt_cur > 2) {
      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        w1[iel] = CS_F_(joulp)->val[iel] * volume[iel];

      if (ielarc >= 1)
        if (cs_glob_elec_option->ixkabe == 2)
          for (cs_lnum_t iel = 0; iel < n_cells; iel++)
            w1[iel] -= CS_F_(radsc)->val[iel] * volume[iel];

      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        smbrs[iel] += w1[iel];

      if (var_cal_opt.verbosity > 0) {
        double valmin = w1[0];
        double valmax = w1[0];

        for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
          valmin = CS_MIN(valmin, w1[iel]);
          valmax = CS_MAX(valmax, w1[iel]);
        }

        cs_parall_min(1, CS_DOUBLE, &valmin);
        cs_parall_max(1, CS_DOUBLE, &valmax);

        bft_printf(" source terms for H min= %14.5E, max= %14.5E\n",
                   valmin, valmax);
      }
    }
  }

  BFT_FREE(w1);
}

/*----------------------------------------------------------------------------
 * compute source terms for vector potential
 *----------------------------------------------------------------------------*/

void
cs_elec_source_terms_v(const cs_mesh_t             *mesh,
                       const cs_mesh_quantities_t  *mesh_quantities,
                       int                          f_id,
                       cs_real_3_t                 *smbrv)
{
  const cs_field_t  *f    = cs_field_by_id(f_id);
  cs_lnum_t          n_cells     = mesh->n_cells;
  const cs_real_t   *volume = mesh_quantities->cell_vol;

  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  /* source term for potential vector */

  if (ielarc >= 2 && f_id == (CS_F_(potva)->id)) {
    cs_real_3_t *cpro_curre = (cs_real_3_t *)(CS_F_(curre)->val);

    if (var_cal_opt.verbosity > 0)
      bft_printf("compute source terms for variable: %s\n", f->name);

    for (cs_lnum_t iel = 0; iel < n_cells; iel++)
      for (int isou = 0; isou < 3; isou++)
        smbrv[iel][isou] += cs_elec_permvi * cpro_curre[iel][isou] * volume[iel];
  }
}

/*----------------------------------------------------------------------------
 * add variables fields
 *----------------------------------------------------------------------------*/

void
cs_elec_add_variable_fields(void)
{
  cs_field_t *f;
  int dim1 = 1;
  int dim3 = 3;

  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");
  const int kivisl = cs_field_key_id("diffusivity_id");

  const cs_data_elec_t  *e_props = cs_glob_elec_properties; /* local name */

  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];
  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];

  {
    int f_id = cs_variable_field_create("enthalpy", "Enthalpy",
                                        CS_MESH_LOCATION_CELLS, dim1);
    f = cs_field_by_id(f_id);
    cs_field_set_key_double(f, kscmin, -cs_math_big_r);
    cs_field_set_key_int(f, kivisl, 0);
    cs_add_model_thermal_field_indexes(f->id);

    /* set thermal model */
    cs_thermal_model_t *thermal = cs_get_glob_thermal_model();
    thermal->thermal_variable = CS_THERMAL_MODEL_ENTHALPY;
  }

  {
    int f_id = cs_variable_field_create("elec_pot_r", "POT_EL_R",
                                        CS_MESH_LOCATION_CELLS, dim1);
    f = cs_field_by_id(f_id);
    cs_field_set_key_double(f, kscmin, -cs_math_big_r);
    cs_field_set_key_double(f, kscmax,  cs_math_big_r);
    cs_field_set_key_int(f, kivisl, 0);
    cs_add_model_field_indexes(f->id);
  }

  if (ieljou == 2 || ieljou == 4) {
    int f_id = cs_variable_field_create("elec_pot_i", "POT_EL_I",
                                        CS_MESH_LOCATION_CELLS, dim1);
    f = cs_field_by_id(f_id);
    cs_field_set_key_double(f, kscmin, -cs_math_big_r);
    cs_field_set_key_double(f, kscmax,  cs_math_big_r);
    cs_field_set_key_int(f, kivisl, 0);
    cs_add_model_field_indexes(f->id);
  }

  if (ielarc > 1) {
    int f_id = cs_variable_field_create("vec_potential", "POT_VEC",
                                        CS_MESH_LOCATION_CELLS, dim3);
    f = cs_field_by_id(f_id);
    //cs_field_set_key_double(f, kscmin, -cs_math_big_r);
    //cs_field_set_key_double(f, kscmax,  cs_math_big_r);
    cs_field_set_key_int(f, kivisl, -1);
    cs_add_model_field_indexes(f->id);
  }

  if (e_props->ngaz > 1) {
    for (int igaz = 0; igaz < e_props->ngaz - 1; igaz++) {
      char *name = NULL;
      char *label = NULL;
      char *suf = NULL;
      BFT_MALLOC(name, strlen("esl_fraction_") + 2 + 1, char);
      BFT_MALLOC(label, strlen("YM_ESL") + 2 + 1, char);
      BFT_MALLOC(suf, 3, char);
      strcpy(name, "esl_fraction_");
      strcpy(label, "YM_ESL");
      sprintf(suf, "%02d", igaz + 1);
      strcat(name, suf);
      strcat(label, suf);

      int f_id = cs_variable_field_create(name, label,
                                          CS_MESH_LOCATION_CELLS, dim1);
      f = cs_field_by_id(f_id);

      cs_field_set_key_double(f, kscmin, 0.);
      cs_field_set_key_double(f, kscmax, 1.);
      cs_field_set_key_int(f, kivisl, 0);
      cs_add_model_field_indexes(f->id);
      BFT_FREE(name);
      BFT_FREE(label);
      BFT_FREE(suf);
    }
  }

  _field_pointer_map_electric_arcs(e_props->ngaz);
}

/*----------------------------------------------------------------------------
 * add properties fields
 *----------------------------------------------------------------------------*/

void
cs_elec_add_property_fields(void)
{
  cs_field_t *f;
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  bool has_previous = false;
  const int klbl   = cs_field_key_id("label");
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");
  const int key_restart_id = cs_field_key_id("restart_file");
  const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];

  {
    f = cs_field_create("temperature",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1, /* dim */
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Temperature");
  }

  {
    f = cs_field_create("joule_power",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1, /* dim */
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "PowJoul");
    cs_field_set_key_int(f, key_restart_id, (int)CS_RESTART_AUXILIARY);
  }

  {
    f = cs_field_create("current_re",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3, /* dim */
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Current_Real");
  }

  {
    f = cs_field_create("electric_field",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,    /* dim */
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Elec_Field");
  }

  /* specific for joule effect */
  if (ieljou == 2 || ieljou == 4) {
    f = cs_field_create("current_im",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3, /* dim */
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Current_Imag");
  }

  /* specific for electric arcs */
  {
    f = cs_field_create("laplace_force",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,    /* dim */
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "For_Lap");

    if (cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > 0)
      cs_field_set_key_int(f, key_restart_id, (int)CS_RESTART_AUXILIARY);

    f = cs_field_create("magnetic_field",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,    /* dim */
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Mag_Field");
  }

  if (cs_glob_elec_option->ixkabe == 1) {
    f = cs_field_create("absorption_coeff",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1, /* dim */
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Coef_Abso");
  }
  else if (cs_glob_elec_option->ixkabe == 2) {
    f = cs_field_create("radiation_source",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1, /* dim */
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "ST_radia");
  }


  _field_pointer_properties_map_electric_arcs();
}

/*----------------------------------------------------------------------------
 * initialize electric fields
 *----------------------------------------------------------------------------*/

void
cs_elec_fields_initialize(const cs_mesh_t   *mesh,
                          int                isuite)
{
  BFT_MALLOC(_elec_option.izreca, mesh->n_i_faces, int);
  for (cs_lnum_t i = 0; i < mesh->n_i_faces; i++)
    _elec_option.izreca[i] = 0;

  cs_lnum_t  n_cells = mesh->n_cells;

  static int ipass = 0;
  ipass += 1;

  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  if (isuite == 0 && ipass == 1) {
    /* enthalpy */
    cs_real_t hinit = 0.;
    if (ielarc > 0) {
      cs_real_t *ym;
      BFT_MALLOC(ym, cs_glob_elec_properties->ngaz, cs_real_t);
      ym[0] = 1.;
      if (cs_glob_elec_properties->ngaz > 1)
        for (int i = 1; i < cs_glob_elec_properties->ngaz; i++)
          ym[i] = 0.;

      cs_real_t tinit = cs_glob_fluid_properties->t0;
      hinit = cs_elec_convert_t_to_h(ym, tinit);
      BFT_FREE(ym);
    }

    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      CS_F_(h)->val[iel] = hinit;
    }

    /* mass fraction of first gas */
    if (cs_glob_elec_properties->ngaz > 1) {
      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        CS_FI_(ycoel, 0)->val[iel] = 1.;
    }
  }
}

/*----------------------------------------------------------------------------
 * scaling electric quantities
 *----------------------------------------------------------------------------*/

void
cs_elec_scaling_function(const cs_mesh_t             *mesh,
                         const cs_mesh_quantities_t  *mesh_quantities,
                         cs_real_t                   *dt)
{
  cs_real_t *volume = mesh_quantities->cell_vol;
  cs_real_t *surfac = mesh_quantities->i_face_normal;
  cs_lnum_t  n_cells   = mesh->n_cells;
  cs_lnum_t  nfac   = mesh->n_i_faces;

  double coepot = 0.;
  double coepoa = 1.;

  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];
  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  if (ielarc >= 1) {
    if (cs_glob_elec_option->modrec == 1) {
      /* standard model */
      double somje = 0.;
      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        somje += CS_F_(joulp)->val[iel] * volume[iel];

      cs_parall_sum(1, CS_DOUBLE, &somje);

      coepot = cs_glob_elec_option->couimp * cs_glob_elec_option->pot_diff
              / CS_MAX(somje, cs_math_epzero);
      coepoa = coepot;

      if (coepot > 1.5)
        coepot = 1.5;
      if (coepot < 0.75)
        coepot = 0.75;

      bft_printf("imposed current / current %14.5e, scaling coef. %14.5e\n",
                 coepoa, coepot);
    }
    else if (cs_glob_elec_option->modrec == 2) {
      /* restrike model */
      cs_gui_elec_model_rec();
      double elcou = 0.;
      cs_real_3_t *cpro_curre = (cs_real_3_t *)(CS_F_(curre)->val);
      if (mesh->halo != NULL)
        cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD,
                                 (cs_real_t *)cpro_curre, 3);
      for (cs_lnum_t ifac = 0; ifac < nfac; ifac++) {
        if (cs_glob_elec_option->izreca[ifac] > 0) {
          bool ok = true;
          for (int idir = 0; idir < 3; idir++)
            if (   fabs(surfac[3 * ifac + idir]) > 0.
                && idir != (cs_glob_elec_option->idreca - 1))
              ok = false;
          if (ok) {
            cs_lnum_t iel = mesh->i_face_cells[ifac][0];
            if (iel < mesh->n_cells)
              elcou += cpro_curre[iel][cs_glob_elec_option->idreca - 1]
                     * surfac[3 * ifac + cs_glob_elec_option->idreca - 1];
          }
        }
      }
      cs_parall_sum(1, CS_DOUBLE, &elcou);
      if (fabs(elcou) > 1.e-6)
        elcou = fabs(elcou);
      else
        elcou = 0.;

      if (fabs(elcou) > 1.e-20)
        coepoa = cs_glob_elec_option->couimp / elcou;

      bft_printf("ELCOU %15.8E\n", elcou);
      _elec_option.elcou = elcou;
    }

    if (   cs_glob_elec_option->modrec == 1
        || cs_glob_elec_option->modrec == 2) {
      double dtj = 1.e15;
      double dtjm = dtj;
      double delhsh = 0.;
      double cdtj = 20.;

      for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
        if (CS_F_(rho)->val[iel] > 0.)
          delhsh = CS_F_(joulp)->val[iel] * dt[iel]
                 / CS_F_(rho)->val[iel];

        if (fabs(delhsh) > 1.e-20)
          dtjm = CS_F_(h)->val[iel] / delhsh;
        else
          dtjm = dtj;
        dtjm = fabs(dtjm);
        dtj = CS_MIN(dtj, dtjm);
      }
      cs_parall_min(1, CS_DOUBLE, &dtj);
      bft_printf("DTJ %15.8E\n", dtj);

      double cpmx = pow(cdtj * dtj, 0.5);
      coepot = cpmx;

      if (cs_glob_time_step->nt_cur > 2) {
        if (coepoa > 1.05)
          coepot = cpmx;
        else
          coepot = coepoa;
      }

      bft_printf(" Cpmx   = %14.5E\n", cpmx);
      bft_printf(" COEPOA   = %14.5E\n", coepoa);
      bft_printf(" COEPOT   = %14.5E\n", coepot);
      bft_printf(" Dpot recale   = %14.5E\n", _elec_option.pot_diff * coepot);

      /* scaling electric fields */
      _elec_option.pot_diff *= coepot;

      /* electric potential (for post treatment) */
      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        CS_F_(potr)->val[iel] *= coepot;

      /* current density */
      if (ielarc > 0) {
        cs_real_3_t *cpro_curre = (cs_real_3_t *)(CS_F_(curre)->val);
        for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
          for (cs_lnum_t i = 0; i < 3; i++)
            cpro_curre[iel][i] *= coepot;
        }
      }

      /* joule effect */
      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        CS_F_(joulp)->val[iel] *= coepot * coepot;
    }
  }

  /* joule effect */
  if (ieljou > 0) {
    /* standard model */
    double somje = 0.;
    for (cs_lnum_t iel = 0; iel < n_cells; iel++)
      somje += CS_F_(joulp)->val[iel] * volume[iel];

    cs_parall_sum(1, CS_DOUBLE, &somje);

    coepot = cs_glob_elec_option->puisim / CS_MAX(somje, cs_math_epzero);
    double coefav = coepot;

    if (coepot > 1.5)
      coepot = 1.5;
    if (coepot < 0.75)
      coepot = 0.75;

    bft_printf("imposed power / sum(jE) %14.5E, scaling coef. %14.5E\n",
               coefav, coepot);

    /* scaling electric fields */
    _elec_option.pot_diff *= coepot;
    _elec_option.coejou *= coepot;

    /* electric potential (for post treatment) */
    if (ieljou != 3 && ieljou != 4)
      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        CS_F_(potr)->val[iel] *= coepot;

    /* current density */
    if (ieljou == 2)
      for (int i = 0; i < 3; i++)
        for (cs_lnum_t iel = 0; iel < n_cells; iel++)
          CS_F_(poti)->val[iel] *= coepot;

    /* joule effect */
    for (cs_lnum_t iel = 0; iel < n_cells; iel++)
      CS_F_(joulp)->val[iel] *= coepot * coepot;
  }

  cs_user_scaling_elec(mesh, mesh_quantities, dt);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert enthalpy to temperature at all boundary faces.
 *
 * This handles both user and model enthalpy conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   h   enthalpy values
 * \param[out]  t   temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_elec_convert_h_to_t_faces(const cs_real_t  h[],
                             cs_real_t        t[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_data_elec_t  *el_p = cs_glob_elec_properties;
  const int n_gasses = el_p->ngaz;

  if (n_gasses == 1) {

    cs_real_t ym[1] = {1.};

    for (cs_lnum_t i = 0; i < n_b_faces; i++)
      t[i] = cs_elec_convert_h_to_t(ym, h[i]);

  }
  else {

    const cs_lnum_t *b_face_cells = m->b_face_cells;

    cs_real_t *ym;
    BFT_MALLOC(ym, n_gasses, cs_real_t);

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

      cs_lnum_t c_id = b_face_cells[f_id];

      ym[n_gasses - 1] = 1.;
      for (int gas_id = 0; gas_id < n_gasses - 1; gas_id++) {
        ym[gas_id] = CS_FI_(ycoel, gas_id)->val[c_id];
        ym[n_gasses - 1] -= ym[gas_id];
      }

      t[f_id] = cs_elec_convert_h_to_t(ym, h[f_id]);

    }

    BFT_FREE(ym);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert single enthalpy value to temperature.
 *
 * \param[in]       ym      mass fraction for each gas
 * \param[in, out]  enthal  enthlapy value
 *
 * \return  temperature value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_elec_convert_h_to_t(const cs_real_t  ym[],
                       cs_real_t        enthal)
{
  int ngaz = cs_glob_elec_properties->ngaz;
  int it   = cs_glob_elec_properties->npoint;

  cs_real_t eh1 = 0.;

  for (int iesp = 0; iesp < ngaz; iesp++)
    eh1 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it-1) + it - 1];

  if (enthal >= eh1) {
    return cs_glob_elec_properties->th[it - 1];
  }

  eh1 = 0.;

  for (int iesp = 0; iesp < ngaz; iesp++)
    eh1 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it-1) + 0];

  if (enthal <= eh1) {
    return cs_glob_elec_properties->th[0];
  }

  for (int itt = 0; itt < cs_glob_elec_properties->npoint - 1; itt++) {
    cs_real_t eh0 = 0.;
    eh1 = 0.;

    for (int iesp = 0; iesp < ngaz; iesp++) {
      eh0 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it-1) + itt];
      eh1 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it-1) + itt+1];
    }

    if (enthal > eh0 && enthal <= eh1) {
      cs_real_t temp = cs_glob_elec_properties->th[itt]
                       + (enthal - eh0) * (  cs_glob_elec_properties->th[itt+1]
                                           - cs_glob_elec_properties->th[itt])
                                        / (eh1 - eh0);
      return temp;
    }
  }

  assert(0);  /* Should not arrive here */
  return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy at all cells.
 *
 * This handles both user and model temperature conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   t   temperature values
 * \param[out]  h   enthalpy values
 */
/*----------------------------------------------------------------------------*/

void
cs_elec_convert_t_to_h_cells(const cs_real_t  t[],
                             cs_real_t        h[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;

  const cs_data_elec_t  *el_p = cs_glob_elec_properties;
  const int n_gasses = el_p->ngaz;

  if (n_gasses == 1) {

    cs_real_t ym[1] = {1.};

    for (cs_lnum_t i = 0; i < n_cells; i++)
      h[i] = cs_elec_convert_t_to_h(ym, t[i]);

  }
  else {

    cs_real_t *ym;
    BFT_MALLOC(ym, n_gasses, cs_real_t);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      ym[n_gasses - 1] = 1.;
      for (int gas_id = 0; gas_id < n_gasses - 1; gas_id++) {
        ym[gas_id] = CS_FI_(ycoel, gas_id)->val[c_id];
        ym[n_gasses - 1] -= ym[gas_id];
      }

      h[c_id] = cs_elec_convert_t_to_h(ym, t[c_id]);

    }

    BFT_FREE(ym);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert temperature to enthalpy at selected boundary faces.
 *
 * This handles both user and model temperature conversions, so can be used
 * safely whenever conversion is needed.
 *
 * \param[in]   n_faces   number of selected faces
 * \param[in]   face_ids  ids of selected faces
 * \param[in]   t         temperature values (defined on all boundary faces)
 * \param[out]  h         enthalpy values (defined on all boundary faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_elec_convert_t_to_h_faces(const cs_lnum_t  n_faces,
                             const cs_lnum_t  face_ids[],
                             const cs_real_t  t[],
                             cs_real_t        h[])
{
  const cs_mesh_t *m = cs_glob_mesh;

  const cs_data_elec_t  *el_p = cs_glob_elec_properties;
  const int n_gasses = el_p->ngaz;

  if (n_gasses == 1) {

    cs_real_t ym[1] = {1.};

    for (cs_lnum_t i = 0; i < n_faces; i++) {
      cs_lnum_t f_id = face_ids[i];
      h[f_id] = cs_elec_convert_t_to_h(ym, t[f_id]);
    }

  }
  else {

    const cs_lnum_t *b_face_cells = m->b_face_cells;

    cs_real_t *ym;
    BFT_MALLOC(ym, n_gasses, cs_real_t);

    for (cs_lnum_t i = 0; i < n_faces; i++) {

      cs_lnum_t f_id = face_ids[i];
      cs_lnum_t c_id = b_face_cells[f_id];
      for (int gas_id = 0; gas_id < n_gasses - 1; gas_id++) {
        ym[gas_id] = CS_FI_(ycoel, gas_id)->val[c_id];
        ym[n_gasses - 1] -= ym[gas_id];
      }

      h[f_id] = cs_elec_convert_t_to_h(ym, t[f_id]);

    }

    BFT_FREE(ym);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert single temperature value to enthalpy.
 *
 * \param[in]       ym    mass fraction for each gas
 * \param[in, out]  temp  temperature value
 *
 * \return  enthalpy values
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_elec_convert_t_to_h(const cs_real_t ym[],
                       cs_real_t       temp)
{
  int ngaz = cs_glob_elec_properties->ngaz;
  int it   = cs_glob_elec_properties->npoint;

  cs_real_t enthal = 0.;

  if (temp >= cs_glob_elec_properties->th[it - 1]) {
    for (int iesp = 0; iesp < ngaz; iesp++)
      enthal += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it-1) + it-1];
  }
  else if (temp <= cs_glob_elec_properties->th[0]) {
    for (int iesp = 0; iesp < ngaz; iesp++)
      enthal += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it-1) + 0];
  }
  else {
    for (int itt = 0; itt < cs_glob_elec_properties->npoint - 1; itt++) {
      if (   temp > cs_glob_elec_properties->th[itt]
          && temp <= cs_glob_elec_properties->th[itt + 1]) {
        cs_real_t eh0 = 0.;
        cs_real_t eh1 = 0.;

        for (int iesp = 0; iesp < ngaz; iesp++) {
          eh0 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it-1) + itt];
          eh1 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it-1) + itt+1];
        }

        enthal = eh0 + (eh1 - eh0) * (temp - cs_glob_elec_properties->th[itt]) /
                       (  cs_glob_elec_properties->th[itt + 1]
                        - cs_glob_elec_properties->th[itt]);

        break;
      }
    }
  }

  return enthal;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create or access function objects specific to
 *        electric arcs models.
 */
/*----------------------------------------------------------------------------*/

void
cs_elec_define_functions(void)
{
  /* For Joule Heating by direct conduction:
     gradient of the imaginary component of the potential */

  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];

  if (ieljou == 2 || ieljou == 4) {
    cs_function_t *f
      = cs_function_define_by_func("elec_pot_gradient_im",
                                   CS_MESH_LOCATION_CELLS,
                                   3,
                                   true,
                                   CS_REAL_TYPE,
                                   _pot_gradient_im_f,
                                   NULL);

    const char label[] = "Pot_Gradient_Im";
    BFT_MALLOC(f->label, strlen(label) + 1, char);
    strcpy(f->label, label);

    f->type = CS_FUNCTION_INTENSIVE;
    f->post_vis = CS_POST_ON_LOCATION;
  }

  /* For Joule heating by direct conduction:
     imaginary component of the current density */

  if (ieljou == 2 || ieljou == 4) {
    cs_function_t *f
      = cs_function_define_by_func("elec_current_im",
                                   CS_MESH_LOCATION_CELLS,
                                   3,
                                   true,
                                   CS_REAL_TYPE,
                                   _current_im_f,
                                   NULL);

    const char label[] = "Current_Im";
    BFT_MALLOC(f->label, strlen(label) + 1, char);
    strcpy(f->label, label);

    f->type = CS_FUNCTION_INTENSIVE;
    f->post_vis = CS_POST_ON_LOCATION;
  }

  /* Calculation of Module of the complex potential */

  if (ieljou == 4) {
    cs_function_t *f
      = cs_function_define_by_func("elec_pot_module",
                                   CS_MESH_LOCATION_CELLS,
                                   1,
                                   true,
                                   CS_REAL_TYPE,
                                   _pot_module_f,
                                   NULL);

    const char label[] = "Pot_Module";
    BFT_MALLOC(f->label, strlen(label) + 1, char);
    strcpy(f->label, label);

    f->type = CS_FUNCTION_INTENSIVE;
    f->post_vis = CS_POST_ON_LOCATION;
  }

  /* Calculation of Argument of the complex potential */

  if (ieljou == 4) {
    cs_function_t *f
      = cs_function_define_by_func("elec_pot_arg",
                                   CS_MESH_LOCATION_CELLS,
                                   1,
                                   true,
                                   CS_REAL_TYPE,
                                   _pot_arg_f,
                                   NULL);

    const char label[] = "Pot_Arg";
    BFT_MALLOC(f->label, strlen(label) + 1, char);
    strcpy(f->label, label);

    f->type = CS_FUNCTION_INTENSIVE;
    f->post_vis = CS_POST_ON_LOCATION;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
