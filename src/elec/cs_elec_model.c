/*============================================================================
 * Base electrical model data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include "cs_log.h"
#include "cs_parall.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_location.h"
#include "cs_time_step.h"
#include "cs_field.h"
#include "cs_parameters.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_field_operator.h"
#include "cs_physical_constants.h"
#include "cs_thermal_model.h"
#include "cs_turbulence_model.h"
#include "cs_gui_specific_physics.h"

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
  \file cs_elec_model.c
        Base electrical model data.
*/

/*----------------------------------------------------------------------------*/

/*! \struct cs_elec_option_t

  \brief option for electric model

  \var  cs_elec_option_t::ieljou
        joule model
        - -1: module not activated
        -  1: use of a real potential
        -  2: use of a complex potential
        -  3: use of real potential and specific boundary conditions
        -  4: use of complex potential and specific boundary conditions
  \var  cs_elec_option_t::ielarc
        electric arc model
        - -1: module not activated
        -  1: determination of the magnetic field by means of the Ampere’ theorem
        -  2: determination of the magnetic field by means of the vector potential
  \var  cs_elec_option_t::ielion
        ionic conduction model
  \var  cs_elec_option_t::ixkabe
        model for radiative properties
        - 0: last column read but not use
        - 1: last column : absorption coefficient
        - 2: last column : radiative TS
  \var  cs_elec_option_t::ntdcla
        first iteration to take into account restrike model
  \var  cs_elec_option_t::irestrike
        indicate if restrike or not
  \var  cs_elec_option_t::restrike_point
        coordinates for restrike point
  \var  cs_elec_option_t::crit_reca
        define plane for scaling
  \var  cs_elec_option_t::ielcor
        indicate if scaling or not
  \var  cs_elec_option_t::modrec
        model for scaling
        - 1: volumic power
        - 2: by plane
        - 3: user function
  \var  cs_elec_option_t::idreca
        direction for current for scaling
  \var  cs_elec_option_t::izreca
        indicator for faces for scaling
  \var  cs_elec_option_t::couimp
        imposed current
  \var  cs_elec_option_t::pot_diff
        potential difference
  \var  cs_elec_option_t::puisim
        imposed power
  \var  cs_elec_option_t::coejou
        coefficient for scaling
  \var  cs_elec_option_t::elcou
        current in scaling plane
  \var  cs_elec_option_t::ficfpp
        data file name
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
        number of gaz in electrical data file
  \var  cs_data_elec_t::npoint
        number of point in electrical data file for each gaz
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

static cs_elec_option_t       *_elec_option     = NULL;
static cs_data_elec_t         *_elec_properties = NULL;
static cs_data_joule_effect_t *_transformer     = NULL;
const cs_elec_option_t        *cs_glob_elec_option;
const cs_data_elec_t          *cs_glob_elec_properties;
const cs_data_joule_effect_t  *cs_glob_transformer;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*!
 * vaccum magnetic permeability constant (H/m). (= 1.2566e-6)
 *
 */
const double cs_elec_permvi = 1.2566e-6;

/*!
 * vaccum permittivity constant (F/m). (= 8.854e-12)
 *
 */
const double cs_elec_epszer = 8.854e-12;

const double epzero = 1.e-12;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

void
cs_f_elec_model_get_pointers(int     **ngazge,
                             int     **ielcor,
                             double  **pot_diff,
                             double  **coejou,
                             double  **elcou,
                             int     **irestrike,
                             int     **ntdcla,
                             double  **restrike_pointX,
                             double  **restrike_pointY,
                             double  **restrike_pointZ);

/*============================================================================
 * Private function definitions
 *============================================================================*/

void
_cs_electrical_model_verify(void)
{
  bool verif = true;

  if (cs_glob_elec_option->ielarc != -1 && cs_glob_elec_option->ielarc !=  2)
    bft_error(__FILE__, __LINE__, 0,
              _("Error for electric arc model\n"
                "only choice -1 or 2 are permitted yet\n"
                "model selected : \"%i\";\n"),
              cs_glob_elec_option->ielarc);

  if (cs_glob_elec_option->ieljou != -1 && cs_glob_elec_option->ieljou !=  1 &&
      cs_glob_elec_option->ieljou !=  2 && cs_glob_elec_option->ieljou !=  3 &&
      cs_glob_elec_option->ieljou !=  4)
    bft_error(__FILE__, __LINE__, 0,
              _("Error for joule model\n"
                "only choice -1, 1, 2, 3 or 4 are permitted yet\n"
                "model selected : \"%i\";\n"),
              cs_glob_elec_option->ieljou);

  if (cs_glob_elec_option->ielion != -1)
    bft_error(__FILE__, __LINE__, 0,
              _("Error for ionic conduction model\n"
                "only choice -1 is permitted yet\n"
                "model selected : \"%i\";\n"),
              cs_glob_elec_option->ielion);

  /* options */
  if (cs_glob_elec_option->ielcor != 0 && cs_glob_elec_option->ielcor != 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Error for scaling model\n"
                "only choice -1 or 2 are permitted yet\n"
                "model selected : \"%i\";\n"),
              cs_glob_elec_option->ielcor);

  if (cs_glob_elec_option->ielcor == 1) {
    if (cs_glob_elec_option->ielarc > 0) {
      if (cs_glob_elec_option->couimp < 0.) {
        bft_printf("value for COUIMP must be strictly positive\n");
        verif = false;
      }
      if (cs_glob_elec_option->pot_diff < 0.) {
        bft_printf("value for DPOT must be strictly positive\n");
        verif = false;
      }
    }
    if (cs_glob_elec_option->ieljou > 0) {
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
  return;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide acces to cs_elec_option
 *----------------------------------------------------------------------------*/

cs_elec_option_t *
cs_get_glob_elec_option(void)
{
  return _elec_option;
}

/*----------------------------------------------------------------------------
 * Provide acces to cs_glob_transformer
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
cs_electrical_model_initialize(cs_int_t ielarc,
                               cs_int_t ieljou,
                               cs_int_t ielion)
{
    BFT_MALLOC(_elec_option, 1, cs_elec_option_t);

    if (ielarc > 0)
      BFT_MALLOC(_elec_properties, 1, cs_data_elec_t);

    if (ieljou >= 3)
      BFT_MALLOC(_transformer, 1, cs_data_joule_effect_t);

    _elec_option->ielarc    = ielarc;
    _elec_option->ieljou    = ieljou;
    _elec_option->ielion    = ielion;
    _elec_option->ixkabe    = 0;
    _elec_option->ntdcla    = 1;
    _elec_option->irestrike = 0;
    for (int i = 0; i < 2; i++)
      _elec_option->restrike_point[i] = 0.;
    _elec_option->izreca    = NULL;
    _elec_option->elcou     = 0.;
    _elec_option->ficfpp    = NULL;
    _elec_option->ielcor    = 0;
    _elec_option->couimp    = 0.;
    _elec_option->puisim    = 0.;
    _elec_option->pot_diff  = 0.;
    _elec_option->coejou    = 1.;
    _elec_option->modrec    = 1;    /* standard model */
    _elec_option->idreca    = 3;
    for (int i = 0; i < 3; i++)
      _elec_option->crit_reca[i] = 0.;
    _elec_option->crit_reca[4] = 0.0002;

    BFT_MALLOC(_elec_option->ficfpp, 7, char);
    strcpy(_elec_option->ficfpp, "dp_ELE");

    cs_glob_elec_option     = _elec_option;
    cs_glob_elec_properties = _elec_properties;
    cs_glob_transformer     = _transformer;

    cs_fluid_properties_t *fluid_properties = cs_get_glob_fluid_properties();
    fluid_properties->icp = 1;
    fluid_properties->irovar = 1;
    fluid_properties->ivivar = 1;

    return;
}

/*----------------------------------------------------------------------------
 * Destroy structures for electrical model
 *
 *----------------------------------------------------------------------------*/

void
cs_electrical_model_finalize(cs_int_t ielarc,
                             cs_int_t ieljou)
{
  if (ielarc > 0) {
    BFT_FREE(_elec_properties->th);
    BFT_FREE(_elec_properties->ehgaz);
    BFT_FREE(_elec_properties->rhoel);
    BFT_FREE(_elec_properties->cpel);
    BFT_FREE(_elec_properties->sigel);
    BFT_FREE(_elec_properties->visel);
    BFT_FREE(_elec_properties->xlabel);
    BFT_FREE(_elec_properties->xkabel);
    BFT_FREE(_elec_properties);
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

  BFT_FREE(_elec_option->ficfpp);
  BFT_FREE(_elec_option->izreca);

  BFT_FREE(_elec_option);
}

/*----------------------------------------------------------------------------
 * Specific initialization for electric arc
 *----------------------------------------------------------------------------*/

void
cs_electrical_model_specific_initialization(      cs_real_t *srrom,
                                                  cs_real_t *visls0,
                                                  cs_real_t *diftl0,
                                                  cs_int_t  *iconv,
                                                  cs_int_t  *istat,
                                                  cs_int_t  *idiff,
                                                  cs_int_t  *idifft,
                                                  cs_int_t  *idircl,
                                                  cs_int_t  *isca,
                                                  cs_real_t *blencv,
                                                  cs_real_t *sigmas,
                                                  cs_int_t  *iwarni,
                                            const cs_int_t   iihmpr)
{
  *srrom = 0.;

  cs_field_t *f = NULL;
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  const int keysca = cs_field_key_id("scalar_id");
  cs_var_cal_opt_t var_cal_opt;

  /* specific initialization for field */
  f = CS_F_(potr);
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
  int id = cs_field_get_key_int(f, keysca) - 1;
  iconv[isca[id] - 1]  = 0;
  istat[isca[id] - 1]  = 0;
  idiff[isca[id] - 1]  = 1;
  idifft[isca[id] - 1] = 0;
  idircl[isca[id] - 1] = 1;

  if (cs_glob_elec_option->ieljou == 2 ||
      cs_glob_elec_option->ieljou == 4) {
    f = CS_F_(poti);
    cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
    id = cs_field_get_key_int(f, keysca) - 1;
    iconv[isca[id] - 1]  = 0;
    istat[isca[id] - 1]  = 0;
    idiff[isca[id] - 1]  = 1;
    idifft[isca[id] - 1] = 0;
    idircl[isca[id] - 1] = 1;
  }

  /* TODO when vector field
  if (cs_glob_elec_option->ielarc > 1) {
    for (int i = 0; i < 3 ; i++) {
      f = CS_FI_(potva, i);
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      var_cal_opt.iconv  = 0;
      var_cal_opt.istat  = 0;
      var_cal_opt.idiff  = 1;
      var_cal_opt.idifft = 0;
      // TODO var_cal_opt.idircl = 1;

      int i = cs_field_get_key_int(fp1, keysca) - 1;
      visls0[i] = 1.;
    }
  }
  */
  if (cs_glob_elec_option->ielarc > 1) {
    cs_field_t  *fp1 = cs_field_by_name_try("vec_potential_01");
    cs_field_t  *fp2 = cs_field_by_name_try("vec_potential_02");
    cs_field_t  *fp3 = cs_field_by_name_try("vec_potential_03");
    cs_field_get_key_struct(fp1, key_cal_opt_id, &var_cal_opt);
    id = cs_field_get_key_int(fp1, keysca) - 1;
    iconv[isca[id] - 1]  = 0;
    istat[isca[id] - 1]  = 0;
    idiff[isca[id] - 1]  = 1;
    idifft[isca[id] - 1] = 0;
    idircl[isca[id] - 1] = 1;
    visls0[id] = 1.;

    cs_field_get_key_struct(fp2, key_cal_opt_id, &var_cal_opt);
    id = cs_field_get_key_int(fp2, keysca) - 1;
    iconv[isca[id] - 1]  = 0;
    istat[isca[id] - 1]  = 0;
    idiff[isca[id] - 1]  = 1;
    idifft[isca[id] - 1] = 0;
    idircl[isca[id] - 1] = 1;
    visls0[id] = 1.;

    cs_field_get_key_struct(fp3, key_cal_opt_id, &var_cal_opt);
    id = cs_field_get_key_int(fp3, keysca) - 1;
    iconv[isca[id] - 1]  = 0;
    istat[isca[id] - 1]  = 0;
    idiff[isca[id] - 1]  = 1;
    idifft[isca[id] - 1] = 0;
    idircl[isca[id] - 1] = 1;
    visls0[id] = 1.;
  }

  /* for all specific field */
  f = CS_F_(h);
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
  id = cs_field_get_key_int(f, keysca) - 1;
  if (iwarni[isca[id] - 1] == -10000)
    iwarni[isca[id] - 1] = 1;
  blencv[isca[id] - 1] = 1.;
  sigmas[id] = 0.7;

  f = CS_F_(potr);
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
  id = cs_field_get_key_int(f, keysca) - 1;
  if (iwarni[isca[id] - 1] == -10000)
    iwarni[isca[id] - 1] = 1;
  blencv[isca[id] - 1] = 1.;
  sigmas[id] = 0.7;

  if (cs_glob_elec_option->ieljou == 2 ||
      cs_glob_elec_option->ieljou == 4) {
    f = CS_F_(poti);
    cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
    id = cs_field_get_key_int(f, keysca) - 1;
    if (iwarni[isca[id] - 1] == -10000)
      iwarni[isca[id] - 1] = 1;
    blencv[isca[id] - 1] = 1.;
    sigmas[id] = 0.7;
  }

  /* TODO when vector field
  if (cs_glob_elec_option->ielarc > 1) {
    for (int i = 0; i < 3; i++) {
      f = CS_FI_(potva, i);
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      if (var_cal_opt.iwarni == -10000)
        var_cal_opt.iwarni = 1;
      var_cal_opt.blencv = 1.;
      var_cal_opt.ischcv = 1;
      var_cal_opt.isstpc = 0;
      var_cal_opt.ircflu = 0;
    }
  }
  */
  if (cs_glob_elec_option->ielarc > 1) {
    cs_field_t  *fp1 = cs_field_by_name_try("vec_potential_01");
    cs_field_t  *fp2 = cs_field_by_name_try("vec_potential_02");
    cs_field_t  *fp3 = cs_field_by_name_try("vec_potential_03");
    cs_field_get_key_struct(fp1, key_cal_opt_id, &var_cal_opt);
    id = cs_field_get_key_int(fp1, keysca) - 1;
    if (iwarni[isca[id] - 1] == -10000)
      iwarni[isca[id] - 1] = 1;
    blencv[isca[id] - 1] = 1.;
    sigmas[id] = 0.7;

    cs_field_get_key_struct(fp2, key_cal_opt_id, &var_cal_opt);
    id = cs_field_get_key_int(fp2, keysca) - 1;
    if (iwarni[isca[id] - 1] == -10000)
      iwarni[isca[id] - 1] = 1;
    blencv[isca[id] - 1] = 1.;
    sigmas[id] = 0.7;

    cs_field_get_key_struct(fp3, key_cal_opt_id, &var_cal_opt);
    id = cs_field_get_key_int(fp3, keysca) - 1;
    if (iwarni[isca[id] - 1] == -10000)
      iwarni[isca[id] - 1] = 1;
    blencv[isca[id] - 1] = 1.;
    sigmas[id] = 0.7;
  }

  if (cs_glob_elec_properties->ngaz > 1) {
    for (int igaz = 0; igaz < cs_glob_elec_properties->ngaz - 1; igaz++) {
      f = CS_FI_(ycoel, igaz);
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
      id = cs_field_get_key_int(f, keysca) - 1;
      if (iwarni[isca[id] - 1] == -10000)
        iwarni[isca[id] - 1] = 1;
      blencv[isca[id] - 1] = 1.;
      sigmas[id] = 0.7;
    }
  }

  if (iihmpr == 1) {
    CS_PROCF(uicpi1,UICPI1) (srrom, diftl0);
    uieli1();
    _elec_option->pot_diff = 1000.;
  }

  _cs_electrical_model_verify();

  return;
}

/*----------------------------------------------------------------------------
 * Read properties file
 *
 *----------------------------------------------------------------------------*/

void
cs_electrical_properties_read(cs_int_t ielarc,
                              cs_int_t ieljou)
{
  if (ielarc <= 0 && ieljou < 3)
    return;

  FILE *file;
  char str[LG_MAX];
  file = fopen(cs_glob_elec_option->ficfpp, "r");

  if (!file)
    bft_error(__FILE__, __LINE__, 0,
              _("Error can not open file \"%s\";\n"),
              cs_glob_elec_option->ficfpp);

  /* Position at the beginning of the file */
  fseek(file, 0, SEEK_SET);

  int nb_line_tot = 0;

  /* read file for electric arc properties */
  if (ielarc > 0) {
    int iesp = 0;
    int it = 0;

    while (fgets(str, LG_MAX, file) != NULL) {
      nb_line_tot++;
      if (nb_line_tot < 8)
        continue;

      /* read number of fluids and number of points */
      if (nb_line_tot == 8)
        sscanf(str, "%d %d",
               &(_elec_properties->ngaz),
               &(_elec_properties->npoint));

      if (_elec_properties->ngaz <= 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("incorrect number of species \"%i\";\n"),
                  _elec_properties->ngaz);

      double size = cs_glob_elec_properties->ngaz * cs_glob_elec_properties->npoint;

      if (nb_line_tot == 8)
      {
        BFT_MALLOC(_elec_properties->th,  cs_glob_elec_properties->npoint, double);
        BFT_MALLOC(_elec_properties->ehgaz,  size, double);
        BFT_MALLOC(_elec_properties->rhoel,  size, double);
        BFT_MALLOC(_elec_properties->cpel,   size, double);
        BFT_MALLOC(_elec_properties->sigel,  size, double);
        BFT_MALLOC(_elec_properties->visel,  size, double);
        BFT_MALLOC(_elec_properties->xlabel, size, double);
        BFT_MALLOC(_elec_properties->xkabel, size, double);
      }

      if (nb_line_tot < 14)
        continue;

      if (nb_line_tot == 14)
        sscanf(str, "%i", &(_elec_option->ixkabe));

      if (cs_glob_elec_option->ixkabe < 0 || cs_glob_elec_option->ixkabe >= 3)
        bft_error(__FILE__, __LINE__, 0,
                  _("incorrect choice for radiative model \"%i\";\n"),
                  cs_glob_elec_option->ixkabe < 0);

      if (nb_line_tot < 22)
        continue;

      if (nb_line_tot >= 22) {
        sscanf(str, "%lf %lf %lf %lf %lf %lf %lf %lf",
               &(_elec_properties->th[it]),
               &(_elec_properties->ehgaz[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties->rhoel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties->cpel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties->sigel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties->visel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties->xlabel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]),
               &(_elec_properties->xkabel[iesp *  (cs_glob_elec_properties->npoint - 1) + it]));
        it++;
        if (it == cs_glob_elec_properties->npoint) {
          iesp++;
          it = 0;
        }
      }
    }
  }

#if 0
  for (int it = 0; it < cs_glob_elec_properties->npoint; it++)
  bft_printf("read ficfpp %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n",
             _elec_properties->th[it],
             _elec_properties->ehgaz[it],
             _elec_properties->rhoel[it],
             _elec_properties->cpel[it],
             _elec_properties->sigel[it],
             _elec_properties->visel[it],
             _elec_properties->xlabel[it],
             _elec_properties->xkabel[it]);
#endif

  /* read file for joule effect */
  if (ieljou >= 3) {
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

        BFT_MALLOC(_transformer->tenspr,  cs_glob_transformer->nbtrf, double);
        BFT_MALLOC(_transformer->rnbs,    cs_glob_transformer->nbtrf, double);
        BFT_MALLOC(_transformer->zr,      cs_glob_transformer->nbtrf, double);
        BFT_MALLOC(_transformer->zi,      cs_glob_transformer->nbtrf, double);
        BFT_MALLOC(_transformer->ibrpr,   cs_glob_transformer->nbtrf, int);
        BFT_MALLOC(_transformer->ibrsec,  cs_glob_transformer->nbtrf, int);

        // alloc for boundary conditions
        BFT_MALLOC(_transformer->uroff,  cs_glob_transformer->nbtrf, double);
        BFT_MALLOC(_transformer->uioff,  cs_glob_transformer->nbtrf, double);
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

      if (nb_line_tot == 7 + cs_glob_transformer->nbtrf * 6)
      {
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
  }

  fclose(file);
}

/*----------------------------------------------------------------------------
 * convert enthalpy-temperature
 *----------------------------------------------------------------------------*/

void
cs_elec_convert_h_t(cs_int_t   mode,
                    cs_real_t *ym,
                    cs_real_t *enthal,
                    cs_real_t *temp)
{
  int ngaz = cs_glob_elec_properties->ngaz;
  int it   = cs_glob_elec_properties->npoint;

  /* convert temperature to enthalpy */
  if (mode == -1) {
    *enthal = 0.;

    if (*temp >= cs_glob_elec_properties->th[it - 1]) {
      for (int iesp = 0; iesp < ngaz; iesp++)
        *enthal += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it - 1) + it - 1];
    }
    else if (*temp <= cs_glob_elec_properties->th[0]) {
      for (int iesp = 0; iesp < ngaz; iesp++)
        *enthal += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it - 1) + 0];
    }
    else {
      for (int itt = 0; itt < cs_glob_elec_properties->npoint - 1; itt++)
      {
        if (*temp > cs_glob_elec_properties->th[itt] &&
            *temp <= cs_glob_elec_properties->th[itt + 1]) {
          double eh0 = 0.;
          double eh1 = 0.;

          for (int iesp = 0; iesp < ngaz; iesp++)
          {
            eh0 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it - 1) + itt    ];
            eh1 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it - 1) + itt + 1];
          }

          *enthal = eh0 + (eh1 - eh0) * (*temp - cs_glob_elec_properties->th[itt]) /
                          (cs_glob_elec_properties->th[itt + 1] - cs_glob_elec_properties->th[itt]);

          break;
        }
      }
    }
    return;
  }
  /* convert enthalpy to temperature */
  else if (mode == 1) {
    double eh1 = 0.;

    for (int iesp = 0; iesp < ngaz; iesp++)
      eh1 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it - 1) + it - 1];

    if (*enthal >= eh1) {
      *temp = cs_glob_elec_properties->th[it - 1];
      return;
    }

    eh1 = 0.;

    for (int iesp = 0; iesp < ngaz; iesp++)
      eh1 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it - 1) + 0];

    if (*enthal <= eh1) {
      *temp = cs_glob_elec_properties->th[0];
      return;
    }

    for (int itt = 0; itt < cs_glob_elec_properties->npoint - 1; itt++)
    {
      double eh0 = 0.;
      eh1 = 0.;

      for (int iesp = 0; iesp < ngaz; iesp++)
      {
        eh0 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it - 1) + itt    ];
        eh1 += ym[iesp] * cs_glob_elec_properties->ehgaz[iesp * (it - 1) + itt + 1];
      }

      if (*enthal > eh0 && *enthal <= eh1)
      {
        *temp = cs_glob_elec_properties->th[itt] +
                  (*enthal - eh0) *
                  (cs_glob_elec_properties->th[itt + 1] -cs_glob_elec_properties->th[itt]) /
                  (eh1 - eh0);
        break;
      }
    }
    return;
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("electric module : \n"
                "bad value for mode (integer equal to -1 or 1 : %i here.\n"),
              mode);
}

/*----------------------------------------------------------------------------
 * compute physical properties
 *----------------------------------------------------------------------------*/

void
cs_elec_physical_properties(const cs_mesh_t *mesh,
                            const cs_mesh_quantities_t *mesh_quantities,
                            cs_real_t srrom)
{
  static long ipass = 0;
  int nt_cur = cs_glob_time_step->nt_cur;
  int isrrom = 0;
  cs_field_t *f = NULL;
  cs_lnum_t  ncel = mesh->n_cells;
  const int keysca = cs_field_key_id("scalar_diffusivity_id");
  int diff_id = cs_field_get_key_int(CS_F_(potr), keysca);
  cs_field_t *c_prop = NULL;
  if (diff_id > -1)
    c_prop = cs_field_by_id(diff_id);
  ipass++;

  if (nt_cur > 1 && srrom > 0.)
    isrrom = 1;

  /* joule effect                  */
  /* law must be specified by user */
  int ifcvsl = cs_field_get_key_int(CS_F_(h), keysca);
  cs_field_t *diff_th;
  if (ifcvsl > 0)
    diff_th = cs_field_by_id(ifcvsl);

  /* electric arc                  */
  if (cs_glob_elec_option->ielarc > 0)
  {
    if (ipass == 1)
      bft_printf("electric arc module : properties read on file\n");

    /* compute temperature from enthalpy */
    int mode = 1;
    int ngaz = cs_glob_elec_properties->ngaz;
    int npt  = cs_glob_elec_properties->npoint;

    double *ym, *yvol, *roesp, *visesp, *cpesp, *sigesp, *xlabes, *xkabes, *coef;
    BFT_MALLOC(ym,     ngaz, double);
    BFT_MALLOC(yvol,   ngaz, double);
    BFT_MALLOC(roesp,  ngaz, double);
    BFT_MALLOC(visesp, ngaz, double);
    BFT_MALLOC(cpesp,  ngaz, double);
    BFT_MALLOC(sigesp, ngaz, double);
    BFT_MALLOC(xlabes, ngaz, double);
    BFT_MALLOC(xkabes, ngaz, double);
    BFT_MALLOC(coef,   ngaz * ngaz, double);

    int ifcsig = cs_field_get_key_int(CS_F_(potr), keysca);

    if (ngaz == 1)
    {
      ym[0] = 1.;

      for (int iel = 0; iel < ncel; iel++)
        cs_elec_convert_h_t(mode, ym,
                          &(CS_F_(h)->val[iel]),
                          &(CS_F_(t)->val[iel]));
    }
    else {

      for (int iel = 0; iel < ncel; iel++)
      {
        ym[ngaz - 1] = 1.;

        for (int ii = 0; ii < ngaz - 1; ii++)
        {
          ym[ii] = CS_FI_(ycoel, ii)->val[iel];
          ym[ngaz - 1] -= ym[ii];
        }

        cs_elec_convert_h_t(mode, ym,
                          &(CS_F_(h)->val[iel]),
                          &(CS_F_(t)->val[iel]));
      }
    }

    /* interpolate properties */
    for (int iel = 0; iel < ncel; iel++)
    {
      // temperature
      double tp = CS_F_(t)->val[iel];

      // determine point
      int it = -1;
      if (tp <= cs_glob_elec_properties->th[0])
        it = 0;
      else if (tp >= cs_glob_elec_properties->th[npt - 1])
        it = npt - 1;
      else
      {
        for (int iiii = 0; iiii < npt - 1; iiii++)
          if (tp > cs_glob_elec_properties->th[iiii] &&
              tp <= cs_glob_elec_properties->th[iiii + 1]) {
            it = iiii;
            break;
          }
      }
      if (it == -1)
        bft_error(__FILE__, __LINE__, 0,
                  _("electric module : properties read on file\n"
                    "Warning : error in cs_elec_physical_properties\n"
                    "Invalid reading with temperature : %f.\n"),
                  tp);

      /* mass frction */
      ym[ngaz - 1] = 1.;

      for (int ii = 0; ii < ngaz - 1; ii++)
      {
        ym[ii] = CS_FI_(ycoel, ii)->val[iel];
        ym[ngaz - 1] -= ym[ii];
      }

      /* density, viscosity, ... for each species */
      if (it == 0)
      {
        for (int ii = 0; ii < ngaz; ii++)
        {
          roesp[ii]  = cs_glob_elec_properties->rhoel[ii * (npt - 1)];
          visesp[ii] = cs_glob_elec_properties->visel[ii * (npt - 1)];
          cpesp[ii]  = cs_glob_elec_properties->cpel[ii * (npt - 1)];
          sigesp[ii] = cs_glob_elec_properties->sigel[ii * (npt - 1)];
          xlabes[ii] = cs_glob_elec_properties->xlabel[ii * (npt - 1)];

          if (cs_glob_elec_option->ixkabe > 0)
            xkabes[ii] = cs_glob_elec_properties->xkabel[ii * (npt - 1)];
        }
      }
      else if (it == npt - 1)
      {
        bft_printf("valeur constante = derniere valeur de la table\n");
        for (int ii = 0; ii < ngaz; ii++)
        {
          roesp[ii]  = cs_glob_elec_properties->rhoel[ii * (npt - 1) + npt - 1];
          visesp[ii] = cs_glob_elec_properties->visel[ii * (npt - 1) + npt - 1];
          cpesp[ii]  = cs_glob_elec_properties->cpel[ii * (npt - 1) + npt - 1];
          sigesp[ii] = cs_glob_elec_properties->sigel[ii * (npt - 1) + npt - 1];
          xlabes[ii] = cs_glob_elec_properties->xlabel[ii * (npt - 1) + npt - 1];

          if (cs_glob_elec_option->ixkabe > 0)
            xkabes[ii] = cs_glob_elec_properties->xkabel[ii * (npt - 1) + npt - 1];
        }
      }
      else
      {
        double delt = cs_glob_elec_properties->th[it + 1] -
                      cs_glob_elec_properties->th[it];

        for (int ii = 0; ii < ngaz; ii++)
        {
          double alpro = (cs_glob_elec_properties->rhoel[ii * (npt - 1) + it + 1] -
                         cs_glob_elec_properties->rhoel[ii * (npt - 1) + it]) / delt;
          roesp[ii]  = cs_glob_elec_properties->rhoel[ii * (npt - 1) + it] +
                       alpro * (tp -cs_glob_elec_properties->th[it]);

          double alpvis = (cs_glob_elec_properties->visel[ii * (npt - 1) + it + 1] -
                          cs_glob_elec_properties->visel[ii * (npt - 1) + it]) / delt;
          visesp[ii] = cs_glob_elec_properties->visel[ii * (npt - 1) + it] +
                       alpvis * (tp -cs_glob_elec_properties->th[it]);

          double alpcp = (cs_glob_elec_properties->cpel[ii * (npt - 1) + it + 1] -
                         cs_glob_elec_properties->cpel[ii * (npt - 1) + it]) / delt;
          cpesp[ii]  = cs_glob_elec_properties->cpel[ii * (npt - 1) + it] +
                       alpcp * (tp -cs_glob_elec_properties->th[it]);

          double alpsig = (cs_glob_elec_properties->sigel[ii * (npt - 1) + it + 1] -
                          cs_glob_elec_properties->sigel[ii * (npt - 1) + it]) / delt;
          sigesp[ii] = cs_glob_elec_properties->sigel[ii * (npt - 1) + it] +
                       alpsig * (tp -cs_glob_elec_properties->th[it]);

          double alplab = (cs_glob_elec_properties->xlabel[ii * (npt - 1) + it + 1] -
                          cs_glob_elec_properties->xlabel[ii * (npt - 1) + it]) / delt;
          xlabes[ii] = cs_glob_elec_properties->xlabel[ii * (npt - 1) + it] +
                       alplab * (tp -cs_glob_elec_properties->th[it]);

          if (cs_glob_elec_option->ixkabe > 0) {
            double alpkab = (cs_glob_elec_properties->xkabel[ii * (npt - 1) + it + 1] -
                            cs_glob_elec_properties->xkabel[ii * (npt - 1) + it]) / delt;
            xkabes[ii] = cs_glob_elec_properties->xkabel[ii * (npt - 1) + it] +
                         alpkab * (tp -cs_glob_elec_properties->th[it]);
          }
        }
      }

      /* compute density */
      double rhonp1 = 0.;

      for (int ii = 0; ii < ngaz; ii++)
        rhonp1 += ym[ii] / roesp[ii];

      rhonp1 = 1. / rhonp1;

      if (isrrom == 1)
        CS_F_(rho)->val[iel] = CS_F_(rho)->val[iel] * srrom +
                              (1. - srrom) * rhonp1;
      else
        CS_F_(rho)->val[iel] = rhonp1;

      for (int ii = 0; ii < ngaz; ii++)
      {
        yvol[ii] = ym[ii] * roesp[ii] / CS_F_(rho)->val[iel];
        if (yvol[ii] <= 0.)
          yvol[ii] = epzero * epzero;
      }

      /* compute molecular viscosity : kg/(m s) */
      for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
        for (int iesp2 = 0; iesp2 < ngaz; iesp2++)
        {
          coef[iesp1 * (ngaz - 1) + iesp2] = 1. + sqrt(visesp[iesp1] / visesp[iesp2]) *
                                             sqrt(sqrt(roesp[iesp2] / roesp[iesp1]));
          coef[iesp1 * (ngaz - 1) + iesp2] *= coef[iesp1 * (ngaz - 1) + iesp2];
          coef[iesp1 * (ngaz - 1) + iesp2] /= (sqrt(1. + roesp[iesp1] / roesp[iesp2]) * sqrt(8.));
        }

      CS_F_(mu)->val[iel] = 0.;

      for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
      {
        double somphi = 0.;
        for (int iesp2 = 0; iesp2 < ngaz; iesp2++)
          if (iesp1 != iesp2)
            somphi += coef[iesp1 * (ngaz - 1) + iesp2] * yvol[iesp2] / yvol[iesp1];

        CS_F_(mu)->val[iel] += visesp[iesp1] / (1. + somphi);
      }

      /* compute specific heat : J/(kg degres) */
      if (cs_glob_fluid_properties->icp > 0)
      {
        CS_F_(cp)->val[iel] = 0.;
        for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
          CS_F_(cp)->val[iel] += ym[iesp1] * cpesp[iesp1];
      }

      /* compute Lambda/Cp : kg/(m s) */
      if (ifcvsl >= 0)
      {
        for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
          for (int iesp2 = 0; iesp2 < ngaz; iesp2++)
          {
            coef[iesp1 * (ngaz - 1) + iesp2] = 1. + sqrt(xlabes[iesp1] / xlabes[iesp2]) *
                                                    sqrt(sqrt(roesp[iesp2] / roesp[iesp1]));
            coef[iesp1 * (ngaz - 1) + iesp2] *= coef[iesp1 * (ngaz - 1) + iesp2];
            coef[iesp1 * (ngaz - 1) + iesp2] /= (sqrt(1. + roesp[iesp1] / roesp[iesp2]) * sqrt(8.));
          }
        /* Lambda */
        diff_th->val[iel] = 0.;

        for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
        {
          double somphi = 0.;
          for (int iesp2 = 0; iesp2 < ngaz; iesp2++)
            if (iesp1 != iesp2)
              somphi += coef[iesp1 * (ngaz - 1) + iesp2] * yvol[iesp2] / yvol[iesp1];

          diff_th->val[iel] += xlabes[iesp1] / (1. + 1.065 * somphi);
        }

        /* Lambda/Cp */
        if (cs_glob_fluid_properties->icp <= 0)
          diff_th->val[iel] /= cs_glob_fluid_properties->cp0;
        else
          diff_th->val[iel] /= CS_F_(cp)->val[iel];
      }

      /* compute electric conductivity : S/m */
      if (ifcsig >= 0)
      {
        c_prop->val[iel] = 0.;
        double val = 0.;

        for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
          val += yvol[iesp1] / sigesp[iesp1];

        c_prop->val[iel] = 1. / val;
      }

      /* compute radiative transfer : W/m3 */
      if (cs_glob_elec_option->ixkabe == 1)
      {
        CS_F_(absco)->val[iel] = 0.;
        double val = 0.;

        for (int iesp1 = 0; iesp1 < ngaz; iesp1++)
          val += yvol[iesp1] * xkabes[iesp1];

        CS_F_(absco)->val[iel] = val;
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
    }

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

  /* not used yet */
  if (cs_glob_elec_option->ielion > 0)
  {
    /* compute density */
    for (int iel = 0; iel < ncel; iel++)
      CS_F_(rho)->val[iel] = 1.;

    /* compute molecular viscosity : kg/(m s) */
    for (int iel = 0; iel < ncel; iel++)
      CS_F_(mu)->val[iel] = 1.e-2;

    /* compute specific heat : J/(kg degres) */
    for (int iel = 0; iel < ncel; iel++)
      CS_F_(cp)->val[iel] = 1000.;

    /* compute Lambda/Cp : kg/(m s) */
    if (ifcvsl >= 0)
    {
      if (cs_glob_fluid_properties->icp <= 0)
        for (int iel = 0; iel < ncel; iel++)
          CS_F_(mu_t)->val[iel] = 1. / cs_glob_fluid_properties->cp0;
      else
        for (int iel = 0; iel < ncel; iel++)
          CS_F_(mu_t)->val[iel] = 1 / CS_F_(cp)->val[iel];
    }

    /* diffusivity for other properties
     * nothing to do
     * no other properties in this case */
  }

  /* now user properties (for joule effect particulary) */
  cs_user_physical_properties(mesh, mesh_quantities);
}


/*----------------------------------------------------------------------------
 * compute specific electric arc fields
 *----------------------------------------------------------------------------*/

void
cs_compute_electric_field(const cs_mesh_t *mesh,
                          const cs_mesh_quantities_t *mesh_quantities,
                          cs_int_t iappel)
{
  cs_lnum_t  ncel   = mesh->n_cells;
  cs_lnum_t  ncelet = mesh->n_cells_with_ghosts;
  const int keysca  = cs_field_key_id("scalar_diffusivity_id");

  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;

  /* if listing printing is needed */
  int modntl = 0;
  //TODO : control listing output

  /* Reconstructed value */
  cs_real_3_t *grad;
  BFT_MALLOC(grad, ncelet, cs_real_3_t);

  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  /* ----------------------------------------------------- */
  /* first call : J, E => J.E                              */
  /* ----------------------------------------------------- */
  if (iappel == 1)
  {
    /* compute grad(potR) */

    /* Get the calculation option from the field */
    cs_field_get_key_struct(CS_F_(potr), key_cal_opt_id, &var_cal_opt);

    cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                               &gradient_type,
                               &halo_type);

    cs_field_gradient_scalar(CS_F_(potr),
                             false, /* use_previous_t */
                             gradient_type,
                             halo_type,
                             1,    /* inc */
                             true, /* recompute_cocg */
                             grad);

    /* compute electric field E = - grad (potR) */

    /* compute current density j = sig E */
    int diff_id = cs_field_get_key_int(CS_F_(potr), keysca);
    cs_field_t *c_prop = NULL;
    if (diff_id > -1)
      c_prop = cs_field_by_id(diff_id);

    if (cs_glob_elec_option->ieljou > 0 || cs_glob_elec_option->ielarc > 0)
    {
      for (int iel = 0; iel < ncel; iel++)
      {
        CS_FI_(curre, 0)->val[iel] = -c_prop->val[iel] * grad[iel][0];
        CS_FI_(curre, 1)->val[iel] = -c_prop->val[iel] * grad[iel][1];
        CS_FI_(curre, 2)->val[iel] = -c_prop->val[iel] * grad[iel][2];
      }
    }

    /* compute joule effect : j . E */
    for (int iel = 0; iel < ncel; iel++)
    {
      CS_F_(joulp)->val[iel] =  c_prop->val[iel] *
                               (grad[iel][0] * grad[iel][0] +
                                grad[iel][1] * grad[iel][1] +
                                grad[iel][2] * grad[iel][2]);
    }

    /* compute min max for E and J */
    if (modntl == 0)
    {
      bft_printf("-----------------------------------------                    \n");
      bft_printf("   Variable         Minimum       Maximum                    \n");
      bft_printf("-----------------------------------------                    \n");

      /* Grad PotR = -E */
      double vrmin, vrmax;

      for (int i = 0; i < 3; i++) {
        vrmin = grad[0][i];
        vrmax = grad[0][i];

        for (int iel = 0; iel < ncel; iel++)
        {
          vrmin = CS_MIN(vrmin, grad[iel][i]);
          vrmax = CS_MAX(vrmax, grad[iel][i]);
        }

        cs_parall_min(1, CS_DOUBLE, &vrmin);
        cs_parall_max(1, CS_DOUBLE, &vrmax);
        if (i == 0)
          bft_printf("v  Gr_PotRX    %12.5E  %12.5E\n", vrmin, vrmax);
        else if (i == 1)
          bft_printf("v  Gr_PotRY    %12.5E  %12.5E\n", vrmin, vrmax);
        else if (i == 2)
          bft_printf("v  Gr_PotRZ    %12.5E  %12.5E\n", vrmin, vrmax);
      }

      /* current real */
      for (int i = 0; i < 3; i++) {
        vrmin = -c_prop->val[0] * grad[0][i];
        vrmax = -c_prop->val[0] * grad[0][i];

        for (int iel = 0; iel < ncel; iel++)
        {
          vrmin = CS_MIN(vrmin, -c_prop->val[iel] * grad[iel][i]);
          vrmax = CS_MAX(vrmax, -c_prop->val[iel] * grad[iel][i]);
        }

        cs_parall_min(1, CS_DOUBLE, &vrmin);
        cs_parall_max(1, CS_DOUBLE, &vrmax);

        if (i == 0)
          bft_printf("v  Cour_ReX    %12.5E  %12.5E\n", vrmin, vrmax);
        else if (i == 1)
          bft_printf("v  Cour_ReY    %12.5E  %12.5E\n", vrmin, vrmax);
        else if (i == 2)
          bft_printf("v  Cour_ReZ    %12.5E  %12.5E\n", vrmin, vrmax);
      }
      bft_printf("-----------------------------------------                    \n");
    }

    if (cs_glob_elec_option->ieljou == 2 || cs_glob_elec_option->ieljou == 4)
    {
      /* compute grad(potI) */

      /* Get the calculation option from the field */
      cs_field_get_key_struct(CS_F_(poti), key_cal_opt_id, &var_cal_opt);

      cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                                &gradient_type,
                                &halo_type);

      cs_field_gradient_scalar(CS_F_(poti),
                               false, /* use_previous_t */
                               gradient_type,
                               halo_type,
                               1,    /* inc */
                               true, /* recompute_cocg */
                               grad);

      /* compute electric field E = - grad (potI) */

      /* compute current density j = sig E */
      int diff_id_i = cs_field_get_key_int(CS_F_(poti), keysca);
      cs_field_t *c_propi = NULL;
      if (diff_id_i > -1)
        c_propi = cs_field_by_id(diff_id_i);

      if (cs_glob_elec_option->ieljou == 4)
      {
        for (int iel = 0; iel < ncel; iel++)
        {
          CS_FI_(curim, 0)->val[iel] = -c_propi->val[iel] * grad[iel][0];
          CS_FI_(curim, 1)->val[iel] = -c_propi->val[iel] * grad[iel][1];
          CS_FI_(curim, 2)->val[iel] = -c_propi->val[iel] * grad[iel][2];
        }
      }

      /* compute joule effect : j . E */
      for (int iel = 0; iel < ncel; iel++)
      {
        CS_F_(joulp)->val[iel] +=  c_propi->val[iel] *
                                  (grad[iel][0] * grad[iel][0] +
                                   grad[iel][1] * grad[iel][1] +
                                   grad[iel][2] * grad[iel][2]);
      }

      /* compute min max for E and J */
      if (modntl == 0)
      {
        /* Grad PotR = -Ei */
        double vrmin, vrmax;

        for (int i = 0; i < 3; i++)
        {
          vrmin = grad[0][i];
          vrmax = grad[0][i];

          for (int iel = 0; iel < ncel; iel++)
          {
            vrmin = CS_MIN(vrmin, grad[iel][0]);
            vrmax = CS_MAX(vrmax, grad[iel][0]);
          }

          cs_parall_min(1, CS_DOUBLE, &vrmin);
          cs_parall_max(1, CS_DOUBLE, &vrmax);

          if (i == 0)
            bft_printf("v  Gr_PotIX    %12.5E  %12.5E\n", vrmin, vrmax);
          else if (i == 1)
            bft_printf("v  Gr_PotIY    %12.5E  %12.5E\n", vrmin, vrmax);
          else if (i == 2)
            bft_printf("v  Gr_PotIZ    %12.5E  %12.5E\n", vrmin, vrmax);
        }

        /* current imaginary */
        for (int i = 0; i < 3; i++)
        {
          vrmin = -c_propi->val[0] * grad[0][i];
          vrmax = -c_propi->val[0] * grad[0][i];

          for (int iel = 0; iel < ncel; iel++)
          {
            vrmin = CS_MIN(vrmin, -c_propi->val[iel] * grad[iel][i]);
            vrmax = CS_MAX(vrmax, -c_propi->val[iel] * grad[iel][i]);
          }

          cs_parall_min(1, CS_DOUBLE, &vrmin);
          cs_parall_max(1, CS_DOUBLE, &vrmax);

          if (i == 0)
            bft_printf("v  Cour_ImX    %12.5E  %12.5E\n", vrmin, vrmax);
          else if (i == 1)
            bft_printf("v  Cour_ImY    %12.5E  %12.5E\n", vrmin, vrmax);
          else if (i == 2)
            bft_printf("v  Cour_ImZ    %12.5E  %12.5E\n", vrmin, vrmax);

        }
      }
    }
  }

  /* ----------------------------------------------------- */
  /* second call : A, B, JXB                               */
  /* ----------------------------------------------------- */
  else if (iappel == 2) {

    double *Bx, *By, *Bz;
    BFT_MALLOC(Bx, ncelet, double);
    BFT_MALLOC(By, ncelet, double);
    BFT_MALLOC(Bz, ncelet, double);

    if (cs_glob_elec_option->ielarc == 2)
    {
      /* compute magnetic field component B */
      //cs_field_get_key_struct(CS_FI_(potva, 0), key_cal_opt_id, &var_cal_opt);
      cs_field_t  *fp1 = cs_field_by_name_try("vec_potential_01");
      cs_field_t  *fp2 = cs_field_by_name_try("vec_potential_02");
      cs_field_t  *fp3 = cs_field_by_name_try("vec_potential_03");
      cs_field_get_key_struct(fp1, key_cal_opt_id, &var_cal_opt);

      cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                                 &gradient_type,
                                 &halo_type);

      /* Ax component */
      //cs_field_gradient_scalar(CS_FI_(potva, 0),
      cs_field_gradient_scalar(fp1,
                               false, /* use_previous_t */
                               gradient_type,
                               halo_type,
                               1,    /* inc */
                               true, /* recompute_cocg */
                               grad);

      for (int iel = 0; iel < ncel; iel++)
      {
        Bx[iel] =  0.;
        By[iel] =  grad[iel][2];
        Bz[iel] = -grad[iel][1];
      }

      /* Ay component */
      //cs_field_gradient_scalar(CS_FI_(potva, 1),
      cs_field_gradient_scalar(fp2,
                               false, /* use_previous_t */
                               gradient_type,
                               halo_type,
                               1,    /* inc */
                               true, /* recompute_cocg */
                               grad);

      for (int iel = 0; iel < ncel; iel++)
      {
        Bx[iel] -=  grad[iel][2];
        By[iel] +=  0.;
        Bz[iel] +=  grad[iel][0];
      }

      /* Az component */
      //cs_field_gradient_scalar(CS_FI_(potva, 2),
      cs_field_gradient_scalar(fp3,
                               false, /* use_previous_t */
                               gradient_type,
                               halo_type,
                               1,    /* inc */
                               true, /* recompute_cocg */
                               grad);

      for (int iel = 0; iel < ncel; iel++)
      {
        Bx[iel] +=  grad[iel][1];
        By[iel] -=  grad[iel][0];
        Bz[iel] +=  0.;
      }
    }
    else if (cs_glob_elec_option->ielarc == 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Error electric arc with ampere theorem not available\n"));

    /* compute laplace effect j x B */
    if (cs_glob_elec_option->ielarc > 0)
      for (int iel = 0; iel < ncel; iel++)
      {
        CS_FI_(laplf, 0)->val[iel] = CS_FI_(curre, 1)->val[iel] * Bz[iel] -
                                     CS_FI_(curre, 2)->val[iel] * By[iel];
        CS_FI_(laplf, 1)->val[iel] = CS_FI_(curre, 2)->val[iel] * Bx[iel] -
                                     CS_FI_(curre, 0)->val[iel] * Bz[iel];
        CS_FI_(laplf, 2)->val[iel] = CS_FI_(curre, 0)->val[iel] * By[iel] -
                                     CS_FI_(curre, 1)->val[iel] * Bx[iel];
      }

    /* compute min max for B */
    if (cs_glob_elec_option->ielarc > 1)
    {
      if (modntl == 0)
      {
        /* Grad PotR = -E */
        double vrmin, vrmax;
        vrmin = Bx[0];
        vrmax = Bx[0];

        for (int iel = 0; iel < ncel; iel++)
        {
          vrmin = CS_MIN(vrmin, Bx[iel]);
          vrmax = CS_MAX(vrmax, Bx[iel]);
        }

        cs_parall_min(1, CS_DOUBLE, &vrmin);
        cs_parall_max(1, CS_DOUBLE, &vrmax);

        bft_printf("v  Ch_MagX    %12.5E  %12.5E\n", vrmin, vrmax);

        vrmin = By[0];
        vrmax = By[0];

        for (int iel = 0; iel < ncel; iel++)
        {
          vrmin = CS_MIN(vrmin, By[iel]);
          vrmax = CS_MAX(vrmax, By[iel]);
        }

        cs_parall_min(1, CS_DOUBLE, &vrmin);
        cs_parall_max(1, CS_DOUBLE, &vrmax);

        bft_printf("v  Ch_MagY    %12.5E  %12.5E\n", vrmin, vrmax);

        vrmin = Bz[0];
        vrmax = Bz[0];

        for (int iel = 0; iel < ncel; iel++)
        {
          vrmin = CS_MIN(vrmin, Bz[iel]);
          vrmax = CS_MAX(vrmax, Bz[iel]);
        }

        cs_parall_min(1, CS_DOUBLE, &vrmin);
        cs_parall_max(1, CS_DOUBLE, &vrmax);

        bft_printf("v  Ch_MagZ    %12.5E  %12.5E\n", vrmin, vrmax);
      }
    }
    BFT_FREE(Bx);
    BFT_FREE(By);
    BFT_FREE(Bz);
  }

  /* Free memory */
  BFT_FREE(grad);
}


/*----------------------------------------------------------------------------
 * compute source terms for energie and vector potential
 *----------------------------------------------------------------------------*/

void
cs_elec_source_terms(const cs_mesh_t *mesh,
                     const cs_mesh_quantities_t *mesh_quantities,
                     const cs_int_t   f_id,
                           cs_real_t *smbrs)
{
  const cs_field_t  *f    = cs_field_by_id(f_id);
  const char        *name = f->name;
  cs_lnum_t  ncel         = mesh->n_cells;
  cs_lnum_t  ncelet       = mesh->n_cells_with_ghosts;
  double    *volume       = mesh_quantities->cell_vol;

  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

  double *w1;
  BFT_MALLOC(w1, ncelet, double);

  /* enthalpy source term */
  if (strcmp(name, "enthalpy") == 0)
  {
    if (var_cal_opt.iwarni > 0)
      bft_printf("compute source terms for variable : %s\n", name);

    if (cs_glob_time_step->nt_cur > 2)
    {
      for (int iel = 0; iel < ncel; iel++)
        w1[iel] = CS_F_(joulp)->val[iel] * volume[iel];

      if (cs_glob_elec_option->ielarc >= 1)
        if (cs_glob_elec_option->ixkabe == 2)
          for (int iel = 0; iel < ncel; iel++)
            w1[iel] -= CS_F_(radsc)->val[iel] * volume[iel];

      for (int iel = 0; iel < ncel; iel++)
        smbrs[iel] += w1[iel];

      if (var_cal_opt.iwarni > 1)
      {
        double valmin = w1[0];
        double valmax = w1[0];

        for (int iel = 0; iel < ncel; iel++)
        {
          valmin = CS_MIN(valmin, w1[iel]);
          valmax = CS_MAX(valmax, w1[iel]);
        }

        cs_parall_min(1, CS_DOUBLE, &valmin);
        cs_parall_max(1, CS_DOUBLE, &valmax);

        bft_printf(" source terms for H min= %14.5E, max= %14.5E\n", valmin, valmax);
      }
    }
  }

  /* source term for potential vector */
  if (cs_glob_elec_option->ielarc >= 2)
  {
      if (strcmp(name, "vec_potential_01") == 0)
      {
        if (var_cal_opt.iwarni > 0)
          bft_printf("compute source terms for variable : %s\n", name);

        for (int iel = 0; iel < ncel; iel++)
          smbrs[iel] += cs_elec_permvi * CS_FI_(curre, 0)->val[iel] * volume[iel];
      }
      else if (strcmp(name, "vec_potential_02") == 0) {
        if (var_cal_opt.iwarni > 0)
          bft_printf("compute source terms for variable : %s\n", name);

        for (int iel = 0; iel < ncel; iel++)
          smbrs[iel] += cs_elec_permvi * CS_FI_(curre, 1)->val[iel] * volume[iel];
      }
      else if (strcmp(name, "vec_potential_03") == 0) {
        if (var_cal_opt.iwarni > 0)
          bft_printf("compute source terms for variable : %s\n", name);

        for (int iel = 0; iel < ncel; iel++)
          smbrs[iel] += cs_elec_permvi * CS_FI_(curre, 2)->val[iel] * volume[iel];
      }
  }

  BFT_FREE(w1);
  return;
}

/*----------------------------------------------------------------------------
 * add variables fields
 *----------------------------------------------------------------------------*/

void
cs_elec_add_variable_fields(const cs_int_t *ielarc,
                            const cs_int_t *ieljou,
                            const cs_int_t *ielion,
                            const cs_int_t *iihmpr)
{
  cs_field_t *f;
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;
  int dim = 1;
  bool has_previous = true;
  bool interleaved = false;
  double grand = 1.e12;

  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");
  const int kivisl = cs_field_key_id("scalar_diffusivity_id");
  //const int keysca = cs_field_key_id("scalar_id");

  {
    int f_id = cs_variable_field_create("enthalpy", "Enthalpy", CS_MESH_LOCATION_CELLS, dim);
    f = cs_field_by_id(f_id);
    cs_field_set_key_double(f, kscmin, -grand);
    cs_field_set_key_int(f, kivisl, 0);
    int isca = cs_add_model_field_indexes(f->id);

    // set thermal model
    cs_thermal_model_t *thermal = cs_get_glob_thermal_model();
    thermal->itherm = 2;
    thermal->iscalt = isca;
  }

  {
    int f_id = cs_variable_field_create("elec_pot_r", "POT_EL_R", CS_MESH_LOCATION_CELLS, dim);
    f = cs_field_by_id(f_id);
    cs_field_set_key_double(f, kscmin, -grand);
    cs_field_set_key_double(f, kscmax,  grand);
    cs_field_set_key_int(f, kivisl, 0);
    int isca = cs_add_model_field_indexes(f->id);
  }

  if (*ieljou == 2 || *ieljou == 4) {
    int f_id = cs_variable_field_create("elec_pot_i", "POT_EL_I", CS_MESH_LOCATION_CELLS, dim);
    f = cs_field_by_id(f_id);
    cs_field_set_key_double(f, kscmin, -grand);
    cs_field_set_key_double(f, kscmax,  grand);
    cs_field_set_key_int(f, kivisl, 0);
    int isca = cs_add_model_field_indexes(f->id);
  }

  if (*ielarc > 1) {
    {
      int f_id = cs_variable_field_create("vec_potential_01", "POT_VEC1", CS_MESH_LOCATION_CELLS, dim);
      f = cs_field_by_id(f_id);
      cs_field_set_key_double(f, kscmin, -grand);
      cs_field_set_key_double(f, kscmax,  grand);
      cs_field_set_key_int(f, kivisl, -1);
      int isca = cs_add_model_field_indexes(f->id);
    }

    {
      int f_id = cs_variable_field_create("vec_potential_02", "POT_VEC2", CS_MESH_LOCATION_CELLS, dim);
      f = cs_field_by_id(f_id);
      cs_field_set_key_double(f, kscmin, -grand);
      cs_field_set_key_double(f, kscmax,  grand);
      cs_field_set_key_int(f, kivisl, -1);
      int isca = cs_add_model_field_indexes(f->id);
    }

    {
      int f_id = cs_variable_field_create("vec_potential_03", "POT_VEC3", CS_MESH_LOCATION_CELLS, dim);
      f = cs_field_by_id(f_id);
      cs_field_set_key_double(f, kscmin, -grand);
      cs_field_set_key_double(f, kscmax,  grand);
      cs_field_set_key_int(f, kivisl, -1);
      int isca = cs_add_model_field_indexes(f->id);
    }
  }

  if (cs_glob_elec_properties->ngaz > 1) {
    for (int igaz = 0; igaz < cs_glob_elec_properties->ngaz - 1; igaz++) {
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

      int f_id = cs_variable_field_create(name, label, CS_MESH_LOCATION_CELLS, dim);
      f = cs_field_by_id(f_id);

      cs_field_set_key_double(f, kscmin, 0.);
      cs_field_set_key_double(f, kscmax, 1.);
      cs_field_set_key_int(f, kivisl, 0);
      int isca = cs_add_model_field_indexes(f->id);
      BFT_FREE(name);
      BFT_FREE(label);
      BFT_FREE(suf);
    }
  }

  int n_gasses = cs_glob_elec_properties->ngaz;
  cs_field_pointer_map_electric_arcs(n_gasses);

  /* Map labels for GUI */
  if (*iihmpr == 1)
    cs_gui_labels_electric_arcs(n_gasses);

  return;
}

/*----------------------------------------------------------------------------
 * add properties fields
 *----------------------------------------------------------------------------*/

void
cs_elec_add_property_fields(const cs_int_t *ielarc,
                            const cs_int_t *ieljou,
                            const cs_int_t *ielion)
{
  cs_field_t *f;
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  int dim = 1;
  bool has_previous = false;
  bool interleaved = false;
  const int klbl   = cs_field_key_id("label");
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");

  {
    f = cs_field_create("temperature",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        dim,
                        interleaved,
                        has_previous);
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Temperature");

    /* Property number and mapping to field and postprocessing */
    CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
  }

  {
    f = cs_field_create("joule_power",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        dim,
                        interleaved,
                        has_previous);
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "PuisJoul");

    /* Property number and mapping to field and postprocessing */
    CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
  }

  {
    f = cs_field_create("current_re_1",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        dim,
                        interleaved,
                        has_previous);
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Cour_re1");

    /* Property number and mapping to field and postprocessing */
    CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
  }

  {
    f = cs_field_create("current_re_2",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        dim,
                        interleaved,
                        has_previous);
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Cour_re2");

    /* Property number and mapping to field and postprocessing */
    CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
  }

  {
    f = cs_field_create("current_re_3",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        dim,
                        interleaved,
                        has_previous);
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Cour_re3");

    /* Property number and mapping to field and postprocessing */
    CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
  }

  /* specific for joule effect */
  if (*ieljou == 2 || *ieljou == 4) {
    {
      f = cs_field_create("current_im_1",
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          dim,
                          interleaved,
                          has_previous);
      cs_field_set_key_int(f, keyvis, 1);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "CouImag1");

      /* Property number and mapping to field and postprocessing */
      CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
    }

    {
      f = cs_field_create("current_im_2",
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          dim,
                          interleaved,
                          has_previous);
      cs_field_set_key_int(f, keyvis, 1);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "CouImag2");

      /* Property number and mapping to field and postprocessing */
      CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
    }

    {
      f = cs_field_create("current_im_3",
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          dim,
                          interleaved,
                          has_previous);
      cs_field_set_key_int(f, keyvis, 1);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "CouImag3");

      /* Property number and mapping to field and postprocessing */
      CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
    }
  }

  /* specific for electric arc */
  if (*ielarc > 0) {
    {
      f = cs_field_create("laplace_force_1",
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          dim,
                          interleaved,
                          has_previous);
      cs_field_set_key_int(f, keyvis, 1);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "For_Lap1");

      /* Property number and mapping to field and postprocessing */
      CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
    }

    {
      f = cs_field_create("laplace_force_2",
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          dim,
                          interleaved,
                          has_previous);
      cs_field_set_key_int(f, keyvis, 1);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "For_Lap2");

      /* Property number and mapping to field and postprocessing */
      CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
    }

    {
      f = cs_field_create("laplace_force_3",
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          dim,
                          interleaved,
                          has_previous);
      cs_field_set_key_int(f, keyvis, 1);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "For_Lap3");

      /* Property number and mapping to field and postprocessing */
      CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
    }

    if (cs_glob_elec_option->ixkabe == 1) {
      f = cs_field_create("absorption_coeff",
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          dim,
                          interleaved,
                          has_previous);
      cs_field_set_key_int(f, keyvis, 1);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "Coef_Abso");

      /* Property number and mapping to field and postprocessing */
      CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
    }
    else if (cs_glob_elec_option->ixkabe == 2) {
      f = cs_field_create("radiation_source",
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          dim,
                          interleaved,
                          has_previous);
      cs_field_set_key_int(f, keyvis, 1);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "TS_radia");

      /* Property number and mapping to field and postprocessing */
      CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
    }
  }

  /* specific for inic conduction */
  if (*ielion > 0) {
    f = cs_field_create("elec_charge",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        dim,
                        interleaved,
                        has_previous);
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Charge");

    /* Property number and mapping to field and postprocessing */
    CS_PROCF(add_property_field_post, ADD_PROPERTY_FIELD_POST)(&(f->id), &dim);
  }

  cs_field_pointer_properties_map_electric_arcs();

  return;
}

/*----------------------------------------------------------------------------
 * initialize electric fields
 *----------------------------------------------------------------------------*/

void
cs_elec_fields_initialize(const cs_mesh_t *mesh,
                          const cs_mesh_quantities_t *mesh_quantities,
                                cs_int_t  isuite,
                                cs_int_t  nvar,
                                cs_int_t  nscal,
                                cs_real_t *dt)
{
  BFT_MALLOC(_elec_option->izreca, mesh->n_i_faces, int);
  for (int i = 0; i < mesh->n_i_faces; i++)
    _elec_option->izreca[i] = 0;

  cs_lnum_t  ncel = mesh->n_cells;

  static int ipass = 0;
  ipass += 1;

  double d2s3 = 2. /3.;

  if (isuite == 0 && ipass == 1) {
    double xkent = 1.e-10;
    double xeent = 1.e-10;

    if (cs_glob_turb_model->itytur == 2) {
      for (int iel = 0; iel < ncel; iel++) {
        CS_F_(k)->val[iel] = xkent;
        CS_F_(eps)->val[iel] = xeent;
      }
    }
    else if (cs_glob_turb_model->itytur == 3) {
      for (int iel = 0; iel < ncel; iel++) {
        CS_F_(r11)->val[iel] = d2s3 * xkent;
        CS_F_(r22)->val[iel] = d2s3 * xkent;
        CS_F_(r33)->val[iel] = d2s3 * xkent;
        CS_F_(r12)->val[iel] = 0.;
        CS_F_(r13)->val[iel] = 0.;
        CS_F_(r23)->val[iel] = 0.;
        CS_F_(eps)->val[iel] = xeent;
      }
    }
    else if (cs_glob_turb_model->itytur == 5) {
      for (int iel = 0; iel < ncel; iel++) {
        CS_F_(k)->val[iel] = xkent;
        CS_F_(eps)->val[iel] = xeent;
        CS_F_(phi)->val[iel] = d2s3;
        CS_F_(f_bar)->val[iel] = 0.;
      }
    }
    else if (cs_glob_turb_model->itytur == 6) {
      for (int iel = 0; iel < ncel; iel++) {
        CS_F_(k)->val[iel] = xkent;
        CS_F_(omg)->val[iel] = xeent / cs_turb_cmu / xkent;
      }
    }
    else if (cs_glob_turb_model->itytur == 7) {
      for (int iel = 0; iel < ncel; iel++) {
        CS_F_(nusa)->val[iel] = cs_turb_cmu * xkent * xkent / xeent;
      }
    }

    /* enthalpy */
    double hinit = 0.;
    if (cs_glob_elec_option->ielarc > 0) {
      double *ym;
      BFT_MALLOC(ym, cs_glob_elec_properties->ngaz, double);
      ym[0] = 1.;
      if (cs_glob_elec_properties->ngaz > 1)
        for (int i = 1; i < cs_glob_elec_properties->ngaz; i++)
          ym[i] = 0.;

      double tinit = cs_glob_fluid_properties->t0;
      cs_elec_convert_h_t(-1, ym, &hinit, &tinit);
      BFT_FREE(ym);
    }

    for (int iel = 0; iel < ncel; iel++) {
      CS_F_(h)->val[iel] = hinit;
    }

    /* volumic fraction */
    if (cs_glob_elec_properties->ngaz > 1) {
      for (int iel = 0; iel < ncel; iel++)
        CS_FI_(ycoel, 0)->val[iel] = 1.;

      for (int igaz = 1; igaz < cs_glob_elec_properties->ngaz - 1; igaz++)
        for (int iel = 0; iel < ncel; iel++)
          CS_FI_(ycoel, igaz)->val[iel] = 0.;
    }

    /* electric potential */
    for (int iel = 0; iel < ncel; iel++)
      CS_F_(potr)->val[iel] = 0.;

    if (cs_glob_elec_option->ieljou == 2 || cs_glob_elec_option->ieljou == 4)
      for (int iel = 0; iel < ncel; iel++)
        CS_F_(poti)->val[iel] = 0.;

    /* potential vector */
    if (cs_glob_elec_option->ielarc > 1) {
      /* TODO when vector field
       * for (int i = 0; i < 3; i++)
       * for (int iel = 0; iel < ncel; iel++)
       *   CS_FI_(potva, i)->val[iel] = 0.;
       */
      cs_field_t  *fp1 = cs_field_by_name_try("vec_potential_01");
      cs_field_t  *fp2 = cs_field_by_name_try("vec_potential_02");
      cs_field_t  *fp3 = cs_field_by_name_try("vec_potential_03");
      for (int iel = 0; iel < ncel; iel++) {
        fp1->val[iel] = 0.;
        fp2->val[iel] = 0.;
        fp3->val[iel] = 0.;
      }
    }

    /* source term */
    for (int iel = 0; iel < ncel; iel++)
      CS_F_(joulp)->val[iel] = 0.;

    if (cs_glob_elec_option->ielarc > 1)
      for (int i = 0; i < 3; i++)
        for (int iel = 0; iel < ncel; iel++)
          CS_FI_(laplf, i)->val[iel] = 0.;
  }

  /* user function */
  if (ipass == 1)
    CS_PROCF(cs_user_initialization, CS_USER_INITIALIZATION)(nvar, nscal, dt);
}

/*----------------------------------------------------------------------------
 * scaling electric quantities
 *----------------------------------------------------------------------------*/

void
cs_elec_scaling_function(const cs_mesh_t *mesh,
                         const cs_mesh_quantities_t *mesh_quantities,
                                cs_real_t *dt)
{
  cs_real_t *volume = mesh_quantities->cell_vol;
  cs_real_t *surfac = mesh_quantities->i_face_normal;
  cs_lnum_t  ncel   = mesh->n_cells;
  cs_lnum_t  nfac   = mesh->n_i_faces;

  double coepot = 0.;
  double coepoa = 1.;

  if (cs_glob_elec_option->ielarc >= 1) {
    if (cs_glob_elec_option->modrec == 1) {
      /* standard model */
      double somje = 0.;
      for (int iel = 0; iel < ncel; iel++)
        somje += CS_F_(joulp)->val[iel] * volume[iel];

      cs_parall_sum(1, CS_DOUBLE, &somje);

      coepot = cs_glob_elec_option->couimp * cs_glob_elec_option->pot_diff
              / CS_MAX(somje, epzero);
      coepoa = coepot;

      if (coepot > 1.5)
        coepot = 1.5;
      if (coepot < 0.75)
        coepot = 0.75;

      bft_printf("imposed current / current %14.5E, scaling coef. %14.5E\n", coepoa, coepot);
    }
    else if (cs_glob_elec_option->modrec == 2) {
      /* restrike model */
      uielrc();
      double elcou = 0.;
      for (int ifac = 0; ifac < nfac; ifac++) {
        if (cs_glob_elec_option->izreca[ifac] > 0) {
          bool ok = true;
          for (int idir = 0; idir < 3; idir++)
            if (fabs(surfac[3 * ifac + idir]) > 0. && idir != (cs_glob_elec_option->idreca - 1))
              ok = false;
          if (ok) {
            int iel = mesh->i_face_cells[ifac][0];
            elcou += CS_FI_(curre, cs_glob_elec_option->idreca - 1)->val[iel]
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
      _elec_option->elcou = elcou;
    }

    if (cs_glob_elec_option->modrec == 1 ||
        cs_glob_elec_option->modrec == 2) {
      double dtj = 1.e15;
      double dtjm = dtj;
      double delhsh = 0.;
      double cdtj = 20.;

      for (int iel = 0; iel < ncel; iel++) {
        if (CS_F_(rho)->val[iel] != 0.)
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
      bft_printf(" Dpot recale   = %14.5E\n", _elec_option->pot_diff * coepot);

      /* scaling electric fields */
      _elec_option->pot_diff *= coepot;

      /* electric potential (for post treatment) */
      for (int iel = 0; iel < ncel; iel++)
        CS_F_(potr)->val[iel] *= coepot;

      /* current density */
      if (cs_glob_elec_option->ielarc > 0)
        for (int i = 0; i < 3; i++)
          for (int iel = 0; iel < ncel; iel++)
            CS_FI_(curre, i)->val[iel] *= coepot;

      /* joule effect */
      for (int iel = 0; iel < ncel; iel++)
        CS_F_(joulp)->val[iel] *= coepot * coepot;
    }
  }

  /* joule effect */
  if (cs_glob_elec_option->ieljou > 0) {
    /* standard model */
    double somje = 0.;
    for (int iel = 0; iel < ncel; iel++)
      somje += CS_F_(joulp)->val[iel] * volume[iel];

    cs_parall_sum(1, CS_DOUBLE, &somje);

    coepot = cs_glob_elec_option->puisim / CS_MAX(somje, epzero);
    double coefav = coepot;

    if (coepot > 1.5)
      coepot = 1.5;
    if (coepot < 0.75)
      coepot = 0.75;

    bft_printf("imposed power / sum(jE) %14.5E, scaling coef. %14.5E\n", coefav, coepot);

    /* scaling electric fields */
    _elec_option->pot_diff *= coepot;
    _elec_option->coejou *= coepot;

    /* electric potential (for post treatment) */
    if (cs_glob_elec_option->ieljou != 3 &&
        cs_glob_elec_option->ieljou != 4)
      for (int iel = 0; iel < ncel; iel++)
        CS_F_(potr)->val[iel] *= coepot;

    /* current density */
    if (cs_glob_elec_option->ieljou == 2)
      for (int i = 0; i < 3; i++)
        for (int iel = 0; iel < ncel; iel++)
          CS_F_(poti)->val[iel] *= coepot;

    /* joule effect */
    for (int iel = 0; iel < ncel; iel++)
      CS_F_(joulp)->val[iel] *= coepot * coepot;
  }

  cs_user_scaling_elec(mesh, mesh_quantities, dt);
  return;
}

/*----------------------------------------------------------------------------*/

void
CS_PROCF (elini1, ELINI1) (      cs_real_t *srrom,
                                 cs_real_t *visls0,
                                 cs_real_t *diftl0,
                                 cs_int_t  *iconv,
                                 cs_int_t  *istat,
                                 cs_int_t  *idiff,
                                 cs_int_t  *idifft,
                                 cs_int_t  *idircl,
                                 cs_int_t  *isca,
                                 cs_real_t *blencv,
                                 cs_real_t *sigmas,
                                 cs_int_t  *iwarni,
                           const cs_int_t  *iihmpr)
{
  /* initialization */
  cs_electrical_model_specific_initialization(srrom, visls0, diftl0, iconv, istat,
                                              idiff, idifft, idircl, isca, blencv,
                                              sigmas, iwarni, *iihmpr);
}

void
CS_PROCF (elflux, ELFLUX) (cs_int_t *iappel)
{
  cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;
  cs_compute_electric_field(mesh, mesh_quantities, *iappel);
}

void
CS_PROCF (elthht, ELTHHT) (cs_int_t  *mode,
                           cs_real_t *ym,
                           cs_real_t *enthal,
                           cs_real_t *temp)
{
  cs_elec_convert_h_t(*mode, ym, enthal, temp);
}

void
CS_PROCF (ellecd, ELLECD) (cs_int_t *ieljou,
                           cs_int_t *ielarc,
                           cs_int_t *ielion)
{
  cs_electrical_model_initialize(*ielarc, *ieljou, *ielion);
  cs_electrical_properties_read(*ielarc, *ieljou);
}

void
CS_PROCF (elphyv, ELPHYV) (cs_real_t *srrom)
{
  cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;
  cs_elec_physical_properties(mesh, mesh_quantities, *srrom);
}

void
CS_PROCF (eltssc, ELTSSC) (const cs_int_t  *isca,
                                 cs_real_t *smbrs)
{
  cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;
  const int keysca = cs_field_key_id("scalar_id");

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) == *isca)
      cs_elec_source_terms(mesh, mesh_quantities, f->id, smbrs);
  }
}

void
CS_PROCF (elvarp, ELVARP) (cs_int_t *ieljou,
                           cs_int_t *ielarc,
                           cs_int_t *ielion,
                           cs_int_t *iihmpr)
{
  cs_elec_add_variable_fields(ielarc, ieljou, ielion, iihmpr);
}

void
CS_PROCF (elprop, ELPROP) (cs_int_t *ieljou,
                           cs_int_t *ielarc,
                           cs_int_t *ielion)
{
  cs_elec_add_property_fields(ielarc, ieljou, ielion);
}

void
CS_PROCF (eliniv, ELINIV) (cs_int_t  *isuite,
                           cs_int_t  *nvar,
                           cs_int_t  *nscal,
                           cs_real_t *dt)
{
  cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;
  cs_elec_fields_initialize(mesh, mesh_quantities, *isuite, *nvar, *nscal, dt);
}

void
CS_PROCF (elreca, ELRECA) (cs_real_t *dt)
{
  cs_mesh_t *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;
  cs_elec_scaling_function(mesh, mesh_quantities, dt);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global electric model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ngazge         --> pointer to cs_glob_elec_properties->ngaz
 *   ielcor         --> pointer to cs_glob_elec_option->ielcor
 *   pot_diff       --> pointer to cs_glob_elec_option->pot_diff
 *   coejou         --> pointer to cs_glob_elec_option->coejou
 *   elcou          --> pointer to cs_glob_elec_option->elcou
 *   irestrike      --> pointer to cs_glob_elec_option->irestrike
 *   ntdcla         --> pointer to cs_glob_elec_option->ntdcla
 *   restrike_point --> pointer to cs_glob_elec_option->restrike_point
 *----------------------------------------------------------------------------*/

void
cs_f_elec_model_get_pointers(int     **ngazge,
                             int     **ielcor,
                             double  **pot_diff,
                             double  **coejou,
                             double  **elcou,
                             int     **irestrike,
                             int     **ntdcla,
                             double  **restrike_pointX,
                             double  **restrike_pointY,
                             double  **restrike_pointZ)
{
  *ngazge          = &(_elec_properties->ngaz);
  *ielcor          = &(_elec_option->ielcor);
  *pot_diff        = &(_elec_option->pot_diff);
  *coejou          = &(_elec_option->coejou);
  *elcou           = &(_elec_option->elcou);
  *irestrike       = &(_elec_option->irestrike);
  *ntdcla          = &(_elec_option->ntdcla);
  *restrike_pointX = &(_elec_option->restrike_point[0]);
  *restrike_pointY = &(_elec_option->restrike_point[1]);
  *restrike_pointZ = &(_elec_option->restrike_point[2]);
}

END_C_DECLS


