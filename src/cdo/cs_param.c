/*============================================================================
 * Routines to handle the definition and usage of material properties
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_mesh_location.h"
#include "cs_field.h"
#include "cs_cdo.h"
#include "cs_cdo_toolbox.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Static local variables
 *============================================================================*/

static int cs_param_n_properties = 0;
static cs_param_pty_t  *cs_param_properties = NULL;
static int cs_param_n_adv_fields = 0;
static cs_param_adv_field_t  *cs_param_adv_fields = NULL;

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_PARAM_PTY_DBG  0

/*============================================================================
 * Local variables
 *============================================================================*/

static const char
cs_param_def_type_name[CS_PARAM_N_DEF_TYPES][CS_CDO_LEN_NAME]=
  { N_("by value"),
    N_("by field"),
    N_("by evaluator"),
    N_("by analytic function"),
    N_("by user function"),
    N_("by law"),
    N_("by file") };

static const char
cs_param_var_type_name[CS_PARAM_N_VAR_TYPES][CS_CDO_LEN_NAME]=
  { N_("scalar"),
    N_("vector"),
    N_("tensor") };

static const char
cs_param_pty_type_name[CS_PARAM_N_PTY_TYPES][CS_CDO_LEN_NAME]=
  { N_("isotropic"),
    N_("orthotropic"),
    N_("anisotropic") };

static const char
cs_param_bc_type_name[CS_PARAM_N_BC_TYPES][CS_CDO_LEN_NAME] =
  { N_("Homogeneous Dirichlet"),
    N_("Dirichlet"),
    N_("Homogeneous Neumann"),
    N_("Neumann"),
    N_("Robin") };

static const char
cs_param_boundary_type_name[CS_PARAM_N_BOUNDARY_TYPES][CS_CDO_LEN_NAME] =
  { N_("Wall"),
    N_("Inlet"),
    N_("Outlet"),
    N_("Symmetry") };

static const char
cs_param_hodge_type_desc[CS_PARAM_N_HODGE_TYPES][CS_CDO_LEN_NAME] =
  { "VpCd",
    "EpFd",
    "FpEd",
    "EdFp",
    "CpVd"  };

static const char
cs_param_hodge_algo_desc[CS_PARAM_N_HODGE_ALGOS][CS_CDO_LEN_NAME] =
  { "Voronoi",
    "Whitney on the Barycentric Subdivision (WBS)",
    "COnsistency-STabilization splitting (COST)" };

static const char
cs_param_source_term_type_name[CS_PARAM_N_SOURCE_TERM_TYPES][CS_CDO_LEN_NAME] =
  { N_("user"),
    N_("mass"),
    N_("head loss") };

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the id related to a material property is valid
 *
 * \param[in]     id     id to check
 */
/*----------------------------------------------------------------------------*/

static void
_check_pty_id(int   id)
{
  if (id < 0 || id >= cs_param_n_properties)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid id (equal to %d and should be between 0 and %d).\n"
                " Cannot access to the definition of this material property.\n"),
              id, cs_param_n_properties);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the id related to an advection field is valid
 *
 * \param[in]     id     id to check
 */
/*----------------------------------------------------------------------------*/

static void
_check_adv_field_id(int   id)
{
  if (id < 0 || id >= cs_param_n_adv_fields)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid id (equal to %d and should be between 0 and %d).\n"
                " Cannot access to the definition of this advection field.\n"),
              id, cs_param_n_adv_fields);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a cs_def_t structure
 *
 * \param[in]      def_type   type of definition (by value, function...)
 * \param[in]      var_type   type of variables (scalar, vector, tensor...)
 * \param[in]      val        value to set
 * \param[in, out] def        pointer to cs_def_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_set_def(cs_param_def_type_t      def_type,
                 cs_param_var_type_t      var_type,
                 const void              *val,
                 cs_def_t                *def)
{
  assert(var_type != CS_PARAM_N_VAR_TYPES);

  switch (def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    if (val == NULL) {
      if (var_type == CS_PARAM_VAR_SCAL)
        def->get.val = 0.0;
      else if (var_type == CS_PARAM_VAR_VECT)
        def->get.vect[0] = def->get.vect[1] = def->get.vect[2] = 0.0;
      else if (var_type == CS_PARAM_VAR_TENS) {
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            def->get.tens[i][j] = 0.0;
      }
      else {
        assert(var_type == CS_PARAM_VAR_SYMTENS);
        for (int i = 0; i < 6; i++)
          def->get.twovects[i] = 0.0;
      }
    }
    else { // val != NULL
      if (var_type == CS_PARAM_VAR_SCAL)
        def->get.val = atof(val);
      else if (var_type == CS_PARAM_VAR_VECT) {
        char s[3][32];
        sscanf(val, "%s %s %s", s[0], s[1], s[2]);
        for (int i = 0; i < 3; i++)
          def->get.vect[i] = atof(s[i]);
      }
      else if (var_type == CS_PARAM_VAR_TENS) {
        char s[9][32];
        sscanf(val, "%s %s %s %s %s %s %s %s %s",
               s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8]);
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            def->get.tens[i][j] = atof(s[3*i+j]);
      }
      else {
        assert(var_type == CS_PARAM_VAR_SYMTENS);
        char s[6][32];
        sscanf(val, "%s %s %s %s %s %s", s[0], s[1], s[2], s[3], s[4], s[5]);
        for (int i = 0; i < 6; i++)
          def->get.twovects[i] = atof(s[i]);
      }
    }
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    if (val == NULL)
      def->analytic = NULL;
    else
      def->analytic = (cs_analytic_func_t *)val;
    break;

  case CS_PARAM_DEF_BY_TIME_FUNCTION:
    if (val == NULL)
      def->time_func = NULL;
    else
      def->time_func = (cs_timestep_func_t *)val;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " This type of definition is not handled yet.\n"
              " Please modify your settings.");
    break;

  } /* end of switch on def_type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a cs_param_pty_t structure from its id
 *
 * \param[in]  pty_id   id related to a property
 *
 * \return a pointer to a cs_param_pty_t
 */
/*----------------------------------------------------------------------------*/

cs_param_pty_t *
cs_param_pty_get(int  pty_id)
{
  _check_pty_id(pty_id);

  return cs_param_properties + pty_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the id related to a property definition from its name
 *
 * \param[in]  ref_name    name of the property to find
 *
 * \return -1 if not found otherwise the associated id
 */
/*----------------------------------------------------------------------------*/

int
cs_param_pty_get_id_by_name(const char  *ref_name)
{
  int  i;

  int  pty_id = -1;
  int  reflen = strlen(ref_name);

  for (i = 0; i < cs_param_n_properties; i++) {

    const cs_param_pty_t  *pty = cs_param_properties + i;
    int  len = strlen(pty->name);

    if (reflen == len) {
      if (strcmp(ref_name, pty->name) == 0) {
        pty_id = i;
        break;
      }
    }

  } /* Loop on mesh locations */

  return pty_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add by default several material properties
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_set_default(void)
{
  int  len;
  cs_param_pty_t *pty = NULL;

  /* Sanity check */
  assert(cs_param_n_properties == 0);
  assert(cs_param_properties == NULL);

  BFT_MALLOC(cs_param_properties, 3, cs_param_pty_t);
  cs_param_n_properties = 3;

  /* Warning: Keep this order of definition.
     This knowledge is used to get a faster access to default properties */

  /* Property: Unity */
  pty = cs_param_properties;
  len = strlen("unity")+1;
  BFT_MALLOC(pty->name, len, char);
  strncpy(pty->name, "unity", len);
  pty->flag = CS_PARAM_FLAG_UNIFORM;
  pty->post_freq = -1;
  pty->field_id = -1;
  pty->type = CS_PARAM_PTY_ISO;
  pty->def_type = CS_PARAM_DEF_BY_VALUE;
  pty->def.get.val = 1;

  /* Property: MassDensity */
  pty = cs_param_properties + 1;
  len = strlen("mass_density")+1;
  BFT_MALLOC(pty->name, len, char);
  strncpy(pty->name, "mass_density", len);
  pty->flag = CS_PARAM_FLAG_UNIFORM;
  pty->type = CS_PARAM_PTY_ISO;
  pty->post_freq = -1;
  pty->field_id = -1;
  pty->def_type = CS_PARAM_DEF_BY_VALUE;
  pty->def.get.val = 1;

  /* Property: LaminarViscosity */
  pty = cs_param_properties + 2;
  len = strlen("laminar_viscosity")+1;
  BFT_MALLOC(pty->name, len, char);
  strncpy(pty->name, "laminar_viscosity", len);
  pty->flag = CS_PARAM_FLAG_UNIFORM;
  pty->post_freq = -1;
  pty->field_id = -1;
  pty->type = CS_PARAM_PTY_ISO;
  pty->def_type = CS_PARAM_DEF_BY_VALUE;
  pty->def.get.val = 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a material property structure
 *
 * \param[in]  name        name of the material property
 * \param[in]  key_type    keyname of the type of property
 * \param[in]  post_freq   -1 (no post-processing), 0 (at the beginning)
 *                         otherwise every post_freq iteration(s)
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_add(const char      *name,
                 const char      *key_type,
                 int              post_freq)
{
  int  pty_id = cs_param_pty_get_id_by_name(name);
  cs_param_pty_t  *pty = NULL;

  if (pty_id > -1) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_(" An existing property already has the same name %s.\n"
                 " Stop adding the material property.\n"), name);
    return;
  }

  /* Allocate and initialize a new material property */
  BFT_REALLOC(cs_param_properties,
              cs_param_n_properties + 1,
              cs_param_pty_t);

  pty = cs_param_properties + cs_param_n_properties;
  cs_param_n_properties++;

  pty->post_freq = post_freq;

  /* Assign a type */
  if (strcmp(key_type, "isotropic") == 0)
    pty->type = CS_PARAM_PTY_ISO;
  else if (strcmp(key_type, "orthotropic") == 0)
    pty->type = CS_PARAM_PTY_ORTHO;
  else if (strcmp(key_type, "anisotropic") == 0)
    pty->type = CS_PARAM_PTY_ANISO;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key %s for setting the type of property.\n"
                " Key is one of the following: isotropic, orthotropic or"
                " anisotropic.\n"
                " Please modify your settings."), key_type);

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(pty->name, len, char);
  strncpy(pty->name, name, len);

  /*Default initialization */
  pty->flag = 0;
  pty->def_type = CS_PARAM_N_DEF_TYPES;
  pty->def.get.val = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a way to compute the value of a material property
 *
 * \param[in]  name      name of the material property
 * \param[in]  def_key   way of defining the value of the property
 * \param[in]  val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_set(const char     *name,
                 const char     *def_key,
                 const void     *val)
{
  int  pty_id = cs_param_pty_get_id_by_name(name);

  if (pty_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop setting the material property %s.\n"
                " Do not find a similar name in the property database.\n"),
              name);

  cs_param_pty_t  *pty = cs_param_properties + pty_id;

  /* Get the type of definition */
  if (strcmp(def_key, "value") == 0) {
    pty->def_type =  CS_PARAM_DEF_BY_VALUE;
    pty->flag |= CS_PARAM_FLAG_UNIFORM;
  }
  else if (strcmp(def_key, "field") == 0)
    pty->def_type = CS_PARAM_DEF_BY_FIELD;
  else if (strcmp(def_key, "evaluator") == 0)
    pty->def_type = CS_PARAM_DEF_BY_EVALUATOR;
  else if (strcmp(def_key, "analytic") == 0)
    pty->def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  else if (strcmp(def_key, "user") == 0)
    pty->def_type = CS_PARAM_DEF_BY_USER_FUNCTION;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the type of definition.\n"
                " Given key: %s\n"
                " Choice among value, field, evaluator, analytic, user, law"
                " or file\n"
                " Please modify your settings."), def_key);

  pty->flag |= CS_PARAM_FLAG_SYMMET;  // A material must be symmetric

  switch (pty->type) {

  case CS_PARAM_PTY_ISO:
    cs_param_set_def(pty->def_type, CS_PARAM_VAR_SCAL, val, &(pty->def));
    break;
  case CS_PARAM_PTY_ORTHO:
    cs_param_set_def(pty->def_type, CS_PARAM_VAR_VECT, val, &(pty->def));
    break;
  case CS_PARAM_PTY_ANISO:
    cs_param_set_def(pty->def_type, CS_PARAM_VAR_TENS, val, &(pty->def));

    if (pty->def_type == CS_PARAM_DEF_BY_VALUE) { /* Check the symmetry */

      cs_get_t  get = pty->def.get;
      if ((get.tens[0][1] - get.tens[1][0]) > cs_get_eps_machine() ||
          (get.tens[0][2] - get.tens[2][0]) > cs_get_eps_machine() ||
          (get.tens[1][2] - get.tens[2][1]) > cs_get_eps_machine())
        bft_error(__FILE__, __LINE__, 0,
                  _(" The definition of the tensor related to the material"
                    " property %s is not symmetric.\n"
                    " This case is not handled. Please check your settings.\n"),
                  pty->name);

    }
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of material property."));
    break;

  } /* switch on property type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query to know if the material property is uniform
 *
 * \param[in]    pty_id    id related to the material property to deal with
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

_Bool
cs_param_pty_is_uniform(int       pty_id)
{
  _check_pty_id(pty_id);

  cs_param_pty_t  *pty = cs_param_properties + pty_id;

  if (pty->flag & CS_PARAM_FLAG_UNIFORM)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query to know the type of material property
 *
 * \param[in]    pty_id    id related to the material property to deal with
 *
 * \return  the type of material property
 */
/*----------------------------------------------------------------------------*/

cs_param_pty_type_t
cs_param_pty_get_type(int       pty_id)
{
  _check_pty_id(pty_id);

  cs_param_pty_t  *pty = cs_param_properties + pty_id;

  return  pty->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of a material property from its id
 *
 * \param[in]    pty_id    id related to the material property to deal with
 *
 * \return  the name of the related property
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_pty_get_name(int            pty_id)
{
  _check_pty_id(pty_id);

  cs_param_pty_t  *pty = cs_param_properties + pty_id;

  return pty->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of all the cs_param_pty_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_summary_all(void)
{
  int  pty_id;

  bft_printf("\n");
  bft_printf("%s", lsepline);
  bft_printf("\tSummary of the definition of material properties\n");
  bft_printf("%s", lsepline);

  for (pty_id = 0; pty_id < cs_param_n_properties; pty_id++) {

    _Bool  is_uniform = false, is_steady = true;

    const cs_param_pty_t  *pty = cs_param_properties + pty_id;

    if (pty->flag & CS_PARAM_FLAG_UNIFORM)
      is_uniform = true;
    if (pty->flag & CS_PARAM_FLAG_UNSTEADY)
      is_steady = false;

    bft_printf(" %s >> uniform [%s], steady [%s], type: %s\n",
               pty->name, cs_base_strtf(is_uniform), cs_base_strtf(is_steady),
               cs_param_pty_type_name[pty->type]);

    switch (pty->def_type) {

    case CS_PARAM_DEF_BY_VALUE:
      {
        const cs_get_t  mat = pty->def.get;

        switch(pty->type) {

        case CS_PARAM_PTY_ISO:
          bft_printf("       value: % 5.3e\n", mat.val);
          break;
        case CS_PARAM_PTY_ORTHO:
          bft_printf("       value: (% 5.3e, % 5.3e, % 5.3e)\n",
                     mat.vect[0], mat.vect[1], mat.vect[2]);
          break;
        case CS_PARAM_PTY_ANISO:
          bft_printf("              |% 5.3e, % 5.3e, % 5.3e|\n",
                     mat.tens[0][0], mat.tens[0][1], mat.tens[0][2]);
          bft_printf("       value: |% 5.3e, % 5.3e, % 5.3e|\n",
                     mat.tens[1][0], mat.tens[1][1], mat.tens[1][2]);
          bft_printf("              |% 5.3e, % 5.3e, % 5.3e|\n",
                     mat.tens[2][0], mat.tens[2][1], mat.tens[2][2]);
          break;
        default:
          bft_error(__FILE__, __LINE__, 0,
                    _(" Invalid type of material property."));
          break;

        } // pty->type
      }
      break; // BY_VALUE

    case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
      bft_printf("         definition by an analytical function\n");
      break;

    case CS_PARAM_DEF_BY_FIELD:
      bft_printf("         definition by a field (id = %d)\n", pty->def.get.id);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of definition for a material property."));
      break;

    } /* switch on def_type */

  } /* Loop on properties */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free structures dedicated to the definition of material properties
 */
/*----------------------------------------------------------------------------*/

void
cs_param_pty_finalize(void)
{
  if (cs_param_properties == NULL)
    return;

  for (int i = 0; i < cs_param_n_properties; i++)
    BFT_FREE(cs_param_properties[i].name);

  BFT_FREE(cs_param_properties);
  cs_param_properties = NULL;
  cs_param_n_properties = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the number of advection fields defined
 *
 * \return the number of advection fields defined
 */
/*----------------------------------------------------------------------------*/

int
cs_param_get_n_adv_fields(void)
{
  return cs_param_n_adv_fields;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query to know if an advection field is uniform
 *
 * \param[in]    adv_id    id related to the advection field to deal with
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_adv_field_is_uniform(int   adv_id)
{
  _check_adv_field_id(adv_id);

  cs_param_adv_field_t  *adv = cs_param_adv_fields + adv_id;

  if (adv->flag & CS_PARAM_FLAG_UNIFORM)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of an advection field from its id
 *
 * \param[in]    adv_id    id related to the advection field to deal with
 *
 * \return  the name of the advection field
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_adv_field_get_name(int    adv_id)
{
  _check_adv_field_id(adv_id);

  cs_param_adv_field_t  *adv = cs_param_adv_fields + adv_id;

  return adv->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a cs_param_adv_field_t structure from its id
 *
 * \param[in]  adv_id   id related to an advection field
 *
 * \return a pointer to a cs_param_adv_field_t
 */
/*----------------------------------------------------------------------------*/

cs_param_adv_field_t *
cs_param_adv_field_get(int  adv_id)
{
  _check_adv_field_id(adv_id);

  return cs_param_adv_fields + adv_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the id related to an advection field definition from its name
 *
 * \param[in]  ref_name    name of the advection field to find
 *
 * \return -1 if not found otherwise the associated id
 */
/*----------------------------------------------------------------------------*/

int
cs_param_adv_get_id_by_name(const char  *ref_name)
{
  int  adv_id = -1;
  int  reflen = strlen(ref_name);

  for (int i = 0; i < cs_param_n_adv_fields; i++) {

    const cs_param_adv_field_t  *adv = cs_param_adv_fields + i;
    int  len = strlen(adv->name);

    if (reflen == len) {
      if (strcmp(ref_name, adv->name) == 0) {
        adv_id = i;
        break;
      }
    }

  } /* Loop on advection fields */

  return adv_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize an advection field structure
 *
 * \param[in]  name        name of the advection field
 * \param[in]  post_freq   -1 (no post-processing), 0 (at the beginning)
 *                         otherwise every post_freq iteration(s)
 */
/*----------------------------------------------------------------------------*/

void
cs_param_adv_field_add(const char      *name,
                       int              post_freq)
{
  int  adv_id = cs_param_adv_get_id_by_name(name);
  cs_param_adv_field_t  *adv = NULL;

  if (adv_id > -1) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_(" An existing advection field has already the same name %s.\n"
                 " Stop adding this advection field.\n"), name);
    return;
  }

  /* Allocate and initialize a new material property */
  BFT_REALLOC(cs_param_adv_fields,
              cs_param_n_adv_fields + 1,
              cs_param_adv_field_t);

  adv = cs_param_adv_fields + cs_param_n_adv_fields;
  cs_param_n_adv_fields++;

  adv->post_freq = post_freq;

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(adv->name, len, char);
  strncpy(adv->name, name, len);

  /*Default initialization */
  adv->flag = 0;
  adv->def_type = CS_PARAM_N_DEF_TYPES;
  for (int k = 0; k < 3; k++)
    adv->def.get.vect[k] = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a way to compute the value of an advection field
 *
 * \param[in]  name      name of the advection field
 * \param[in]  def_key   way of defining the value of the advection field
 * \param[in]  val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_param_adv_field_set(const char   *name,
                       const char   *def_key,
                       const void   *val)
{
  int  adv_id = cs_param_adv_get_id_by_name(name);

  if (adv_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop setting the material property %s.\n"
                " Do not find a similar name in the property database.\n"),
              name);

  cs_param_adv_field_t  *adv = cs_param_adv_fields + adv_id;

  /* Get the type of definition */
  if (strcmp(def_key, "value") == 0) {
    adv->def_type =  CS_PARAM_DEF_BY_VALUE;
    adv->flag |= CS_PARAM_FLAG_UNIFORM;
  }
  else if (strcmp(def_key, "field") == 0)
    adv->def_type = CS_PARAM_DEF_BY_FIELD;
  else if (strcmp(def_key, "evaluator") == 0)
    adv->def_type = CS_PARAM_DEF_BY_EVALUATOR;
  else if (strcmp(def_key, "analytic") == 0)
    adv->def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  else if (strcmp(def_key, "user") == 0)
    adv->def_type = CS_PARAM_DEF_BY_USER_FUNCTION;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the type of definition.\n"
                " Given key: %s\n"
                " Choice among value, field, evaluator, analytic, user, law"
                " or file\n"
                " Please modify your settings."), def_key);

  cs_param_set_def(adv->def_type, CS_PARAM_VAR_VECT, val, &(adv->def));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free structures dedicated to the definition of advection fields
 */
/*----------------------------------------------------------------------------*/

void
cs_param_adv_field_finalize(void)
{
  if (cs_param_adv_fields == NULL)
    return;

  for (int i = 0; i < cs_param_n_adv_fields; i++)
    BFT_FREE(cs_param_adv_fields[i].name);

  BFT_FREE(cs_param_adv_fields);
  cs_param_adv_fields = NULL;
  cs_param_n_adv_fields = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field related to a material property and/or advection
 *         fields
 */
/*----------------------------------------------------------------------------*/

void
cs_param_add_fields(void)
{
  int  id, dim;

  for (id = 0; id < cs_param_n_properties; id++) {

    cs_param_pty_t  *pty = cs_param_properties + id;

    if (pty->post_freq > -1) { // Post-processing is requested
      if (pty->def_type != CS_PARAM_DEF_BY_FIELD) { // not already a field

        _Bool has_previous = (pty->flag & CS_PARAM_FLAG_UNSTEADY) ? true : false;
        int  field_mask = CS_FIELD_PROPERTY;

        /* Define dim */
        switch (pty->type) {
        case CS_PARAM_PTY_ISO:
          dim = 1;
          break;
        case CS_PARAM_PTY_ORTHO:
          dim = 3;
          break;
        case CS_PARAM_PTY_ANISO:
          dim = 9;
          break;
        default:
          dim = 0; // avoid a warning
          bft_error(__FILE__, __LINE__, 0,
                    _(" Type of property for %s is invalid with the creation"
                      " of field.\n"), pty->name);
        }

        cs_field_t  *fld = cs_field_create(pty->name,
                                           field_mask,
                                           CS_MESH_LOCATION_CELLS,
                                           dim,
                                           true,          // interleave
                                           has_previous);

        pty->field_id = cs_field_id_by_name(pty->name);

        /* Allocate and initialize values */
        cs_field_allocate_values(fld);

      }
    }

  } // Loop on material properties

  for (id = 0; id < cs_param_n_adv_fields; id++) {

    cs_param_adv_field_t  *adv = cs_param_adv_fields + id;

    if (adv->post_freq > -1) { // Post-processing is requested
      if (adv->def_type != CS_PARAM_DEF_BY_FIELD) { // not already a field

        _Bool has_previous = (adv->flag & CS_PARAM_FLAG_UNSTEADY) ? true : false;
        int  field_mask = CS_FIELD_PROPERTY;

        cs_field_t  *fld = cs_field_create(adv->name,
                                           field_mask,
                                           CS_MESH_LOCATION_VERTICES,
                                           3,    // always a vector-valued field
                                           true, // interleave
                                           has_previous);

        adv->field_id = cs_field_id_by_name(adv->name);

        /* Allocate and initialize values */
        cs_field_allocate_values(fld);

      }
    }

  } // Loop on advection fields

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_param_bc_t structure
 *
 * \param[in]  default_bc     default boundary condition
 *
 * \return a pointer to the new structure (free with cs_equation_param_t)
 */
/*----------------------------------------------------------------------------*/

cs_param_bc_t *
cs_param_bc_create(cs_param_bc_type_t  default_bc)
{
  cs_param_bc_t  *bc = NULL;

  BFT_MALLOC(bc, 1, cs_param_bc_t);

  bc->default_bc = default_bc;
  /* Initialization by default */
  bc->enforcement = CS_PARAM_BC_ENFORCE_STRONG;
  bc->quad_type = CS_QUADRATURE_BARY;
  bc->use_subdiv = false;

  bc->n_defs = 0;
  bc->defs = NULL;

  return bc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a cs_param_bc_def_t structure
 *
 * \param[in, out] bcpd       pointer to cs_param_bc_def_t struct. to set
 * \param[in]      loc_id     id related to a cs_mesh_location_t
 * \param[in]      bc_type    generic type of admissible boundary conditions
 * \param[in]      var_type   type of variables (scalar, vector, tensor...)
 * \param[in]      def_type   by value, function...
 * \param[in]      coef1      access to the value of the first coef
 * \param[in]      coef2      access to the value of the second coef (optional)
 */
/*----------------------------------------------------------------------------*/

void
cs_param_bc_def_set(cs_param_bc_def_t      *bcpd,
                    int                     loc_id,
                    cs_param_bc_type_t      bc_type,
                    cs_param_var_type_t     var_type,
                    cs_param_def_type_t     def_type,
                    const void             *coef1,
                    const void             *coef2)
{
  if (bcpd == NULL)
    return;

  /* Sanity checks */
  assert(def_type != CS_PARAM_N_DEF_TYPES);
  assert(bc_type != CS_PARAM_N_BC_TYPES);

  bcpd->loc_id = loc_id;
  bcpd->var_type = var_type;
  bcpd->bc_type = bc_type;
  bcpd->def_type = def_type;

  cs_param_set_def(def_type, var_type, coef1, &(bcpd->def_coef1));
  cs_param_set_def(def_type, var_type, coef2, &(bcpd->def_coef2));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new reaction term. The strcuture related to this reaction
 *         term has already been allocated among the list of reaction terms
 *         associated to an equation
 *
 * \param[in, out] rp         pointer to cs_param_reaction_t structure
 * \param[in]      r_name     name of the reaction term
 * \param[in]      pty_id     id of the associated property
 * \param[in]      h_type     type of discrete Hodge op. associated to this term
 * \param[in]      h_algo     algorithm used to build the discrete Hodge op.
 * \param[in]      r_type     type of reaction term
 */
/*----------------------------------------------------------------------------*/

void
cs_param_reaction_term_add(cs_param_reaction_t          *rp,
                           const char                   *r_name,
                           int                           pty_id,
                           cs_param_hodge_type_t         h_type,
                           cs_param_hodge_algo_t         h_algo,
                           cs_param_source_term_type_t   r_type)
{
  if (rp == NULL)
    return;

  rp->type = r_type;
  rp->do_lumping = false;   // No lumping by default

  /* Name of the reaction term */
  int len = strlen(r_name)+1;
  BFT_MALLOC(rp->name, len, char);
  strncpy(rp->name, r_name, len);

  /* Initialiaze the related discrete Hodge operator */
  rp->hodge.pty_id = pty_id;
  rp->hodge.inv_pty = false;     // inverse property ?
  rp->hodge.type = h_type;
  rp->hodge.algo = h_algo;
  rp->hodge.coef = 1/3.;         // Not used by default but set to DGA method
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a source term. This source term is added to the list of
 *         source terms associated to an equation
 *
 * \param[in, out] stp        pointer to cs_param_source_term_t structure
 * \param[in]      st_name    name of the source term
 * \param[in]      ml_id      id of the related to a cs_mesh_location_t struct.
 * \param[in]      type       type of source term
 * \param[in]      var_type   type of variables (scalar, vector, tensor...)
 * \param[in]      quad_type  type of quadrature rule to use
 * \param[in]      def_type   type of definition (by value, function...)
 * \param[in]      val        access to the definition of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_param_source_term_add(cs_param_source_term_t       *stp,
                         const char                   *st_name,
                         int                           ml_id,
                         cs_param_source_term_type_t   type,
                         cs_param_var_type_t           var_type,
                         cs_quadra_type_t              quad_type,
                         cs_param_def_type_t           def_type,
                         const void                   *val)
{
  if (stp == NULL)
    return;

  int len = strlen(st_name)+1;
  BFT_MALLOC(stp->name, len, char);
  strncpy(stp->name, st_name, len);

  stp->ml_id = ml_id;
  stp->type = type;
  stp->var_type = var_type;
  stp->quad_type = quad_type;
  stp->def_type = def_type;
  stp->use_subdiv = false; // default behaviour

  cs_param_set_def(def_type, var_type, val, &(stp->def));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a reaction term
 *
 * \param[in] r_info     cs_param_reaction_t structure
 *
 * \return the name of the reaction term
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_reaction_get_name(const cs_param_reaction_t   r_info)
{
  return r_info.name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of reaction term
 *
 * \param[in] r_info     set of parameters related to a reaction term
 *
 * \return the name associated with this type of reaction term
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_reaction_get_type_name(cs_param_reaction_t  r_info)
{
  switch (r_info.type) {
  case CS_PARAM_REACTION_TYPE_LINEAR:
    return  "Linear";
    break;
  case CS_PARAM_N_REACTION_TYPES:
    return "Not set";
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of reaction term. Stop execution."));
  }

  return "NULL";
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a source term
 *
 * \param[in] st_info     cs_param_source_term_t structure
 *
 * \return the name of the source term
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_source_term_get_name(const cs_param_source_term_t   st_info)
{
  return st_info.name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of type of a given source term structure
 *
 * \param[in] st_info     cs_param_source_term_t structure
 *
 * \return  the name of the type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_source_term_get_type_name(const cs_param_source_term_t   st_info)
{
  return cs_param_source_term_type_name[st_info.type];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a type of variable
 *
 * \param[in] type     cs_param_var_type_t
 *
 * \return the name associated to this type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_var_type_name(const cs_param_var_type_t   type)
{
  return cs_param_var_type_name[type];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a type of definition
 *
 * \param[in] type     cs_param_def_type_t
 *
 * \return the name associated to this type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_def_type_name(const cs_param_def_type_t   type)
{
  return cs_param_def_type_name[type];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of algorithm related to a discrete Hdoge operator
 *
 * \param[in] h_info     cs_param_hodge_t structure
 *
 * \return the name of the algorithm
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_hodge_get_algo_name(const cs_param_hodge_t   h_info)
{
  return cs_param_hodge_algo_desc[h_info.algo];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the type of discrete Hodge operator
 *
 * \param[in] h_info     cs_param_hodge_t structure
 *
 * \return the name of the type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_hodge_get_type_name(const cs_param_hodge_t   h_info)
{
  return cs_param_hodge_type_desc[h_info.type];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the solver
 *
 * \param[in] solver     type of iterative solver
 *
 * \return the associated solver name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_solver_name(cs_param_itsol_type_t  solver)
{
  switch (solver) {
  case CS_PARAM_ITSOL_CG:
    return  "CG";
    break;
  case CS_PARAM_ITSOL_BICG:
    return "BiCGstab";
    break;
  case CS_PARAM_ITSOL_GMRES:
    return "GMRES";
    break;
  case CS_PARAM_ITSOL_AMG:
    return "Algebraic.Multigrid";
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid solver. Stop execution."));
  }

  return "NULL";
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the preconditioner
 *
 * \param[in] precond     type of preconditioner
 *
 * \return the associated preconditioner name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_precond_name(cs_param_precond_type_t  precond)
{
  switch (precond) {
  case CS_PARAM_PRECOND_DIAG:
    return  "Diagonal";
    break;
  case CS_PARAM_PRECOND_POLY1:
    return  "Neumann.Poly.O1";
    break;
  case CS_PARAM_PRECOND_SSOR:
    return  "SSOR";
    break;
  case CS_PARAM_PRECOND_ILU0:
    return  "ILU0";
    break;
  case CS_PARAM_PRECOND_ICC0:
    return  "ICC0";
    break;
  case CS_PARAM_PRECOND_AMG:
    return  "Algebraic.MultiGrid";
    break;
  case CS_PARAM_PRECOND_AS:
    return  "Additive.Schwarz";
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid preconditioner. Stop execution."));
  }

  return "NULL";
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of boundary condition
 *
 * \param[in] bc          type of boundary condition
 *
 * \return the associated bc name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_bc_name(cs_param_bc_type_t  bc)
{
  switch(bc) {

  case CS_PARAM_BC_HMG_DIRICHLET:
    return "Homogeneous Dirichlet";
    break;
  case CS_PARAM_BC_DIRICHLET:
    return "Dirichlet";
    break;
  case CS_PARAM_BC_HMG_NEUMANN:
    return "Homogeneous Neumann";
  case CS_PARAM_BC_NEUMANN:
    return "Neumann";
  case CS_PARAM_BC_ROBIN:
    return "Robin";
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid BC type. Stop execution."));
  }

  return "NULL"; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of enforcement of the boundary condition
 *
 * \param[in] type          type of enforcement of boundary conditions
 *
 * \return the associated name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_bc_enforcement_name(cs_param_bc_enforce_t  type)
{
  switch(type) {

  case CS_PARAM_BC_ENFORCE_STRONG:
    return "strong";
    break;
  case CS_PARAM_BC_ENFORCE_WEAK_PENA:
    return "weak with a big penalization coefficient";
    break;
  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    return "weak using the Nitsche method";
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    return "weak using the symmetrized Nitsche method";
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of enforcement. Stop execution."));
  }

  return "NULL"; // avoid a warning
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
