/*============================================================================
 * Management of the GUI parameters file: boundary conditions
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_selector.h"

#include "cs_ale.h"
#include "cs_atmo.h"
#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_conditions_priv.h"
#include "cs_boundary_zone.h"
#include "cs_cf_thermo.h"
#include "cs_coal_boundary_conditions.h"
#include "cs_combustion_model.h"
#include "cs_equation_param.h"
#include "cs_parameters.h"
#include "cs_gui_util.h"
#include "cs_gui.h"
#include "cs_gui_specific_physics.h"
#include "cs_ht_convert.h"
#include "cs_meg_prototypes.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_physical_model.h"
#include "cs_thermal_model.h"
#include "cs_timer.h"
#include "cs_tree.h"
#include "cs_turbulence_bc.h"
#include "cs_turbulence_model.h"
#include "cs_parall.h"
#include "cs_elec_model.h"
#include "cs_physical_model.h"
#include "cs_vof.h"
#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_boundary_conditions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Enum for boundary conditions */

typedef struct {
  int        read_data;    /* 1 if profile is calculated from data            */
  int        automatic;    /* 1 if nature of the boundary is automatic        */
} cs_meteo_t;

typedef struct {

  int            n_fields; /* number of handled fields */
  int            n_zones;  /* number of associated zones */

  const char   **label;    /* pointer to label for each boundary zone */
  const char   **nature;   /* nature for each boundary zone */
  int           *bc_num;   /* associated number */

  int           *ientfu;   /* 1 for a fuel flow inlet (gas combustion - D3P) */
  int           *ientox;   /* 1 for an air flow inlet (gas combustion - D3P) */
  int           *ientgb;   /* 1 for burned gas inlet (gas combustion) */
  int           *ientgf;   /* 1 for unburned gas inlet (gas combustion) */
  double        *tkent;    /* inlet temperature (gas combustion) */
  double        *fment;    /* Mean Mixture Fraction at Inlet (gas combustion) */
  int           *itype;    /* type of inlet/outlet (compressible model) */
  double        *prein;    /* inlet pressure (compressible model) */
  double        *rhoin;    /* inlet density  (compressible model) */
  double        *tempin;   /* inlet temperature (compressible model) */
  int          **type_code;  /* type of boundary for each variable */
  double        *rough;    /* roughness size */
  bool          *head_loss_e;  /* formula for head loss (free inlet/outlet) */

  cs_meteo_t    *meteo;     /* inlet or outlet info for atmospheric flow */

} cs_gui_boundary_t;

/* xdef contexts associated to various cases
   ----------------------------------------- */

/*! Arguments passed by context pointer using "per zone" values */

typedef struct {

  const  cs_zone_t    *zone;        /*<! Pointer to zone */

  cs_real_t            val;         /*<! Associated value */

} cs_gui_boundary_const_context_t;

/*! Arguments passed by context pointer for velocity inlet */

typedef struct {

  const  cs_zone_t    *zone;        /*<! Pointer to zone */

  int                  norm_type;   /*<!  0: direct value,
                                          1: direct mass flow rate,
                                          2: direct volumic flow rate,
                                         10: MEG for value,
                                         11: MEG for mass flow rate,
                                         12: MEG for volumic flow rate */
  int                  dir_type;    /*<! 0: surface normal,
                                         1: prescribed direction,
                                         2: MEG function */

  cs_real_t            v;           /*<! value when constant */
  cs_real_t            dir_v[3];    /*<! direction value when constant and
                                      prescribed */

  cs_real_t            c_pr;        /*<! imposed pressure (compressible) */
  cs_real_t            c_tk;        /*<! imposed temperature (compressible) */


} cs_gui_boundary_vel_context_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer on the main boundaries structure */

static cs_gui_boundary_t *boundaries = NULL;

static int     _n_b_contexts = 0;
static void  **_b_contexts = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add new constant value-based cs_dof_func_t context info.
 *
 * \param[in]  zone      pointer to associated zone
 * \param[in]  val       associated value
 *
 * \return: pointer to cs_dof_func_t context info
 */
/*----------------------------------------------------------------------------*/

static cs_gui_boundary_const_context_t *
_add_boundary_const_context(const  cs_zone_t   *zone,
                            cs_real_t           val)
{
  BFT_REALLOC(_b_contexts,
              _n_b_contexts+1,
              void *);

  cs_gui_boundary_const_context_t  *c = NULL;
  BFT_MALLOC(c, 1, cs_gui_boundary_const_context_t);

  c->zone = zone;
  c->val = val;

  /* Now set in structure */

  _b_contexts[_n_b_contexts] = c;
  _n_b_contexts += 1;

  return c;
}

/*----------------------------------------------------------------------------
 * Return a node associated with a given zone's boundary condition definition.
 *
 * parameters:
 *   tn <-- first node associated with search
 *
 * return:
 *   pointer to matching child string
 *----------------------------------------------------------------------------*/

static cs_tree_node_t *
_get_zone_bc_node(cs_tree_node_t *tn_start,
                  int             izone)
{
  cs_tree_node_t *tn = tn_start;

  /* if the start BC node is of a different type, search from parent */

  if (strcmp(tn->name, boundaries->nature[izone]))
    tn = cs_tree_node_get_child(tn_start->parent,
                                boundaries->nature[izone]);

  /* Now searh from siblings */

  tn = cs_tree_node_get_sibling_with_tag
        (tn, "label", boundaries->label[izone]);

  return tn;
}

/*-----------------------------------------------------------------------------
 * Get status of data for inlet or outlet information.
 *
 * parameters:
 *   nature      <--  nature of the boundary
 *   label       <--  label of the inlet or the outlet
 *   tag         <--  name of researched data
 *   data       -->   value associated to the data
 *----------------------------------------------------------------------------*/

static void
_boundary_status(const char  *nature,
                 const char  *label,
                 const char  *tag,
                 int         *data)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "boundary_conditions");
  tn = cs_tree_get_node(tn, nature);
  tn = cs_tree_node_get_sibling_with_tag(tn, "label", label);
  tn = cs_tree_get_node(tn, tag);

  cs_gui_node_get_status_int(tn, data);
}

/*-----------------------------------------------------------------------------
 * Get status of data for velocity_pressure node of inlet or outlet information.
 *
 * parameters:
 *   nature      <--  nature of the boundary
 *   label       <--  label of the inlet or the outlet
 *   tag         <--  name of researched data
 *   data       -->   value associated to the data
 *----------------------------------------------------------------------------*/

static void
_boundary_status_vp(const char  *nature,
                    const char  *label,
                    const char  *tag,
                    int         *data)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "boundary_conditions");
  tn = cs_tree_get_node(tn, nature);
  tn = cs_tree_node_get_sibling_with_tag(tn, "label", label);
  tn = cs_tree_get_node(tn, "velocity_pressure");
  tn = cs_tree_get_node(tn, tag);
  cs_gui_node_get_status_int(tn, data);
}

/*-----------------------------------------------------------------------------
 * Check if a zone uses a mapped inlet, and define the associated mapping
 * if this is the case.
 *
 * parameters:
 *   label    <-- label of wall boundary condition
 *    z       <-- pointer to boundary zone
 *----------------------------------------------------------------------------*/

static void
_check_and_add_mapped_inlet(const char       *label,
                            const cs_zone_t  *z)
{
  int mapped_inlet = 0;

  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "boundary_conditions");
  tn = cs_tree_get_node(tn, "inlet");
  tn = cs_tree_node_get_sibling_with_tag(tn, "label", label);

  tn = cs_tree_get_node(tn, "mapped_inlet");
  cs_gui_node_get_status_int(tn, &mapped_inlet);

  if (mapped_inlet) {
    cs_real_t coord_shift[3] = {0., 0., 0.};
    const char *tname[] = {"translation_x",
                           "translation_y",
                           "translation_z"};

    for (int i = 0; i < 3; i++) {
      cs_tree_node_t *node = cs_tree_get_node(tn, tname[i]);

      const  cs_real_t *v = NULL;
      v = cs_tree_node_get_values_real(node);
      if (v != NULL )
        coord_shift[i] = v[0];
    }

    cs_boundary_conditions_add_map(z->location_id,
                                   CS_MESH_LOCATION_CELLS,
                                   coord_shift,
                                   0.1);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_eval_at_location_t function to compute the flow rate at
 *        boundary faces using a MEG generated value.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection.
 *
 * \param[in]       location_id   base associated mesh location id
 * \param[in]       n_elts        number of elements to consider
 * \param[in]       elt_ids       list of elements ids
 * \param[in]       input         pointer to cs_boundary_conditions_inlet_t
 * \param[in, out]  vals          resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_meg_flow_rate(int               location_id,
               cs_lnum_t         n_elts,
               const cs_lnum_t  *elt_ids,
               void             *input,
               void             *vals_p)
{
  assert(location_id == CS_MESH_LOCATION_NONE);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *face_cen = (const cs_real_3_t *)mq->b_face_cog;

  cs_boundary_conditions_open_t  *c
    = (cs_boundary_conditions_open_t *)input;

  cs_real_t  *vals = (cs_real_t *)vals_p;

  /* Prescribed mass flow, depending on options */

  const char *q_meg_formula_type[] = {"flow1_formula",
                                      "flow2_formula"};

  const char *formula_type = NULL;
  if (c->vel_rescale == CS_BC_VEL_RESCALE_MASS_FLOW_RATE)
    formula_type = q_meg_formula_type[0];
  else if (c->vel_rescale == CS_BC_VEL_RESCALE_VOLUME_FLOW_RATE)
    formula_type = q_meg_formula_type[1];

  cs_meg_boundary_function(c->zone->name,
                           n_elts,
                           elt_ids,
                           face_cen,
                           "velocity",
                           formula_type,
                           vals);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_eval_at_location_t function to compute the velocity at
 *        boundary faces using a MEG generated norm and direction.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection.
 *
 * \param[in]       location_id   base associated mesh location id
 * \param[in]       n_elts        number of elements to consider
 * \param[in]       elt_ids       list of elements ids
 * \param[in]       input         pointer to cs_boundary_conditions_inlet_t
 * \param[in, out]  vals          resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_vel_profile(int               location_id,
             cs_lnum_t         n_elts,
             const cs_lnum_t  *elt_ids,
             void             *input,
             void             *vals_p)
{
  assert(   cs_mesh_location_get_type(location_id)
         == CS_MESH_LOCATION_BOUNDARY_FACES);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *f_n = (const cs_real_3_t *)mq->b_face_u_normal;

  cs_boundary_conditions_open_t  *c
    = (cs_boundary_conditions_open_t *)input;

  int  normalization = 0;
  if (c->vel_rescale == CS_BC_VEL_RESCALE_MASS_FLOW_RATE)
    normalization = 1;
  else if (c->vel_rescale == CS_BC_VEL_RESCALE_VOLUME_FLOW_RATE)
    normalization = 2;

  int dir_type = 2;
  if (c->vel_flags & CS_BC_OPEN_UNIFORM_DIRECTION)
    dir_type = 1;
  else if (c->vel_flags & CS_BC_OPEN_NORMAL_DIRECTION)
    dir_type = 0;

  assert(dir_type != 2);  /* Handled by other function */

  cs_real_t  *vals = (cs_real_t *)vals_p;

  cs_real_t v = c->vel_values[3];

  if (normalization > 0)   /* For flow rates, rescaling is done later */
    v = 1.0;

  switch(dir_type) {
  case 0:
    {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
        for (cs_lnum_t k = 0; k < 3; k++)
          vals[i*3 + k] = -f_n[elt_id][k] * v;
      }
    }
    break;
  case 1:
    {
      cs_real_t dir_v[3];
      cs_math_3_normalize(c->vel_values, dir_v);
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        for (cs_lnum_t k = 0; k < 3; k++)
          vals[i*3 + k] = dir_v[k] * v;
      }
    }
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_eval_at_location_t function to compute the velocity at
 *        boundary faces using a MEG generated norm and possibly direction.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection.
 *
 * \param[in]       location_id   base associated mesh location id
 * \param[in]       n_elts        number of elements to consider
 * \param[in]       elt_ids       list of elements ids
 * \param[in]       input         pointer to cs_boundary_conditions_inlet_t
 * \param[in, out]  vals          resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_vel_profile_by_meg_norm(int               location_id,
                         cs_lnum_t         n_elts,
                         const cs_lnum_t  *elt_ids,
                         void             *input,
                         void             *vals_p)
{
  assert(   cs_mesh_location_get_type(location_id)
         == CS_MESH_LOCATION_BOUNDARY_FACES);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *f_n = (const cs_real_3_t *)mq->b_face_u_normal;
  const cs_real_3_t *face_cen = (const cs_real_3_t *)mq->b_face_cog;

  cs_boundary_conditions_open_t  *c
    = (cs_boundary_conditions_open_t *)input;

  int dir_type = 2;
  if (c->vel_flags & CS_BC_OPEN_UNIFORM_DIRECTION)
    dir_type = 1;
  else if (c->vel_flags & CS_BC_OPEN_NORMAL_DIRECTION)
    dir_type = 0;

  cs_real_t  *vals = (cs_real_t *)vals_p;

  /* Local velocity norm */

  cs_real_t *v_loc = NULL;
  BFT_MALLOC(v_loc, n_elts, cs_real_t);
  cs_meg_boundary_function(c->zone->name,
                           n_elts,
                           elt_ids,
                           face_cen,
                           "velocity",
                           "norm_formula",
                           v_loc);

  switch(dir_type) {
  case 0:
    {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
        for (cs_lnum_t k = 0; k < 3; k++)
          vals[i*3 + k] = -f_n[elt_id][k] * v_loc[i];
      }
    }
    break;
  case 1:
    {
      cs_real_t dir_v[3];
      cs_math_3_normalize(c->vel_values, dir_v);
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        for (cs_lnum_t k = 0; k < 3; k++)
          vals[i*3 + k] = dir_v[k] * v_loc[i];
      }
    }
    break;
  case 2:
    {
      cs_real_t *v_dir = NULL;
      BFT_MALLOC(v_dir, 3*n_elts, cs_real_t);
      cs_meg_boundary_function(c->zone->name,
                               n_elts,
                               elt_ids,
                               face_cen,
                               "direction",
                               "formula",
                               v_dir);
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        /* Note: cs_meg_boundary_function output is not interleaved */
        for (cs_lnum_t k = 0; k < 3; k++)
          vals[i*3 + k] = v_dir[k*n_elts + i] * v_loc[i];
      }
      BFT_FREE(v_dir);
    }
    break;
  }

  BFT_FREE(v_loc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_eval_at_location_t function to compute the velocity at
 *        boundary faces using a MEG generated direction and uniform norm.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection.
 *
 * \param[in]       location_id   base associated mesh location id
 * \param[in]       n_elts        number of elements to consider
 * \param[in]       elt_ids       list of elements ids
 * \param[in]       input         pointer to cs_boundary_conditions_inlet_t
 * \param[in, out]  vals          resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_vel_profile_by_meg_dir(int               location_id,
                        cs_lnum_t         n_elts,
                        const cs_lnum_t  *elt_ids,
                        void             *input,
                        void             *vals_p)
{
  assert(   cs_mesh_location_get_type(location_id)
         == CS_MESH_LOCATION_BOUNDARY_FACES);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *face_cen = (const cs_real_3_t *)mq->b_face_cog;

  cs_boundary_conditions_open_t  *c
    = (cs_boundary_conditions_open_t *)input;

  int  normalization = 0;
  if (c->vel_rescale == CS_BC_VEL_RESCALE_MASS_FLOW_RATE)
    normalization = 1;
  else if (c->vel_rescale == CS_BC_VEL_RESCALE_VOLUME_FLOW_RATE)
    normalization = 2;

  cs_real_t  *vals = (cs_real_t *)vals_p;

  cs_real_t *v_dir = NULL;
  BFT_MALLOC(v_dir, 3 * n_elts, cs_real_t);
  cs_meg_boundary_function(c->zone->name,
                           n_elts,
                           elt_ids,
                           face_cen,
                           "direction",
                           "formula",
                           v_dir);

  cs_real_t v = c->vel_values[3];

  if (normalization > 0) {   /* For flow rates, rescaling is done later */
    v = 1.0;
  }

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    /* Note: cs_meg_boundary_function output is not interleaved */
    for (cs_lnum_t k = 0; k < 3; k++)
      vals[i*3 + k] = v_dir[k*n_elts + i] * v;
  }

  BFT_FREE(v_dir);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign cs_dof_func_t function for Dirichlet velocity condition.
 *
 * Builds and adds matching context to list.
 *
 * \param[in]  tn_vp   tree node associated with BC velocity/pressure
 * \param[in]  z       pointer to associated zone
 */
/*----------------------------------------------------------------------------*/

static void
_set_vel_profile(cs_tree_node_t    *tn_vp,
                 const  cs_zone_t  *z)
{
  cs_equation_param_t *eqp = cs_gui_get_equation_param("velocity");

  if (cs_equation_find_bc(eqp, z->name) != NULL)   /* Ignore if already set */
    return;                                        /* (priority) */

  /* Find or add associated context */

  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(z);

  c->bc_pm_zone_num = boundaries->bc_num[z->id - 1];

  /* Initialize context values */

  c->zone = z;

  int norm_type = 0;
  int dir_type = 0;
  cs_real_t  v = 0;
  cs_real_t  dir_v[] = {0, 0, 0};

  c->c_pr = -1;
  c->c_tk = -1;

  const char *choice_v = cs_gui_node_get_tag(tn_vp, "choice");
  const char *choice_d = cs_gui_node_get_tag(tn_vp, "direction");

  /* Set norm from tree values */

  if (cs_gui_strcmp(choice_v, "norm")) {
    norm_type = 0;
    cs_gui_node_get_child_real(tn_vp, choice_v, &v);
  }
  else if (cs_gui_strcmp(choice_v, "flow1")) {
    norm_type = 1;
    cs_gui_node_get_child_real(tn_vp, choice_v, &v);
  }
  else if (cs_gui_strcmp(choice_v, "flow2")) {
    norm_type = 2;
    cs_gui_node_get_child_real(tn_vp, choice_v, &v);
  }
  else if (cs_gui_strcmp(choice_v, "norm_formula")) {
    if (cs_tree_node_get_child_value_str(tn_vp, choice_v) != NULL)
      norm_type = 10;
  }
  else if (cs_gui_strcmp(choice_v, "flow1_formula")) {
    if (cs_tree_node_get_child_value_str(tn_vp, choice_v) != NULL)
      norm_type = 11;
  }
  else if (cs_gui_strcmp(choice_v, "flow2_formula")) {
    if (cs_tree_node_get_child_value_str(tn_vp, choice_v) != NULL)
      norm_type = 12;
  }
  else if (choice_v != NULL)
    bft_error
      (__FILE__, __LINE__, 0,
       _("%s: unexpected \"choice\" type \"%s\" for zone \"%s\"."),
       __func__, choice_v, z->name);

  /* Special cases for compressible model, where values may require
     special scalings */

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
    cs_tree_node_t *tn = cs_tree_get_node(tn_vp, "compressible_type");
    const char *choice_c = cs_gui_node_get_tag(tn, "choice");

  /* For "imposed_inlet" type BC with mass flow, and prescribed pressure
     and temperature, save values to allow computing density.
     The same values are also read for the relevant boundary conditions,
     but duplicated here to avoid inducing dependencies */

    if (cs_gui_strcmp(choice_c, "imposed_inlet")) {
      bool status_p = false, status_t = false;

      cs_tree_node_t *tn_p = cs_tree_node_get_child(tn_vp, "pressure");
      cs_tree_node_t *tn_k = cs_tree_node_get_child(tn_vp, "temperature");
      cs_gui_node_get_status_bool(tn_p, &status_p);
      cs_gui_node_get_status_bool(tn_k, &status_t);

      if (status_p && status_t) {
        cs_gui_node_get_real(tn_p, &(c->c_pr));
        cs_gui_node_get_real(tn_k, &(c->c_tk));
      }
    }

    /* For compressible "subsonic_inlet_PH" type BC, flow is recomputed,
       so only direction counts
       (In the future, assign specific normalization function upstream of
       classical icodcl/rcodcl for this). */

    else if (cs_gui_strcmp(choice_c, "subsonic_inlet_PH")) {
      norm_type = 0;
      v = 1.;
    }
  }

  /* Set direction from tree values */

  else if (cs_gui_strcmp(choice_d, "normal"))
    dir_type = 0;
  else if (cs_gui_strcmp(choice_d, "coordinates")) {
    dir_type = 1;
    cs_gui_node_get_child_real(tn_vp, "direction_x", dir_v);
    cs_gui_node_get_child_real(tn_vp, "direction_y", dir_v+1);
    cs_gui_node_get_child_real(tn_vp, "direction_z", dir_v+2);
  }
  else if (cs_gui_strcmp(choice_d, "formula")) {
    if (cs_tree_node_get_child_value_str(tn_vp, "direction_formula") != NULL)
      dir_type = 2;
  }
  else if (choice_d != NULL)
    bft_error
      (__FILE__, __LINE__, 0,
       _("%s: unexpected \"direction\" type \"%s\" for zone \"%s\"."),
       __func__, choice_d, z->name);

  /* Now associate cs_dof_func_t function for velocity */

  bool eval_set = false;

  if (dir_type == 0) {
    c->vel_values[3] = v;
    if (norm_type == 0 || norm_type == 2) {
      cs_boundary_conditions_open_set_velocity_by_normal_value(z, v);
      eval_set = true;
    }
    else if (norm_type == 1) {
      cs_boundary_conditions_open_set_mass_flow_rate_by_value(z, v);
      eval_set = true;
    }
  }

  else if (dir_type == 1) {
    /* Direction assigned to context in call cases */
    cs_math_3_normalize(dir_v, dir_v);
    if (norm_type == 0) {
      for (int i = 0; i < 3; i++)
        dir_v[i] *= v;
    }

    if (norm_type != 10) {
      cs_boundary_conditions_open_set_velocity_by_value(z, dir_v);
      eval_set = true;
    }
  }

  if (eval_set == false) {
    cs_eval_at_location_t *func = _vel_profile;
    if (norm_type == 10)
      func = _vel_profile_by_meg_norm;
    else if (dir_type == 2)
      func = _vel_profile_by_meg_dir;

    cs_boundary_conditions_open_set_velocity_by_func(z, func, c);

    if (norm_type == 0 && dir_type == 2) {
      assert(func == _vel_profile_by_meg_dir);
      c->vel_values[3] = v;  /* Needed by that function */
    }

    eval_set = true;

    /* Force some flags in case they are not already set
       (as they may be used by functions assigned here) */

    if (dir_type == 0) {
      c->vel_flags |= CS_BC_OPEN_NORMAL_DIRECTION;
    }
    else if (dir_type == 1) {
      c->vel_flags |= CS_BC_OPEN_UNIFORM_DIRECTION;
      for (int i = 0; i < 3; i++)
        c->vel_values[i] = dir_v[i];
    }

  }

  /* Rescaling for flow rate handled by standard mechanism */

  if (norm_type == 1) {
    c->vel_values[3] = v;
    if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1)
      cs_coal_boundary_conditions_inlet_set_air_mass_flow_rate_by_value(z, v);
    else
      cs_boundary_conditions_open_set_mass_flow_rate_by_value(z, v);
  }
  else if (norm_type == 2) {
    c->vel_values[3] = v;
    cs_boundary_conditions_open_set_volume_flow_rate_by_value(z, v);
  }
  else if (norm_type == 11) {
    if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1)
      cs_coal_boundary_conditions_inlet_set_air_mass_flow_rate_by_func
        (z, _meg_flow_rate, c);
    else
      cs_boundary_conditions_open_set_mass_flow_rate_by_func
        (z, _meg_flow_rate, c);
  }
  else if (norm_type == 12) {
    cs_boundary_conditions_open_set_volume_flow_rate_by_func(z,
                                                             _meg_flow_rate,
                                                             c);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute the Enthalpy at boundary faces
 *        using a constant zone temperature.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to cs_gui_boundary_const_context_t
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_const_t2h(cs_lnum_t         n_elts,
               const cs_lnum_t  *elt_ids,
               bool              dense_output,
               void             *input,
               cs_real_t        *retval)
{
  cs_gui_boundary_const_context_t  *c
    = (cs_gui_boundary_const_context_t *)input;

  assert(n_elts == c->zone->n_elts && elt_ids == c->zone->elt_ids);

  cs_real_t t_val = c->val;

  if (dense_output) {
    cs_real_t *t_l;
    BFT_MALLOC(t_l, n_elts, cs_real_t);

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      t_l[i] = t_val;
    }

    cs_ht_convert_t_to_h_faces_z(c->zone, t_l, retval);

    BFT_FREE(t_l);
  }

  else {
    const cs_mesh_t *m = cs_glob_mesh;

    cs_real_t *t_b;
    BFT_MALLOC(t_b, m->n_b_faces, cs_real_t);

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
      t_b[elt_id] = t_val;
    }

    cs_ht_convert_t_to_h_faces_l(n_elts,
                                 elt_ids,
                                 t_b,
                                 retval);
    BFT_FREE(t_b);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute the Enthalpy at boundary faces
 *        using a temperature MEG generated profile.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to cs_gui_boundary_meg_context_t
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_meg_t2h(cs_lnum_t         n_elts,
             const cs_lnum_t  *elt_ids,
             bool              dense_output,
             void             *input,
             cs_real_t        *retval)
{
  cs_gui_boundary_meg_context_t  *c
    = (cs_gui_boundary_meg_context_t *)input;

  assert(strcmp(c->name, "enthalpy") == 0);

  assert(n_elts == c->zone->n_elts && elt_ids == c->zone->elt_ids);

  const cs_real_3_t *face_cen = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;

  cs_real_t *t_loc = NULL;
  BFT_MALLOC(t_loc, n_elts, cs_real_t);
  cs_meg_boundary_function(c->zone->name,
                           n_elts,
                           elt_ids,
                           face_cen,
                           "temperature",
                           c->condition,
                           t_loc);

  if (dense_output)
    cs_ht_convert_t_to_h_faces_z(c->zone, t_loc, retval);

  else {

    const cs_mesh_t *m = cs_glob_mesh;

    cs_real_t *t_b;
    BFT_MALLOC(t_b, m->n_b_faces, cs_real_t);

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
      t_b[elt_id] = t_loc[i];
    }

    cs_ht_convert_t_to_h_faces_l(n_elts,
                                 elt_ids,
                                 t_b,
                                 retval);

    BFT_FREE(t_b);

  }

  BFT_FREE(t_loc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute a rescaled electric potential
 *        at boundary faces using a constant zone value.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to cs_gui_boundary_const_context_t
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_const_elec_rescaled(cs_lnum_t         n_elts,
                         const cs_lnum_t  *elt_ids,
                         bool              dense_output,
                         void             *input,
                         cs_real_t        *retval)
{
  cs_gui_boundary_const_context_t  *c
    = (cs_gui_boundary_const_context_t *)input;

  assert(n_elts == c->zone->n_elts && elt_ids == c->zone->elt_ids);
  assert(cs_glob_elec_option->ielcor == 1);

  cs_real_t pot_val = c->val * cs_glob_elec_option->coejou;

  if (dense_output || elt_ids == NULL) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      retval[i] = pot_val;
    }
  }

  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      retval[elt_ids[i]] = pot_val;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute a rescaled electric potential
 *        at boundary faces using a MEG generated profile.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to cs_gui_boundary_meg_context_t
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_meg_elec_rescaled(cs_lnum_t         n_elts,
                       const cs_lnum_t  *elt_ids,
                       bool              dense_output,
                       void             *input,
                       cs_real_t        *retval)
{
  cs_gui_boundary_meg_context_t  *c
    = (cs_gui_boundary_meg_context_t *)input;

  assert(n_elts == c->zone->n_elts && elt_ids == c->zone->elt_ids);
  assert(cs_glob_physical_model_flag[CS_JOULE_EFFECT] > -1);
  assert(cs_glob_elec_option->ielcor == 1);

  const cs_real_t joule_coef = cs_glob_elec_option->coejou;

  const cs_real_3_t *face_cen
    = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;

  cs_real_t *v_loc = NULL;
  BFT_MALLOC(v_loc, n_elts, cs_real_t);
  cs_meg_boundary_function(c->zone->name,
                           n_elts,
                           elt_ids,
                           face_cen,
                           c->name,
                           c->condition,
                           v_loc);

  if (dense_output || elt_ids == NULL) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      retval[i] = v_loc[i] * joule_coef;
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      retval[elt_ids[i]] = v_loc[i] * joule_coef;
    }
  }

  BFT_FREE(v_loc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute the electic potential using an
 *        implicitely-defined Dirichlet value.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to field structure (cs_field_t *)
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_dirichlet_implicit_elec_pot(cs_lnum_t         n_elts,
                                 const cs_lnum_t  *elt_ids,
                                 bool              dense_output,
                                 void             *input,
                                 cs_real_t        *retval)
{
  const cs_field_t  *f = (const cs_field_t *)input;

  assert (cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > -1);
  assert(f == CS_F_(potr) || f == CS_F_(poti));

  cs_real_t value = cs_glob_elec_option->pot_diff;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
    cs_lnum_t j = dense_output ? i : elt_id;

    retval[j] = value;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute field values at boundary faces
 *        using an implicit Neumann condition (converted to Dirichlet)
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to field structure (cs_field_t *)
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_neumann_implicit(cs_lnum_t         n_elts,
                      const cs_lnum_t  *elt_ids,
                      bool              dense_output,
                      void             *input,
                      cs_real_t        *retval)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const cs_field_t  *f = (const cs_field_t *)input;

  assert(f->location_id == CS_MESH_LOCATION_CELLS);
  const cs_lnum_t dim = f->dim;
  const cs_real_t *val_pre = f->val_pre;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
    cs_lnum_t j = dense_output ? i : elt_id;
    cs_lnum_t c_id = b_face_cells[elt_id];

    for (cs_lnum_t k = 0; k < dim; k++)
      retval[j*dim + k] = val_pre[c_id*dim + k];
  }
}

/*-----------------------------------------------------------------------------
 * Value of velocity for sliding wall.
 *
 * parameters:
 *   tn_vp  <-- tree node associated with BC velocity/pressure
 *   z_name <-- zone name
 *----------------------------------------------------------------------------*/

static void
_sliding_wall(cs_tree_node_t   *tn_vp,
              const char       *z_name)
{
  const char f_name[] = "velocity";
  cs_equation_param_t *eqp = cs_gui_get_equation_param(f_name);

  if (cs_equation_find_bc(eqp, z_name) != NULL)  /* Ignore if already set */
    return;

  cs_field_t  *f = cs_field_by_name(f_name);

  cs_real_t value[3] = {0, 0, 0};

  for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_vp, "dirichlet");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {
    const char *name = cs_gui_node_get_tag(tn, "name");
    int c_id = -1;
    cs_gui_node_get_child_int(tn, "component", &c_id);
    if (strcmp("velocity", name) == 0 && c_id > -1 && c_id < f->dim) {
      const  cs_real_t *v = cs_tree_node_get_values_real(tn);
      if (v != NULL)
        value[c_id] = v[0];
    }
  }

  cs_equation_add_bc_by_value(eqp,
                              CS_PARAM_BC_DIRICHLET,
                              z_name,
                              value);
}

/*-----------------------------------------------------------------------------
 * Values for turbulence variable for the current inlet.
 *
 * parameters:
 *   tn_bc <-- tree node associated with inlet BC
 *   izone <-- associated BC zone id
 *----------------------------------------------------------------------------*/

static void
_inlet_turbulence(cs_tree_node_t  *tn_bc,
                  int              izone)
{
  cs_tree_node_t  *tn_t = cs_tree_node_get_child(tn_bc, "turbulence");
  const char *choice = cs_tree_node_get_tag(tn_t, "choice");
  if (choice == NULL)
    return;

  cs_bc_turbulence_compute_t bc_compute = CS_BC_TURB_NONE;

  if (cs_gui_strcmp(choice, "hydraulic_diameter"))
    bc_compute = CS_BC_TURB_BY_HYDRAULIC_DIAMETER;
  else if(cs_gui_strcmp(choice, "turbulent_intensity"))
    bc_compute = CS_BC_TURB_BY_TURBULENT_INTENSITY;
  else if(cs_gui_strcmp(choice, "formula")) {
    return;
  }
  else
    return;

  const cs_zone_t *z = cs_boundary_zone_by_name(boundaries->label[izone]);

  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(z);

  cs_gui_node_get_child_real(tn_t, "hydraulic_diameter",
                             &(c->hyd_diameter));

  if (cs_gui_strcmp(choice, "turbulent_intensity")) {
    const cs_real_t *v
      = cs_tree_node_get_child_values_real(tn_t, "turbulent_intensity");
    if (v != NULL)
      c->turb_intensity = v[0] * 0.01;
  }

  if (bc_compute == CS_BC_TURB_BY_HYDRAULIC_DIAMETER)
    cs_boundary_conditions_inlet_set_turbulence_hyd_diam(c->zone,
                                                         c->hyd_diameter);

  else if (bc_compute == CS_BC_TURB_BY_TURBULENT_INTENSITY)
    cs_boundary_conditions_inlet_set_turbulence_intensity(c->zone,
                                                          c->turb_intensity);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute a profile at boundary faces
 *        using a MEG generated function for exchange coefficients.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces
 *
 * The retval values can be decomposed into the alpha, u0, and g values in:
 *
 * K du/dn + alpha*(u - u0) = g
 *
 * For multidimentsional cases, we assume scalar alpha, variable dimension u0,
 * and dimension^2 g (g = 0 here but storage required), so we have a stride
 * of 1 + dim + dim*dim, with values alpha, u0, and g in order.
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_meg_exchange_coefficient_profile(cs_lnum_t         n_elts,
                                      const cs_lnum_t  *elt_ids,
                                      bool              dense_output,
                                      void             *input,
                                      cs_real_t        *retval)
{
  cs_gui_boundary_meg_context_t  *c
    = (cs_gui_boundary_meg_context_t *)input;

  const cs_lnum_t dim = c->dim;
  const cs_lnum_t stride = 1 + dim + dim*dim;

  const cs_real_3_t *face_cen
    = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;

  if (dim > 3)
    bft_error(__FILE__, __LINE__, 0,
              _("In %s, variable dimension > 3 (%s) not handled yet\n"
                "for exchange coefficient boundary coefficients."),
                __func__, c->name);

  cs_real_t *v_loc = NULL;
  if (dim == 1)
    BFT_MALLOC(v_loc, 2 * n_elts, cs_real_t);
  else
    BFT_MALLOC(v_loc, dim * n_elts, cs_real_t);

  cs_meg_boundary_function(c->zone->name,
                           n_elts,
                           elt_ids,
                           face_cen,
                           c->name,
                           c->condition,
                           v_loc);

  /* Exchange coefficient first, Dirichlet values second */

  if (dense_output) {  /* common/expected case */

    if (dim == 1) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        retval[i*3]   = -v_loc[n_elts + i];
        retval[i*3+1] = v_loc[i];
        retval[i*3+2] = 0;
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        retval[i*stride] = -v_loc[dim*n_elts + i];
        for (cs_lnum_t k = 0; k < dim; k++) {
          retval[i*stride + k + 1] = v_loc[k*n_elts + i];
        }
        for (cs_lnum_t k = 1+dim; k < stride; k++) {
          retval[i*stride + k] = 0.;
        }
      }
    }

  }
  else { /* sparse/indirect case */

    if (dim == 1) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
        retval[elt_id*3]   = -v_loc[n_elts + i];
        retval[elt_id*3+1] = v_loc[i];
        retval[elt_id*3+2] = 0;
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
        retval[elt_id*stride] = -v_loc[dim*n_elts + i];
        for (cs_lnum_t k = 0; k < dim; k++) {
          retval[elt_id*stride + k + 1] = v_loc[k*n_elts + i];
        }
        for (cs_lnum_t k = 1+dim; k < stride; k++) {
          retval[elt_id*stride + k] = 0.;
        }
      }
    }

  }

  BFT_FREE(v_loc);
}

/*-----------------------------------------------------------------------------
 * Handle specific scalar types for electric potential.
 *
 * parameters:
 *   tn_bc  <-- tree node associated with "scalar" boundary condition
 *   z      <-- associated zone
 *   choice <-- BC type string in tree
 *   f      <-> pointer to associated field
 *   eqp    <-> pointer to associated equation (field) parameters
 *
 * return:
 *   true if a special case was handled, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_boundary_elec_potential(cs_tree_node_t       *tn_s,
                         const cs_zone_t      *z,
                         const char           *choice,
                         cs_field_t           *f,
                         cs_equation_param_t  *eqp)
{
  bool special_case = true;
  int rescale = cs_glob_elec_option->ielcor == 1;

  if (cs_equation_find_bc(eqp, z->name) != NULL)  /* Ignore if already set */
    return special_case;

  cs_param_bc_type_t bc_wall_prescribed
    = (cs_glob_turb_model->iturb > CS_TURB_NONE) ?
      CS_PARAM_BC_WALL_PRESCRIBED : CS_PARAM_BC_DIRICHLET;

  /* BC definition type ? */

  if (strcmp(choice, "dirichlet") == 0) {

    if (rescale) {
      assert(eqp->dim == 1);

      cs_real_t value[1] = {0};
      const cs_real_t *v = cs_tree_node_get_child_values_real(tn_s, choice);
      if (v != NULL) {
        value[0] = *v;
      }

      cs_gui_boundary_const_context_t  *c
        = _add_boundary_const_context(z, value[0]);
      cs_equation_add_bc_by_dof_func(eqp,
                                     bc_wall_prescribed,
                                     z->name,
                                     cs_flag_boundary_face,
                                     _dof_const_elec_rescaled,
                                     c);
    }
    else
      special_case = false;  /* Without rescaling, standard conditions apply */

  }

  else if (rescale && strcmp(choice, "dirichlet_formula") == 0) {

    if (rescale) {
      const char *s = cs_tree_node_get_child_value_str(tn_s, choice);

      if (s != NULL) {
        cs_gui_boundary_meg_context_t  *c
          = cs_gui_boundary_add_meg_context(z, f->name, choice, f->dim);

        cs_equation_add_bc_by_dof_func(eqp,
                                       bc_wall_prescribed,
                                       z->name,
                                       cs_flag_boundary_face,
                                       _dof_meg_elec_rescaled,
                                       c);
      }
    }
    else
      special_case = false;  /* Without rescaling, standard conditions apply */

  }

  else if (cs_gui_strcmp(choice, "dirichlet_implicit")) {

    /* Currently used only real potential (according to theory documentation,
       it seems this should also apply to the imaginary part of potential and
       vector potential ?) */

    assert(f == CS_F_(potr));

    cs_equation_add_bc_by_dof_func(eqp,
                                   bc_wall_prescribed,
                                   z->name,
                                   cs_flag_boundary_face,
                                   _dof_dirichlet_implicit_elec_pot,
                                   f);
  }

  else if (cs_gui_strcmp(choice, "neumann_implicit")) {

    /* Currently used only for vector potential */

    assert(f == CS_F_(potva));

    /* Actually use Diriclet BC with values based on adjacent cell values,
       rather than actual Neumann BC. */

    cs_equation_add_bc_by_dof_func(eqp,
                                   bc_wall_prescribed,
                                   z->name,
                                   cs_flag_boundary_face,
                                   _dof_neumann_implicit,
                                   f);
  }

  else
    special_case = false;

  return special_case;
}

/*-----------------------------------------------------------------------------
 * get scalar's values
 *
 * parameters:
 *   tn_bc <-- tree node associated with boundary condition
 *   z     <-- associated zone
 *   f     <-- associated field
 *----------------------------------------------------------------------------*/

static void
_boundary_scalar(cs_tree_node_t   *tn_bc,
                 const cs_zone_t  *z,
                 cs_field_t       *f)
{
  const int dim = f->dim;

  cs_tree_node_t *tn_s = cs_tree_node_get_child(tn_bc, "scalar");
  tn_s = cs_tree_node_get_sibling_with_tag(tn_s, "name", f->name);

  cs_equation_param_t *eqp = cs_gui_get_equation_param(f->name);

  /* Some specific models or the user may have associated boundary
     conditions already, so if the value here is the default, it should
     probably be ignored to ensure the appropriate priority */

  if (cs_equation_find_bc(eqp, z->name) != NULL)
    return;

  /* BC definition type ? */

  const char *choice = cs_tree_node_get_tag(tn_s, "choice");
  const char *cnv = cs_tree_node_get_tag(tn_s, "convert");

  if (cnv != NULL) {
    if (f != CS_F_(h) || strcmp(cnv, "temperature") != 0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: conversion for field %s from variable %s not handled."),
         __func__, f->name, cnv);
  }

  if (choice == NULL)
    return;

  /* Handle special cases for specific physical models: */

  bool special_case = false;

  if (   cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > -1
      || cs_glob_physical_model_flag[CS_JOULE_EFFECT] > -1) {
    if (f == CS_F_(potr) || f == CS_F_(poti) || f == CS_F_(potva))
      special_case = _boundary_elec_potential(tn_s, z, choice, f, eqp);
  }

  if (special_case)
    return;

  cs_param_bc_type_t bc_wall_prescribed
    = (cs_glob_turb_model->iturb > CS_TURB_NONE) ?
      CS_PARAM_BC_WALL_PRESCRIBED : CS_PARAM_BC_DIRICHLET;

  /* Now handle standard scalar BC types */

  if (   strcmp(choice, "dirichlet") == 0
      || strcmp(choice, "neumann") == 0) {

    /* Read associated values */

    cs_real_t value[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    assert(dim <= 9);

    if (dim == 1) {
      const cs_real_t *v = cs_tree_node_get_child_values_real(tn_s, choice);
      if (v != NULL) {
        value[0] = *v;
      }
    }

    else { /* dim > 1, not produced by the GUI yet */
      for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_s, choice);
           tn != NULL;
           tn = cs_tree_node_get_next_of_name(tn)) {
        int c_id = -1;
        cs_gui_node_get_child_int(tn, "component", &c_id);
        if (c_id > -1 && c_id < f->dim) {
          const  cs_real_t *v = cs_tree_node_get_values_real(tn);
          if (v != NULL)
            value[c_id] = v[0];
        }
      }
    }

    /* Assign definition */

    if (! strcmp(choice, "dirichlet")) {
      if (f == CS_F_(h) && cs_gui_strcmp(cnv, "temperature")) {
        cs_gui_boundary_const_context_t  *c
          = _add_boundary_const_context(z, value[0]);
        cs_equation_add_bc_by_dof_func(eqp,
                                       bc_wall_prescribed,
                                       z->name,
                                       cs_flag_boundary_face,
                                       _dof_const_t2h,
                                       c);
      }

      else
        cs_equation_add_bc_by_value(eqp,
                                    bc_wall_prescribed,
                                    z->name,
                                    value);
    }

    else if (! strcmp(choice, "neumann"))
      cs_equation_add_bc_by_value(eqp,
                                  CS_PARAM_BC_NEUMANN,
                                  z->name,
                                  value);

  }

  else if (! strcmp(choice, "dirichlet_formula")) {

    const char *s = cs_tree_node_get_child_value_str(tn_s, choice);

    if (s != NULL) {
      cs_gui_boundary_meg_context_t  *c
        = cs_gui_boundary_add_meg_context(z, f->name, choice, dim);

      if (f == CS_F_(h) && cs_gui_strcmp(cnv, "temperature"))
        cs_equation_add_bc_by_dof_func(eqp,
                                       bc_wall_prescribed,
                                       z->name,
                                       cs_flag_boundary_face,
                                       _dof_meg_t2h,
                                       c);

      else {
        cs_equation_add_bc_by_dof_func(eqp,
                                       bc_wall_prescribed,
                                       z->name,
                                       cs_flag_boundary_face,
                                       cs_gui_boundary_conditions_dof_func_meg,
                                       c);
      }

    }

  }

  else if (! strcmp(choice, "neumann_formula")) {

    const char *s = cs_tree_node_get_child_value_str(tn_s, choice);
    if (s != NULL) {
      cs_gui_boundary_meg_context_t  *c
        = cs_gui_boundary_add_meg_context(z, f->name, choice, dim);

      cs_equation_add_bc_by_dof_func(eqp,
                                     CS_PARAM_BC_NEUMANN,
                                     z->name,
                                     cs_flag_boundary_face,
                                     cs_gui_boundary_conditions_dof_func_meg,
                                     c);
    }

  }

  else if (! strcmp(choice, "exchange_coefficient_formula")) {

    const char *s = cs_tree_node_get_child_value_str(tn_s, choice);
    if (s != NULL) {
      cs_gui_boundary_meg_context_t  *c
        = cs_gui_boundary_add_meg_context(z, f->name, choice, dim);

      cs_equation_add_bc_by_dof_func(eqp,
                                     CS_PARAM_BC_ROBIN,
                                     z->name,
                                     cs_flag_boundary_face,
                                     _dof_meg_exchange_coefficient_profile,
                                     c);
    }
  }

  else if (! strcmp(choice, "exchange_coefficient")) {

    if (dim > 3)
      bft_error(__FILE__, __LINE__, 0,
                _("In %s, variable dimension > 3 (%s) not handled yet\n"
                  "for exchange coefficient boundary coefficients."),
                __func__, f->name);

    cs_tree_node_t *tn = NULL;
    cs_real_t value[1 + 3 + 9];

    for (int i = 0; i < 1+3+9; i++)
      value[i] = 0.;

    for (tn = cs_tree_node_get_child(tn_s, "dirichlet");
         tn != NULL;
         tn = cs_tree_node_get_next_of_name(tn)) {
      int c_id = (dim == 1) ? 0 : -1;
      cs_gui_node_get_child_int(tn, "component", &c_id);
      if (c_id > -1 && c_id < f->dim) {
        const  cs_real_t *v = cs_tree_node_get_values_real(tn);
        if (v != NULL)
          value[1 + c_id] = v[0];
      }
    }

    tn = cs_tree_node_get_child(tn_s, choice);
    if (tn != NULL) {
      const  cs_real_t *v = cs_tree_node_get_values_real(tn);
      if (v != NULL)
        value[0] = - v[0];
    }

    cs_equation_add_bc_by_value(eqp,
                                CS_PARAM_BC_ROBIN,
                                z->name,
                                value);
  }

  else if (! strcmp(choice, "syrthes_coupling")) {
    /* Handled locally, using legacy BC's for now.
       To switch to an xdef-based definition, it will be necessary
       to store the exchanged values in an array, as a coupling
       can involve more than 1 zone. */
  }

  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: zone %s, BC type %s for variable %s not handled."),
              __func__, z->name, f->name, choice);
}

/*-----------------------------------------------------------------------------
 * Get coal's data for inlet. Check if the current zone is an inlet only
 * for an oxydant, or for oxydant and coal.
 *
 * parameters:
 *   tn_vp <-- tree node associated with velocity and pressure
 *   z     <-- pointer to associated zone
 *----------------------------------------------------------------------------*/

static void
_inlet_coal(cs_tree_node_t   *tn_vp,
            const cs_zone_t  *z)
{
  const int n_coals = cs_glob_combustion_model->coal->n_coals;
  const int *nclpch = cs_glob_combustion_model->coal->n_classes_per_coal;

  cs_coal_bc_inlet_t *ci = cs_coal_boundary_conditions_get_inlet(z);

  cs_gui_node_get_child_real
    (tn_vp, "temperature", &(ci->t_air));
  cs_gui_node_get_child_int
    (tn_vp, "oxydant", &(ci->inmoxy));

  int _n_coals = 0;

  /* Count coal definitions */

  for (cs_tree_node_t *tn0 = cs_tree_node_get_child(tn_vp, "coal");
       tn0 != NULL;
       tn0 = cs_tree_node_get_next_of_name(tn0), _n_coals++) {

    const char *name = cs_tree_node_get_tag(tn0, "name");
    if (name == NULL)
      continue;

    int icoal = -1;
    if (sscanf(name + strlen("coal"), "%d", &icoal) != 1)
      continue;
    icoal -= 1;
    if (icoal +1 > n_coals)
      continue;

    /* mass flow rate of coal */
    const cs_real_t *v = cs_tree_node_get_child_values_real(tn0, "flow1");
    if (v != NULL)
      ci->qimpcp[icoal] = v[0];

    /* temperature of coal */
    v = cs_tree_node_get_child_values_real(tn0, "temperature");
    if (v != NULL) {
      ci->timpcp[icoal] = v[0];
    }

    /* loop on number of class by coal for ratio (%) stored in distch */

    for (int iclass = 0; iclass < nclpch[icoal]; iclass++) {

      char classname[32];
      snprintf(classname, 31, "class%2.2i", iclass+1);

      cs_tree_node_t *tn1
        = cs_tree_get_node_with_tag(tn0, "ratio", "name", classname);
      v = cs_tree_node_get_values_real(tn1);
      if (v != NULL)
        ci->distch[icoal][iclass] = v[0];

    }

  }

  /* if there is no coal, it is an inlet only for oxydant */
  if (_n_coals == 0) {
    ci->ientat = 1;
    ci->ientcp = 0;
  }
  else {
    ci->ientat = 0;
    ci->ientcp = 1;
    if (_n_coals != n_coals)
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid number of coal: %i xml: %i\n"),
                n_coals, _n_coals);
  }
}

/*-----------------------------------------------------------------------------
 * Get gas combustion's data for inlet.
 *
 * Check if the current zone is an inlet only for an oxydant,
 * or for oxydant and coal.
 *
 * parameters:
 *   tn_vp <-- tree node associated with velocity and pressure
 *   izone <-- associated zone id
 *----------------------------------------------------------------------------*/

static void
_inlet_gas(cs_tree_node_t  *tn_vp,
           int              izone)
{
  cs_tree_node_t *tn = cs_tree_get_node(tn_vp, "gas_type");
  const char *choice = cs_gui_node_get_tag(tn, "choice");

  if (cs_gui_strcmp(choice, "oxydant"))
    boundaries->ientox[izone] = 1;

  else if (cs_gui_strcmp(choice, "fuel"))
    boundaries->ientfu[izone] = 1;

  else if (cs_gui_strcmp(choice, "unburned")) {

    boundaries->ientgf[izone] = 1;

    cs_gui_node_get_child_real(tn_vp, "temperature",
                               &boundaries->tkent[izone]);
    cs_gui_node_get_child_real(tn_vp, "fraction",
                               &boundaries->fment[izone]);

  }
  else if (cs_gui_strcmp(choice, "burned")) {

    boundaries->ientgb[izone] = 1;

    cs_gui_node_get_child_real(tn_vp, "temperature",
                               &boundaries->tkent[izone]);
    cs_gui_node_get_child_real(tn_vp, "fraction",
                               &boundaries->fment[izone]);

  }
}

/*-----------------------------------------------------------------------------
 * Get compressible data for inlet.
 *
 * parameters:
 *   tn_vp <-- tree node associated with velocity and pressure
 *   izone <-- associated zone id
 *----------------------------------------------------------------------------*/

static void
_inlet_compressible(cs_tree_node_t  *tn_vp,
                    int              izone)
{
  const cs_zone_t *z = cs_boundary_zone_by_id(izone + 1);

  bool status;

  cs_tree_node_t *tn = cs_tree_get_node(tn_vp, "compressible_type");
  const char *choice = cs_gui_node_get_tag(tn, "choice");

  if (cs_gui_strcmp(choice, "imposed_inlet")) {

    cs_real_t te_in = cs_math_infinite_r;

    boundaries->itype[izone] = CS_ESICF;

    tn = cs_tree_node_get_child(tn_vp, "pressure");
    status = false;
    cs_gui_node_get_status_bool(tn, &status);
    if (status)
      cs_gui_node_get_real(tn, &boundaries->prein[izone]);

    tn = cs_tree_node_get_child(tn_vp, "density");
    status = false;
    cs_gui_node_get_status_bool(tn, &status);
    if (status)
      cs_gui_node_get_real(tn, &boundaries->rhoin[izone]);

    tn = cs_tree_node_get_child(tn_vp, "temperature");
    status = false;
    cs_gui_node_get_status_bool(tn, &status);
    if (status)
      cs_gui_node_get_real(tn, &boundaries->tempin[izone]);

    tn = cs_tree_node_get_child(tn_vp, "energy");
    status = false;
    cs_gui_node_get_status_bool(tn, &status);
    if (status)
      cs_gui_node_get_real(tn, &te_in);

    cs_equation_param_t *eqp = cs_gui_get_equation_param("total_energy");
    cs_equation_remove_bc(eqp, z->name);
    cs_equation_add_bc_by_value(eqp,
                                CS_PARAM_BC_DIRICHLET,
                                z->name,
                                &te_in);

  }
  else if (cs_gui_strcmp(choice, "subsonic_inlet_PH")) {

    boundaries->itype[izone] = CS_EPHCF;

    cs_gui_node_get_child_real
      (tn_vp, "total_pressure", &boundaries->prein[izone]);

    cs_real_t h_in = cs_math_infinite_r;
    cs_gui_node_get_child_real(tn_vp, "enthalpy", &h_in);

    cs_equation_param_t *eqp = cs_gui_get_equation_param("total_energy");
    cs_equation_remove_bc(eqp, z->name);
    cs_equation_add_bc_by_value(eqp,
                                CS_PARAM_BC_DIRICHLET,
                                z->name,
                                &h_in);

  }
}

/*-----------------------------------------------------------------------------
 * Get compressible data for inlet.
 *
 * parameters:
 *   tn_bc <-- tree node associated with boundary zone
 *   izone <-- associated zone id
 *----------------------------------------------------------------------------*/

static void
_outlet_compressible(cs_tree_node_t  *tn_bc,
                     int              izone)
{
  const char *z_name = boundaries->label[izone];

  cs_tree_node_t *tn = cs_tree_node_get_child(tn_bc, "compressible_type");

  const char *choice = cs_tree_node_get_tag(tn, "choice");

  if (cs_gui_strcmp(choice, "supersonic_outlet")) {
    boundaries->itype[izone] = CS_SSPCF;
  }
  else if (cs_gui_strcmp(choice, "subsonic_outlet")) {
    boundaries->itype[izone] = CS_SOPCF;

    const char name[] = "pressure";
    cs_equation_param_t *eqp = cs_gui_get_equation_param(name);

    tn = cs_tree_node_get_child(tn_bc, "dirichlet");
    tn = cs_tree_node_get_sibling_with_tag(tn, "name", name);

    if (tn != NULL) {
      cs_real_t value = 0;
      const  cs_real_t *v = cs_tree_node_get_values_real(tn);
      if (v != NULL)
        value = v[0];

      if (cs_equation_find_bc(eqp, z_name) == NULL)  /* Ignore if already set */
        cs_equation_add_bc_by_value(eqp,
                                    CS_PARAM_BC_DIRICHLET,
                                    z_name,
                                    &value);
    }
  }
}

/*-----------------------------------------------------------------------------
 * Get pressure value for darcy (inlet/outlet/groundwater).
 *
 * parameters:
 *   tn_bc <-- tree node associated with boundary conditions
 *   z     <-- associated zone
 *----------------------------------------------------------------------------*/

static void
_boundary_darcy(cs_tree_node_t   *tn_bc,
                const cs_zone_t  *z)
{
  cs_tree_node_t *tn = cs_tree_node_get_child(tn_bc, "hydraulicHead");
  const char *choice = cs_gui_node_get_tag(tn, "choice");

  tn = cs_tree_node_get_child(tn_bc, choice);
  tn = cs_tree_node_get_sibling_with_tag(tn, "name", "hydraulic_head");

  cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(head));
  if (eqp == NULL)
    eqp = cs_gui_get_equation_param("pressure_head"); /* CDO version */

  if (cs_equation_find_bc(eqp, z->name) != NULL)  /* Ignore if already set */
    return;

  if (cs_gui_strcmp(choice, "dirichlet")) {
    cs_real_t value = 0;
    cs_gui_node_get_real(tn, &value);
    cs_equation_add_bc_by_value(eqp,
                                CS_PARAM_BC_DIRICHLET,
                                z->name,
                                &value);
  }
  else if (cs_gui_strcmp(choice, "neumann")) {
    /* Vector values per component for CDO, scalar (1st component) for legacy */
    cs_real_t value[3] = {0, 0, 0};
    cs_gui_node_get_real(tn, value);
    cs_equation_add_bc_by_value(eqp,
                                CS_PARAM_BC_NEUMANN,
                                z->name,
                                value);
  }
  else if (cs_gui_strcmp(choice, "dirichlet_formula")) {
    if (tn == NULL) { /* compatibility with inconsistant tag */
      tn = cs_tree_node_get_child(tn_bc, choice);
      tn = cs_tree_node_get_sibling_with_tag(tn, "name", "hydraulicHead");
    }
    const char *formula = cs_tree_node_get_child_value_str(tn, "formula");
    if (formula != NULL) {
      cs_gui_boundary_meg_context_t  *c
        = cs_gui_boundary_add_meg_context(z, "hydraulic_head", choice, 1);
      cs_equation_add_bc_by_dof_func(eqp,
                                     CS_PARAM_BC_DIRICHLET,
                                     z->name,
                                     cs_flag_boundary_face,
                                     cs_gui_boundary_conditions_dof_func_meg,
                                     c);
   }
    else {
      bft_printf("Warning : groundwater flow boundary conditions\n"
                 "          without formula for hydraulic head.\n");
    }
  }
}

/*-----------------------------------------------------------------------------
 * Get pressure value for imposed pressure boundary
 *
 * parameters:
 *   tn_bc  <-- tree node associated with boundary conditions
 *   z_name <-- id of the current zone
 *----------------------------------------------------------------------------*/

static void
_boundary_imposed_pressure(cs_tree_node_t  *tn_bc,
                           const char      *z_name)
{
  const char name[] = "pressure";
  cs_tree_node_t *tn = cs_tree_node_get_child(tn_bc, "dirichlet");
  tn = cs_tree_node_get_sibling_with_tag(tn, "name", name);

  cs_real_t value = 0;
  cs_gui_node_get_real(tn, &value);

  cs_equation_param_t *eqp = cs_gui_get_equation_param(name);

  if (cs_equation_find_bc(eqp, z_name) == NULL)  /* Ignore if already set */
    cs_equation_add_bc_by_value(eqp,
                                CS_PARAM_BC_DIRICHLET,
                                z_name,
                                &value);
}

/*-----------------------------------------------------------------------------
 * Helper to get the face list for the izone
 *
 * parameters:
 *   label     <--  boundary label
 *   n_faces   -->  number of faces
 *
 * returns:
 *   pointer to face list
 *----------------------------------------------------------------------------*/

static const cs_lnum_t *
_get_boundary_faces(const char   *label,
                    cs_lnum_t    *n_faces)
{
  const cs_lnum_t *face_ids = NULL;

  const cs_zone_t *z = cs_boundary_zone_by_name(label);

  *n_faces = z->n_elts;
  face_ids = z->elt_ids;

  return face_ids;
}

/*----------------------------------------------------------------------------
 * Boundary conditions treatment: global structure initialization
 *----------------------------------------------------------------------------*/

static void
_init_boundaries(void)
{
  assert(boundaries == NULL);
  int n_fields = cs_field_n_fields();

  int n_zones = cs_tree_get_node_count(cs_glob_tree,
                                       "boundary_conditions/boundary");

  bool solid_fuels = false;
  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1)
    solid_fuels = true;
  bool gas_combustion = false;
  for (cs_physical_model_type_t m_type = CS_COMBUSTION_3PT;
       m_type <= CS_COMBUSTION_LW;
       m_type++) {
    if (cs_glob_physical_model_flag[m_type] > -1)
      gas_combustion = true;
  }
  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1)
    solid_fuels = true;

  cs_boundary_conditions_create_legacy_zone_data();

  BFT_MALLOC(boundaries, 1, cs_gui_boundary_t);

  boundaries->n_fields = n_fields;
  boundaries->n_zones = n_zones;

  BFT_MALLOC(boundaries->label,     n_zones,    const char *);
  BFT_MALLOC(boundaries->nature,    n_zones,    const char *);
  BFT_MALLOC(boundaries->bc_num,    n_zones,    int);

  boundaries->ientfu = NULL;
  boundaries->ientox = NULL;
  boundaries->ientgb = NULL;
  boundaries->ientgf = NULL;

  boundaries->tkent  = NULL;
  boundaries->fment  = NULL;
  boundaries->itype = NULL;
  boundaries->prein = NULL;
  boundaries->rhoin = NULL;
  boundaries->tempin = NULL;

  BFT_MALLOC(boundaries->type_code, n_fields,   int *);

  BFT_MALLOC(boundaries->rough,     n_zones,    double);

  BFT_MALLOC(boundaries->head_loss_e, n_zones,  bool);

  boundaries->meteo = NULL;

  if (gas_combustion) {
    BFT_MALLOC(boundaries->ientfu,  n_zones, int);
    BFT_MALLOC(boundaries->ientox,  n_zones, int);
    BFT_MALLOC(boundaries->ientgb,  n_zones, int);
    BFT_MALLOC(boundaries->ientgf,  n_zones, int);
    BFT_MALLOC(boundaries->tkent,   n_zones, double);
    BFT_MALLOC(boundaries->fment,   n_zones, double);
  }
  else if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
    BFT_MALLOC(boundaries->itype,   n_zones, int);
    BFT_MALLOC(boundaries->prein,   n_zones, double);
    BFT_MALLOC(boundaries->rhoin,   n_zones, double);
    BFT_MALLOC(boundaries->tempin,  n_zones, double);
  }

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1)
    BFT_MALLOC(boundaries->meteo, n_zones, cs_meteo_t);
  else
    boundaries->meteo = NULL;

  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);

    if (f->type & CS_FIELD_VARIABLE) {
      BFT_MALLOC(boundaries->type_code[f->id], n_zones, int);
    }
  }

  /* initialize for each zone */
  for (int izone = 0; izone < n_zones; izone++) {
    boundaries->rough[izone]     = -999;
    boundaries->head_loss_e[izone] = false;

    if (gas_combustion) {
      boundaries->ientfu[izone]  = 0;
      boundaries->ientox[izone]  = 0;
      boundaries->ientgb[izone]  = 0;
      boundaries->ientgf[izone]  = 0;
      boundaries->tkent[izone]   = 0;
      boundaries->fment[izone]   = 0;
    }

    else if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
      boundaries->itype[izone]     = 0;
      boundaries->prein[izone]     = cs_math_infinite_r;
      boundaries->rhoin[izone]     = 0;
      boundaries->tempin[izone]    = cs_math_infinite_r;
    }

    else if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
      boundaries->meteo[izone].read_data = 0;
      boundaries->meteo[izone].automatic = 0;
    }
  }

  /* Initialization of boundary->type_code and boundary->values */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);

    if (f->type & CS_FIELD_VARIABLE) {
      int i = f->id;
      for (int izone = 0; izone < n_zones; izone++) {
        boundaries->type_code[i][izone] = -1;
      }
    }
  }

  /* filling of the "boundaries" structure */

  cs_tree_node_t *tn_b0 = cs_tree_get_node(cs_glob_tree,
                                           "boundary_conditions");

  int izone = 0;

  /* Complete boundary zone definitions */

  for (cs_tree_node_t *tn = cs_tree_get_node(tn_b0, "boundary");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), izone++) {

    /* nature, label and description of the ith boundary zone;
       zones are shifted by 1, as default zone 0 is defined
       first, before GUI-based definitions (and non-GUI-based
       user definitions come last). */

    const char *label = cs_tree_node_get_tag(tn, "label");
    const char *nature = cs_tree_node_get_tag(tn, "nature");

    int bc_num = izone+1;

    const int *vi = cs_tree_node_get_child_values_int(tn, "name");
    if (vi != NULL)
      bc_num = vi[0];

    /* label of the ith initialization zone */

    const cs_zone_t *z = cs_boundary_zone_by_id(izone + 1);

    assert(strcmp(label, z->name) == 0);

    /* Note: boundary nature is determined again later for most cases,
       but at least symmetry is only defined here, and ALE can define
       "wall" sections instead of the appropriate type to handle
       fixed sections, so also predefine the boundary nature here */

    boundaries->label[izone] = z->name;
    boundaries->nature[izone] = nature;
    boundaries->bc_num[izone] = bc_num;

  }

  cs_wall_functions_t *wall_fnt = cs_get_glob_wall_functions();

  /* Now loop on boundary condition definitions proper */

  cs_tree_node_t *tn_b1 = (tn_b0 != NULL) ? tn_b0->children : tn_b0;

  for (cs_tree_node_t *tn = tn_b1; tn != NULL; tn = tn->next) {

    if (cs_gui_strcmp(tn->name, "boundary")) /* handled in previous loop */
      continue;

    const cs_zone_t *z = NULL;

    const char *nature = tn->name;
    const char *label = cs_tree_node_get_tag(tn, "label");
    if (label != NULL)
      z = cs_boundary_zone_by_name_try(label);

    if (z == NULL)  /* may occur when "dead" leaves or present in tree */
      continue;

    izone = z->id - 1;
    assert(izone >= 0);

    /* Note: ALE may define mesh BCs as "wall" zones even where this is not
       appropriate, so skip it here (as this is handled elsewhere) */
    if (strcmp(boundaries->nature[izone], nature))
        continue;

    if (cs_gui_strcmp(nature, "inlet")) {

      _check_and_add_mapped_inlet(label, z);

      cs_tree_node_t *tn_vp
        = cs_tree_node_get_child(tn, "velocity_pressure");

      bool read_inlet_data = false;

      /* Inlet: options for atmospheric flows */
      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
        if (cs_glob_atmo_option->meteo_profile > 0) {
          cs_gui_node_get_child_status_int
            (tn_vp, "meteo_data", &boundaries->meteo[izone].read_data);
          cs_gui_node_get_child_status_int
            (tn_vp, "meteo_automatic", &boundaries->meteo[izone].automatic);
          if (boundaries->meteo[izone].read_data == 1)
            read_inlet_data = true;
        }
      }

      /* Inlet: velocity */
      if (CS_F_(vel) != NULL) {
        if (cs_glob_physical_model_flag[CS_GROUNDWATER] < 0) {
          if (read_inlet_data == false)
            _set_vel_profile(tn_vp, z);
        }
        else
          _boundary_darcy(tn, z);
      }

      /* Inlet: data for coal combustion */
      if (solid_fuels)
        _inlet_coal(tn_vp, z);

      /* Inlet: data for gas combustion */
      if (gas_combustion)
        _inlet_gas(tn_vp, izone);

      /* Inlet: complete data for compressible model */
      if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1)
        _inlet_compressible(tn_vp, izone);

      /* Inlet: turbulence */
      _inlet_turbulence(tn, izone);

    }
    else if (cs_gui_strcmp(nature, "wall")) {

      /* sliding wall: velocity */

      cs_tree_node_t *tn_vp
        = cs_tree_node_get_child(tn, "velocity_pressure");

      if (tn_vp != NULL) {
        /* Wall: ROUGH */
        if (   wall_fnt->iwallf != CS_WALL_F_DISABLED
            && wall_fnt->iwallf != CS_WALL_F_1SCALE_POWER
            && wall_fnt->iwallf != CS_WALL_F_SCALABLE_2SCALES_LOG
            && wall_fnt->iwallf != CS_WALL_F_2SCALES_CONTINUOUS) {
          cs_gui_node_get_child_real(tn_vp, "roughness",
                                     &boundaries->rough[izone]);
        }
      }
    }

    else if (cs_gui_strcmp(nature, "outlet")) {
      /* Outlet: data for ATMOSPHERIC FLOWS */
      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
        _boundary_status_vp("outlet", label, "meteo_data",
                            &boundaries->meteo[izone].read_data);
        _boundary_status_vp("outlet", label, "meteo_automatic",
                            &boundaries->meteo[izone].automatic);
      }

      /* Outlet: data for COMPRESSIBLE MODEL */
      if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1)
        _outlet_compressible(tn, izone);

      /* Inlet: data for darcy */
      if (cs_glob_physical_model_flag[CS_GROUNDWATER] > -1)
        _boundary_darcy(tn, z);
    }

    else if (cs_gui_strcmp(nature, "free_inlet_outlet")) {
      cs_tree_node_t *tn_hlf = cs_tree_get_node(tn, "headLoss/formula");
      const char *hl_formula = cs_tree_node_get_value_str(tn_hlf);
      if (hl_formula != NULL)
        boundaries->head_loss_e[izone] = true;
      else {
        bft_printf("Warning : free inlet outlet boundary conditions\n"
                   "          without external head loss definition\n");
      }
    }
    else if (cs_gui_strcmp(nature, "imposed_p_outlet")) {
      _boundary_imposed_pressure(tn, label);
    }
    else if (cs_gui_strcmp(nature, "groundwater")) {
      _boundary_darcy(tn, z);
    }

    /* for each zone */
    if (!cs_gui_strcmp(nature, "symmetry")) {

      /* Thermal scalar */

      cs_field_t *f_tm = cs_thermal_model_field();

      if (f_tm != NULL) {
        if (boundaries->meteo == NULL)
          _boundary_scalar(tn, z, f_tm);
        else if (boundaries->meteo[izone].read_data == 0)
          _boundary_scalar(tn, z, f_tm);
      }

      const char *scalar_sections[]
        =  {"thermophysical_models/atmospheric_flows/variable",
            "thermophysical_models/joule_effect/variable"};

      /* Meteo scalars only if required */
      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] < 0)
        scalar_sections[0] = NULL;
      else {
        if (boundaries->meteo[izone].read_data != 0)
          scalar_sections[0] = NULL;
      }

      /* Electric arc scalars only if required */
      if (cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] < 0)
        scalar_sections[1] = NULL;

      /* Loop on possible specific model scalar sections */
      for (int s_id = 0; s_id < 2; s_id++) {
        if (scalar_sections[s_id] == NULL)
          continue;
        for (cs_tree_node_t *tn_sv
               = cs_tree_get_node(cs_glob_tree, scalar_sections[s_id]);
             tn_sv != NULL;
             tn_sv = cs_tree_node_get_next_of_name(tn_sv)) {
          const char *_name = cs_gui_node_get_tag(tn_sv, "name");
          cs_field_t *f = cs_field_by_name_try(_name);
          if (f != NULL)
            _boundary_scalar(tn, z, f);
        }
      }

      /* VoF void fraction */
      if (cs_glob_vof_parameters->vof_model > 0) {
        _boundary_scalar(tn, z, CS_F_(void_f));
      }

      /* User scalars */
      for (int f_id = 0; f_id < n_fields; f_id++) {
        cs_field_t  *f = cs_field_by_id(f_id);
        if (   (f->type & CS_FIELD_VARIABLE)
            && (f->type & CS_FIELD_USER))
          _boundary_scalar(tn, z, f);
      }
    }

  }  /* for izones */
}

/*----------------------------------------------------------------------------
 * Boundary conditions treatment: initialize and check zone info
 *
 * parameters:
 *   n_b_faces            <-- number of boundary faces
 *   izfppp               <-- zone number for each boundary face
 *----------------------------------------------------------------------------*/

static void
_init_zones(const cs_lnum_t   n_b_faces,
            int              *izfppp)
{
  const int nozppm = CS_MAX_BC_PM_ZONE_NUM;

  assert(boundaries != NULL);

  int n_zones = cs_tree_get_node_count(cs_glob_tree,
                                       "boundary_conditions/boundary");

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++)
    izfppp[ifac] = 0;

  for (int izone = 0; izone < n_zones; izone++) {

    int zone_nbr = boundaries->bc_num[izone];

    if (zone_nbr > nozppm)
      bft_error(__FILE__, __LINE__, 0,
                _("zone's label number %i is greater than %i,"
                  " the maximum allowed \n"), zone_nbr, nozppm);

    cs_lnum_t n_faces = 0;
    const cs_lnum_t *face_ids
      = _get_boundary_faces(boundaries->label[izone], &n_faces);

    /* check if faces are already marked with a zone number */

    for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {
      cs_lnum_t ifbr = face_ids[f_id];
      izfppp[ifbr] = zone_nbr;
    }

  } /*  for izone */

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Boundary conditions treatment
 *
 * Remember: rcodcl[k][j][i] = rcodcl[k*dim1*dim2 + j*dim1 + i]
 *
 * Fortran Interface:
 *
 * subroutine uiclim
 * *****************
 *
 * integer          itypfb   <-- type of boundary for each face
 *----------------------------------------------------------------------------*/

void cs_gui_boundary_conditions_processing(int *itypfb)
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  const cs_real_3_t *face_cen = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;

  cs_boundary_condition_pm_info_t *bc_pm_info = cs_glob_bc_pm_info;

  bool gas_combustion = false;
  for (cs_physical_model_type_t m_type = CS_COMBUSTION_3PT;
       m_type <= CS_COMBUSTION_LW;
       m_type++) {
    if (cs_glob_physical_model_flag[m_type] > -1)
      gas_combustion = true;
  }

  /* First pass only: initialize izfppp */

  static bool initialized = false;
  if (initialized == false) {
    _init_zones(n_b_faces, bc_pm_info->izfppp);
    initialized = true;
  }

  cs_tree_node_t *tn_bc = cs_tree_get_node(cs_glob_tree,
                                           "boundary_conditions/boundary");

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--boundary zones count: %i\n", boundaries->n_zones);
#endif

  /* At each time-step, loop on boundary face zones:
     set itypfb, rcodcl and icodcl thanks to the arrays
     of the structures defined in the first part of the function */

  for (int izone = 0; izone < boundaries->n_zones; izone++) {

    int zone_nbr = boundaries->bc_num[izone];
    const cs_zone_t *bz = cs_boundary_zone_by_id(zone_nbr);

#if _XML_DEBUG_
    bft_printf("\n---zone %i label: %s\n", zone_nbr, boundaries->label[izone]);
    bft_printf("---zone %i nature: %s\n", zone_nbr, boundaries->nature[izone]);
    bft_printf("---zone %i number of faces: %i\n", zone_nbr, bz->n_elts);
#endif

    /* Boundary conditions by boundary type
       ------------------------------------ */

    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {

      tn_bc = _get_zone_bc_node(tn_bc, izone);

      /* Update the zone's arrays (iqimp, qimp,...)
         because they are re-initialized at each time step
         in PRECLI and PPPRCL routines */

      /* data by zone */

      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1)
        bc_pm_info->iprofm[zone_nbr] = boundaries->meteo[izone].read_data;
      else if (gas_combustion) {
        bc_pm_info->ientfu[zone_nbr] = boundaries->ientfu[izone];
        bc_pm_info->ientox[zone_nbr] = boundaries->ientox[izone];
        bc_pm_info->ientgb[zone_nbr] = boundaries->ientgb[izone];
        bc_pm_info->ientgf[zone_nbr] = boundaries->ientgf[izone];
        bc_pm_info->tkent [zone_nbr] = boundaries->tkent[izone];
        bc_pm_info->fment [zone_nbr] = boundaries->fment[izone];
      }
      else if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {

        if (  boundaries->itype[izone] == CS_ESICF
            ||boundaries->itype[izone] == CS_EPHCF) {
          const cs_field_t  *fp = cs_field_by_name("pressure");
          if (fp->bc_coeffs != NULL) {
            cs_real_t *rcodcl1 = fp->bc_coeffs->rcodcl1;

            for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
              cs_lnum_t face_id = bz->elt_ids[elt_id];
              rcodcl1[face_id] = boundaries->prein[izone];
            }
          }
        }

        if (boundaries->itype[izone] == CS_ESICF) {
          cs_field_t *b_rho = cs_field_by_name_try("boundary_density");
          const cs_field_t  *ft = cs_field_by_name("temperature");
          if (ft->bc_coeffs != NULL) {
            cs_real_t *rcodcl1 = ft->bc_coeffs->rcodcl1;

            for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
              cs_lnum_t face_id = bz->elt_ids[elt_id];
              rcodcl1[face_id] = boundaries->tempin[izone];
              b_rho->val[face_id] = boundaries->rhoin[izone];
            }
          }
        }
      }

      /* data by boundary faces */

      int inlet_type = CS_INLET;

      if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1)
        inlet_type = boundaries->itype[izone];
      else {
        int convective_inlet = 0;
        _boundary_status("inlet", boundaries->label[izone],
                         "convective_inlet", &convective_inlet);
        if (convective_inlet)
          inlet_type = CS_CONVECTIVE_INLET;
      }

      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];

        /* zone number and nature of boundary */
        bc_pm_info->izfppp[face_id] = zone_nbr;
        itypfb[face_id] = inlet_type;
      }

      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
         bc_pm_info->iprofm[zone_nbr] = boundaries->meteo[izone].read_data;
        if (boundaries->meteo[izone].automatic) {
          for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
            cs_lnum_t face_id = bz->elt_ids[elt_id];
             bc_pm_info->iautom[face_id] = 1;
          }
        }
      }

      /* turbulent inlet, with formula */

      cs_turb_model_type_t iturb = CS_TURB_NONE;
      if (cs_glob_turb_model != NULL)
        iturb = cs_glob_turb_model->iturb;

      if (iturb != CS_TURB_NONE) {

        cs_tree_node_t *tn_t = cs_tree_node_get_child(tn_bc, "turbulence");
        const char *formula = cs_tree_node_get_child_value_str(tn_t, "formula");

        if (formula != NULL) {

          const char *model = cs_gui_get_thermophysical_model("turbulence");
          if (model == NULL)
            continue;

          /* Check we have not selected another computation type even if
             MEG formula information is present. */

          cs_boundary_conditions_open_t *c
            = cs_boundary_conditions_open_find(bz);
          if (c != NULL) {
            if (c->turb_compute != CS_BC_TURB_NONE)
              continue;
          }

          int n_bnd_vals = 1;
          if (   cs_gui_strcmp(model, "k-epsilon")
              || cs_gui_strcmp(model, "k-epsilon-PL")
              || cs_gui_strcmp(model, "k-omega-SST")) {
            n_bnd_vals = 2;
          }
          else if (   cs_gui_strcmp(model, "Rij-epsilon")
                   || cs_gui_strcmp(model, "Rij-SSG")) {
            n_bnd_vals = 7;
          }
          else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
            n_bnd_vals = 8;
          }
          else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
            n_bnd_vals = 4;
          }
          else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
            n_bnd_vals = 1;
          }
          else
            bft_error(__FILE__, __LINE__, 0,
                      _("Invalid turbulence model: %s.\n"), model);


          cs_real_t *bnd_vals = NULL;
          BFT_MALLOC(bnd_vals, n_bnd_vals * bz->n_elts, cs_real_t);

          if (   cs_gui_strcmp(model, "k-epsilon")
              || cs_gui_strcmp(model, "k-epsilon-PL")) {

            cs_meg_boundary_function(bz->name,
                                     bz->n_elts,
                                     bz->elt_ids,
                                     face_cen,
                                     "turbulence_ke",
                                     "formula",
                                     bnd_vals);

            cs_field_t *c_k   = cs_field_by_name("k");
            cs_field_t *c_eps = cs_field_by_name("epsilon");

            if (c_k->bc_coeffs != NULL && c_eps->bc_coeffs != NULL) {
              cs_real_t *rcodcl1_k = c_k->bc_coeffs->rcodcl1;
              cs_real_t *rcodcl1_eps = c_eps->bc_coeffs->rcodcl1;

              for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
                cs_lnum_t face_id = bz->elt_ids[elt_id];
                rcodcl1_k[face_id]   = bnd_vals[0 * bz->n_elts + elt_id];
                rcodcl1_eps[face_id] = bnd_vals[1 * bz->n_elts + elt_id];
              }
            }

          }
          else if (   cs_gui_strcmp(model, "Rij-epsilon")
                   || cs_gui_strcmp(model, "Rij-SSG")) {

            cs_meg_boundary_function(bz->name,
                                     bz->n_elts,
                                     bz->elt_ids,
                                     face_cen,
                                     "turbulence_rije",
                                     "formula",
                                     bnd_vals);

            cs_field_t *c_rij = cs_field_by_name("rij");
            cs_field_t *c_eps = cs_field_by_name("epsilon");

            if (c_rij->bc_coeffs != NULL && c_eps->bc_coeffs != NULL) {
              cs_real_t *rcodcl1_rij = c_rij->bc_coeffs->rcodcl1;
              cs_real_t *rcodcl1_eps = c_eps->bc_coeffs->rcodcl1;

              for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
                cs_lnum_t face_id = bz->elt_ids[elt_id];

                /* Values are stored for rij components then epsilon */
                for (cs_lnum_t ii = 0; ii < 6; ii++)
                  rcodcl1_rij[ii*n_b_faces + face_id]
                    = bnd_vals[bz->n_elts*ii + elt_id];

                rcodcl1_eps[face_id] = bnd_vals[bz->n_elts*6 + elt_id];
              }
            }

          }
          else if (cs_gui_strcmp(model, "Rij-EBRSM")) {

            cs_meg_boundary_function(bz->name,
                                     bz->n_elts,
                                     bz->elt_ids,
                                     face_cen,
                                     "turbulence_rij_ebrsm",
                                     "formula",
                                     bnd_vals);

            cs_field_t *c_rij = cs_field_by_name("rij");
            cs_field_t *c_eps = cs_field_by_name("epsilon");
            cs_field_t *c_a   = cs_field_by_name("alpha");

            if (   c_rij->bc_coeffs != NULL
                && c_eps->bc_coeffs != NULL
                && c_a->bc_coeffs != NULL) {
              cs_real_t *rcodcl1_rij = c_rij->bc_coeffs->rcodcl1;
              cs_real_t *rcodcl1_eps = c_eps->bc_coeffs->rcodcl1;
              cs_real_t *rcodcl1_a = c_a->bc_coeffs->rcodcl1;

              for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
                cs_lnum_t face_id = bz->elt_ids[elt_id];

                /* Values are stored for rij components then epsilon and alpha*/
                for (cs_lnum_t ii = 0; ii < 6; ii++)
                  rcodcl1_rij[ii*n_b_faces + face_id]
                    = bnd_vals[bz->n_elts*ii + elt_id];

                rcodcl1_eps[face_id] = bnd_vals[bz->n_elts*6 + elt_id];

                rcodcl1_a[face_id] = bnd_vals[bz->n_elts * 7 + elt_id];
              }
            }

          }
          else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
            cs_meg_boundary_function(bz->name,
                                     bz->n_elts,
                                     bz->elt_ids,
                                     face_cen,
                                     "turbulence_v2f",
                                     "formula",
                                     bnd_vals);

            cs_field_t *c_k   = cs_field_by_name("k");
            cs_field_t *c_eps = cs_field_by_name("epsilon");
            cs_field_t *c_phi = cs_field_by_name("phi");
            cs_field_t *c_a   = cs_field_by_name("alpha");

            if (     c_k->bc_coeffs != NULL && c_eps->bc_coeffs != NULL
                && c_phi->bc_coeffs != NULL && c_a->bc_coeffs != NULL) {
              cs_real_t *rcodcl1_k   = c_k->bc_coeffs->rcodcl1;
              cs_real_t *rcodcl1_eps = c_eps->bc_coeffs->rcodcl1;
              cs_real_t *rcodcl1_phi = c_phi->bc_coeffs->rcodcl1;
              cs_real_t *rcodcl1_a   = c_a->bc_coeffs->rcodcl1;

              for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
                cs_lnum_t face_id = bz->elt_ids[elt_id];
                rcodcl1_k[face_id]   = bnd_vals[0 * bz->n_elts + elt_id];
                rcodcl1_eps[face_id] = bnd_vals[1 * bz->n_elts + elt_id];
                rcodcl1_phi[face_id] = bnd_vals[2 * bz->n_elts + elt_id];
                rcodcl1_a[face_id]   = bnd_vals[3 * bz->n_elts + elt_id];
              }
            }

          }
          else if (cs_gui_strcmp(model, "k-omega-SST")) {
            cs_meg_boundary_function(bz->name,
                                     bz->n_elts,
                                     bz->elt_ids,
                                     face_cen,
                                     "turbulence_kw",
                                     "formula",
                                     bnd_vals);

            cs_field_t *c_k = cs_field_by_name("k");
            cs_field_t *c_o = cs_field_by_name("omega");

            if (c_k->bc_coeffs != NULL && c_o->bc_coeffs != NULL) {
              cs_real_t *rcodcl1_k = c_k->bc_coeffs->rcodcl1;
              cs_real_t *rcodcl1_o = c_o->bc_coeffs->rcodcl1;

              for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
                cs_lnum_t face_id = bz->elt_ids[elt_id];
                rcodcl1_k[face_id] = bnd_vals[0 * bz->n_elts + elt_id];
                rcodcl1_o[face_id] = bnd_vals[1 * bz->n_elts + elt_id];
              }
            }

          }
          else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
            cs_meg_boundary_function(bz->name,
                                     bz->n_elts,
                                     bz->elt_ids,
                                     face_cen,
                                     "turbulence_spalart",
                                     "formula",
                                     bnd_vals);

            cs_field_t *c_nu = cs_field_by_name("nu_tilda");

            if (c_nu->bc_coeffs != NULL) {
              cs_real_t *rcodcl1_nu = c_nu->bc_coeffs->rcodcl1;

              for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
                cs_lnum_t face_id = bz->elt_ids[elt_id];
                rcodcl1_nu[face_id] = bnd_vals[elt_id];
              }
            }

          }
          else
            bft_error(__FILE__, __LINE__, 0,
                      _("Invalid turbulence model: %s.\n"), model);

          BFT_FREE(bnd_vals);
        }
      }

#if _XML_DEBUG_
      if (gas_combustion) {
        bft_printf("-----iqimp=%i \n",
                   bc_pm_info->iqimp[zone_nbr]);
        bft_printf("-----ientox=%i, ientfu=%i, ientgf=%i, ientgb=%i \n",
                   bc_pm_info->ientox[zone_nbr],
                   bc_pm_info->ientfu[zone_nbr],
                   bc_pm_info->ientgf[zone_nbr],
                   bc_pm_info->ientgb[zone_nbr]);
      }
      else if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
        if (boundaries->itype[izone] == CS_ESICF) {
          bft_printf("-----imposed_inlet\n");
          bft_printf("-----premin=%g \n",boundaries->prein[zone_nbr-1]);
          bft_printf("-----rhoin=%g \n",boundaries->rhoin[zone_nbr-1]);
          bft_printf("-----tempin=%g \n",boundaries->tempin[zone_nbr-1]);
        }
        if (boundaries->itype[izone] == CS_EPHCF) {
          bft_printf("-----subsonic_inlet_PH\n");
          bft_printf("-----prein=%g \n",boundaries->prein[zone_nbr-1]);
        }
      }
      else {
        bft_printf("-----iqimp=%i, qimp=%12.5e \n",
                   bc_pm_info->iqimp[zone_nbr],
                   bc_pm_info->qimp[zone_nbr]);
      }

      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
        bft_printf("-----iprofm=%i, automatic=%i \n",
                   bc_pm_info->iprofm[zone_nbr],
                   boundaries->meteo[izone].automatic);
      }
#endif

    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "wall")) {
      int iwall = CS_SMOOTHWALL;

      if (boundaries->rough[izone] >= 0.0) {

        iwall = CS_ROUGHWALL;
        cs_field_t *f_roughness = cs_field_by_name_try("boundary_roughness");
        cs_field_t *f_roughness_t
          = cs_field_by_name_try("boundary_thermal_roughness");

        /* Roughness value (z0) */
        for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
          cs_lnum_t face_id = bz->elt_ids[elt_id];
          assert(f_roughness != NULL);
          f_roughness->val[face_id] = boundaries->rough[izone];

          /* Thermal Roughness value.
           * In this case thermal roughness is equal to the roughness. */
          if (f_roughness_t != NULL)
            f_roughness_t->val[face_id] = boundaries->rough[izone];
        }
      }

      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        bc_pm_info->izfppp[face_id] = zone_nbr;
        itypfb[face_id] = iwall;
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "outlet")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        bc_pm_info->izfppp[face_id] = zone_nbr;
        if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1)
          itypfb[face_id] = boundaries->itype[izone];
        else
          itypfb[face_id] = CS_OUTLET;
      }

      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
         bc_pm_info->iprofm[zone_nbr] = boundaries->meteo[izone].read_data;
        if (boundaries->meteo[izone].automatic) {
          for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
            cs_lnum_t face_id = bz->elt_ids[elt_id];
             bc_pm_info->iautom[face_id] = 1;
          }
        }
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "imposed_p_outlet")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        bc_pm_info->izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_OUTLET;
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "symmetry")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        bc_pm_info->izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_SYMMETRY;
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "free_inlet_outlet")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        bc_pm_info->izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_FREE_INLET;
      }

      if (boundaries->head_loss_e[izone]) {
        cs_real_t *bnd_vals = NULL;
        BFT_MALLOC(bnd_vals, bz->n_elts, cs_real_t);
        cs_meg_boundary_function(bz->name,
                                 bz->n_elts,
                                 bz->elt_ids,
                                 face_cen,
                                 "head_loss",
                                 "formula",
                                 bnd_vals);

        const cs_field_t  *fp = cs_field_by_name("pressure");
        if (fp->bc_coeffs != NULL) {
          cs_real_t *rcodcl2_p   = fp->bc_coeffs->rcodcl2;

          for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
            cs_lnum_t face_id = bz->elt_ids[elt_id];
            rcodcl2_p[face_id] = bnd_vals[elt_id];
          }
        }
        BFT_FREE(bnd_vals);
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "free_surface")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        bc_pm_info->izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_FREE_SURFACE;
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "groundwater")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        bc_pm_info->izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_INDEF;
      }

      /* set velocity to 0 */
      const cs_field_t  *fp2 = cs_field_by_name_try("velocity");
      if (fp2 != NULL) {
        if (fp2->bc_coeffs != NULL) {
          int *icodcl = fp2->bc_coeffs->icodcl;
          cs_real_t *rcodcl1 = fp2->bc_coeffs->rcodcl1;

          for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
            cs_lnum_t face_id = bz->elt_ids[elt_id];
            for (cs_lnum_t i = 0; i < 3; i++) {
              icodcl[i*n_b_faces + face_id] = 3;
              rcodcl1[i*n_b_faces + face_id] = 0.;
            }
          }
        }
      }

    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "undefined")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        bc_pm_info->izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_INDEF;
      }

    }

    else {
      bft_error(__FILE__, __LINE__, 0,
                _("boundary nature %s is unknown \n"),
                boundaries->nature[izone]);
    }

#if _XML_DEBUG_
    if (bz->n_elts > 0) {
      cs_lnum_t face_id = bz->elt_ids[0];

      for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
        const cs_field_t  *f = cs_field_by_id(f_id);
        if (f->type & CS_FIELD_VARIABLE) {
          bft_printf("------%s: icodcl=%i, "
                     "rcodcl1=%12.5e, rcodcl1=%12.5e, rcodcl3=%12.5e\n",
                     f->name,
                     f->bc_coeffs->icodcl[face_id],
                     f->bc_coeffs->rcodcl1[face_id],
                     f->bc_coeffs->rcodcl2[face_id],
                     f->bc_coeffs->rcodcl3[face_id]);
          for (cs_lnum_t j = 1; j < f->dim; j++) {
            bft_printf("rcodcl1_%d=%12.5e, rcodcl2_%d=%12.5e, rcodcl3_%d=%12.5e\n",
                       j, f->bc_coeffs->rcodcl1[face_id],
                       j, f->bc_coeffs->rcodcl2[face_id],
                       j, f->bc_coeffs->rcodcl3[face_id]);
          }
        }
      }
    }
#endif

  } /*  for (izone=0; izone < boundaries->n_zones; izone++) */

  /* Define boundary conditions based on cs_equation_param_t structures */

  cs_boundary_conditions_compute(itypfb);
}

/*----------------------------------------------------------------------------
 * Boundary conditions input verification
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICLVE
 * *****************
 *----------------------------------------------------------------------------*/

//void CS_PROCF (uiclve, UICLVE)(void)
void
cs_gui_boundary_conditions_verify(void)
{
  int inature = -1;

  for (int izone = 0; izone < boundaries->n_zones; izone++) {
    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {
      inature = CS_INLET;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "wall")) {
      inature = CS_ROUGHWALL;
      if (boundaries->rough[izone] < 0.0)
        inature = CS_SMOOTHWALL;
    }
    else if (   cs_gui_strcmp(boundaries->nature[izone], "outlet")
             || cs_gui_strcmp(boundaries->nature[izone], "imposed_p_outlet")) {
      inature = CS_OUTLET;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "symmetry")) {
      inature = CS_SYMMETRY;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "free_inlet_outlet")) {
      inature = CS_FREE_INLET;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "free_surface")
             && cs_glob_ale != CS_ALE_NONE) {
      inature = CS_FREE_SURFACE;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "undefined")) {
      inature = CS_INDEF;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "groundwater")) {
      inature = CS_INDEF;
    }

    if (inature < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("boundary nature %s is unknown \n"),
                boundaries->nature[izone]);

    int zone_nbr = boundaries->bc_num[izone];

    if (zone_nbr > CS_MAX_BC_PM_ZONE_NUM)
      bft_error(__FILE__, __LINE__, 0,
                _("zone's label number %i is greater than %i,"
                  " the maximum allowed \n"), zone_nbr, CS_MAX_BC_PM_ZONE_NUM);

  } /*  for izone */
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define boundary conditions based on setup file.
 *
 * \param[in, out]  bdy   boundaries structure to update
 *                        (if NULL, default to cs_glob_domain->boundaries)
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_define(cs_boundary_t  *bdy)
{
  if (bdy == NULL)
    bdy = cs_glob_domain->boundaries;

  cs_tree_node_t *tn_b0 = cs_tree_get_node(cs_glob_tree,
                                           "boundary_conditions");

  int izone = 0;

  /* Wall function info to filter roughness */
  cs_wall_functions_t *wall_fnt = cs_get_glob_wall_functions();

  /* Build boundary zone definitions */

  for (cs_tree_node_t *tn = cs_tree_get_node(tn_b0, "boundary");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), izone++) {

    /* nature, label and description of the ith boundary zone;
       zones are shifted by 1, as default zone 0 is defined
       first, before GUI-based definitions (and non-GUI-based
       user definitions come last). */

    const char *label = cs_tree_node_get_tag(tn, "label");

    int bc_num = izone+1;

    const int *vi = cs_tree_node_get_child_values_int(tn, "name");
    if (vi != NULL)
      bc_num = vi[0];

    /* label of the ith initialization zone */

    const cs_zone_t *z = cs_boundary_zone_by_id(izone + 1);

    if (strcmp(label, z->name) != 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Mismatch between GUI-defined zone %d (%s)\n"
                  "and boundary condition %d (%s), number %d."),
                z->id, z->name, izone+1, label, bc_num);

    /* Note: boundary nature is determined again later for most cases,
       but at least symmetry is only defined here, and ALE can define
       "wall" sections instead of the appropriate type to handle
       fixed sections, so also predefine the boundary nature here */

    /* Now loop on boundary condition definitions proper */

    cs_tree_node_t *tn_bc = NULL;

    cs_tree_node_t *tn_b1 = (tn_b0 != NULL) ? tn_b0->children : tn_b0;
    for (tn_bc = tn_b1; tn_bc != NULL; tn_bc = tn_bc->next) {

      if (cs_gui_strcmp(tn_bc->name, "boundary")) /* handled in parent loop */
        continue;

      const char *c_label = cs_tree_node_get_tag(tn_bc, "label");
      if (c_label != NULL) {
        /* Search for label matching boundary */
        if (strcmp(c_label, label) == 0)
          break;
      }

    }

    if (tn_bc == NULL)
      continue;

    z = cs_boundary_zone_by_name_try(label);

    if (z == NULL)  /* may occur when "dead" leaves or present in tree */
      continue;

    const char *nature = tn_bc->name;

    cs_boundary_type_t bc_type = 0;

    izone = z->id - 1;
    assert(izone >= 0);

    if (cs_gui_strcmp(nature, "inlet")) {

      bc_type |= CS_BOUNDARY_INLET;

      cs_tree_node_t *tn_vp
        = cs_tree_node_get_child(tn_bc, "velocity_pressure");

      if (cs_glob_physical_model_flag[CS_GROUNDWATER] < 0)
        bc_type |= CS_BOUNDARY_IMPOSED_VEL;

      if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {

        cs_tree_node_t *tnc = cs_tree_get_node(tn_vp, "compressible_type");
        const char *choice = cs_gui_node_get_tag(tnc, "choice");

        if (cs_gui_strcmp(choice, "imposed_inlet"))
          bc_type |= CS_BOUNDARY_INLET_QH;
        else if (cs_gui_strcmp(choice, "subsonic_inlet_PH"))
          bc_type |= CS_BOUNDARY_INLET_SUBSONIC_PH;

      }

    }
    else if (cs_gui_strcmp(nature, "wall")) {

      bc_type |= CS_BOUNDARY_WALL;

      /* sliding wall: velocity */

      cs_tree_node_t *tn_vp
        = cs_tree_node_get_child(tn_bc, "velocity_pressure");

      if (tn_vp != NULL) {
        const char *choice = cs_gui_node_get_tag(tn_vp, "choice");
        if (cs_gui_strcmp(choice, "on")) {
          bc_type |= CS_BOUNDARY_SLIDING_WALL;
          _sliding_wall(tn_vp, label);
        }

        /* check for roughness */
        if (   wall_fnt->iwallf != CS_WALL_F_DISABLED
            && wall_fnt->iwallf != CS_WALL_F_1SCALE_POWER
            && wall_fnt->iwallf != CS_WALL_F_SCALABLE_2SCALES_LOG
            && wall_fnt->iwallf != CS_WALL_F_2SCALES_CONTINUOUS) {
          cs_real_t roughness = -1.;
          cs_gui_node_get_child_real(tn_vp, "roughness", &roughness);
          if (roughness > 0) {
            bc_type |= CS_BOUNDARY_ROUGH_WALL;
            /* Create roughness field if needed */
            cs_field_find_or_create("boundary_roughness",
                                    CS_FIELD_INTENSIVE + CS_FIELD_PROPERTY,
                                    CS_MESH_LOCATION_BOUNDARY_FACES,
                                    1, /* dim */
                                    false); /* has previous */
          }
        }

      }
    }

    else if (cs_gui_strcmp(nature, "outlet")) {

      bc_type |= CS_BOUNDARY_OUTLET;

      if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {

        cs_tree_node_t *tnc = cs_tree_get_node(tn_bc, "compressible_type");
        const char *choice = cs_gui_node_get_tag(tnc, "choice");

        if (cs_gui_strcmp(choice, "supersonic_outlet"))
          bc_type |= CS_BOUNDARY_SUPERSONIC;
        else if (cs_gui_strcmp(choice, "subsonic_outlet"))
          bc_type |= CS_BOUNDARY_SUBSONIC;

      }

    }

    else if (cs_gui_strcmp(nature, "free_inlet_outlet")) {
      bc_type = bc_type | CS_BOUNDARY_INLET | CS_BOUNDARY_OUTLET;
    }

    else if (cs_gui_strcmp(nature, "imposed_p_outlet")) {
      bc_type = bc_type | CS_BOUNDARY_OUTLET;
      bc_type = bc_type | CS_BOUNDARY_IMPOSED_P;
    }

    else if (!cs_gui_strcmp(nature, "symmetry")) {
      bc_type = bc_type | CS_BOUNDARY_SYMMETRY;
    }

    cs_boundary_add(bdy, bc_type, z->name);
  }

  /* Definition of the boundaries structure
     and some equation parameters */

  if (boundaries == NULL)
    _init_boundaries();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free GUI boundary condition structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_free_memory(void)
{
  /* clean memory for global private structure boundaries */

  if (boundaries != NULL) {

    for (int f_id = 0; f_id < boundaries->n_fields; f_id++) {
      const cs_field_t  *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        if (boundaries->type_code != NULL)
          BFT_FREE(boundaries->type_code[f->id]);
      }
    }

    /* Gas combustion */
    BFT_FREE(boundaries->ientfu);
    BFT_FREE(boundaries->ientox);
    BFT_FREE(boundaries->ientgb);
    BFT_FREE(boundaries->ientgf);
    BFT_FREE(boundaries->tkent);
    BFT_FREE(boundaries->fment);

    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
      BFT_FREE(boundaries->itype);
      BFT_FREE(boundaries->prein);
      BFT_FREE(boundaries->rhoin);
      BFT_FREE(boundaries->tempin);
    }
    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1)
      BFT_FREE(boundaries->meteo);

    BFT_FREE(boundaries->label);
    BFT_FREE(boundaries->nature);
    BFT_FREE(boundaries->bc_num);

    BFT_FREE(boundaries->type_code);
    BFT_FREE(boundaries->rough);
    BFT_FREE(boundaries->head_loss_e);

    BFT_FREE(boundaries);
  }

  /* Clean MEG contexts */

  for (int i = 0; i < _n_b_contexts; i++)
    BFT_FREE(_b_contexts[i]);

  BFT_FREE(_b_contexts);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add new MEG-based cs_dof_func_t context info.
 *
 * \param[in]  zone       pointer to associated zone
 * \param[in]  name       name of associated field or array
 * \param[in]  condition  associated condition type
 * \param[in]  dim        associated dimension
 *
 * \return: pointer to cs_dof_func_t context info
 */
/*----------------------------------------------------------------------------*/

cs_gui_boundary_meg_context_t *
cs_gui_boundary_add_meg_context(const  cs_zone_t   *zone,
                                const  char        *name,
                                const  char        *condition,
                                int                 dim)
{
  BFT_REALLOC(_b_contexts,
              _n_b_contexts+1,
              void *);

  cs_gui_boundary_meg_context_t  *c = NULL;
  BFT_MALLOC(c, 1, cs_gui_boundary_meg_context_t);

  c->zone = zone;
  c->name = name;
  c->condition = condition;
  c->dim = dim;

  /* Now set in structure */

  _b_contexts[_n_b_contexts] = c;
  _n_b_contexts += 1;

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute boundary condition values
 *        using a MEG generated function.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_dof_func_meg(cs_lnum_t         n_elts,
                                        const cs_lnum_t  *elt_ids,
                                        bool              dense_output,
                                        void             *input,
                                        cs_real_t        *retval)
{
  cs_gui_boundary_meg_context_t  *c
    = (cs_gui_boundary_meg_context_t *)input;

  const cs_real_3_t *face_cen = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;

  const cs_lnum_t dim = c->dim;

  cs_real_t *v_loc = NULL;
  BFT_MALLOC(v_loc, dim * n_elts, cs_real_t);

  cs_meg_boundary_function(c->zone->name,
                           n_elts,
                           elt_ids,
                           face_cen,
                           c->name,
                           c->condition,
                           v_loc);

  if (dense_output) {  /* common/expected case */

    if (dim == 1) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        retval[i] = v_loc[i];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        for (cs_lnum_t k = 0; k < dim; k++)
          retval[i*dim + k] = v_loc[k*n_elts + i];
      }
    }

  }
  else { /* sparse/indirect case */

    if (dim == 1) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
        retval[elt_id] = v_loc[i];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
        for (cs_lnum_t k = 0; k < dim; k++)
          retval[elt_id*dim + k] = v_loc[k*n_elts + i];
      }
    }

  }

  BFT_FREE(v_loc);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
