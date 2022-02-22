/*============================================================================
 * Management of the GUI parameters file: boundary conditions
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include "cs_boundary_zone.h"
#include "cs_cf_thermo.h"
#include "cs_combustion_model.h"
#include "cs_equation_param.h"
#include "cs_parameters.h"
#include "cs_gui_util.h"
#include "cs_gui.h"
#include "cs_gui_specific_physics.h"
#include "cs_ht_convert.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_default.h"
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
#include "cs_prototypes.h"
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

typedef enum {
  BY_XDEF = -1,           /* to mark usage of newer system */
  HYDRAULIC_DIAMETER,
  TURBULENT_INTENSITY,
} cs_boundary_value_t;

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

  int           *iqimp;    /* 1 if a flow rate is applied */
  int           *ientfu;   /* 1 for a fuel flow inlet (gas combustion - D3P) */
  int           *ientox;   /* 1 for an air flow inlet (gas combustion - D3P) */
  int           *ientgb;   /* 1 for burned gas inlet (gas combustion) */
  int           *ientgf;   /* 1 for unburned gas inlet (gas combustion) */
  int           *ientat;   /* 1 if inlet for oxydant (coal combustion)  */
  int           *ientcp;   /* 1 if inlet for oxydant+coal (coal combustion) */
  int           *icalke;   /* automatic boundaries for turbulent variables */
  double        *qimp;     /* oxydant flow rate (coal combustion) */
  int           *inmoxy;   /* oxydant number (coal combustion) */
  double        *timpat;   /* inlet temperature of oxydant (coal combustion) */
  double        *tkent;    /* inlet temperature (gas combustion) */
  double       **qimpcp;   /* inlet coal flow rate (coal combustion) */
  double       **timpcp;   /* inlet coal temperature (coal combustion)  */
  double        *fment;    /* Mean Mixture Fraction at Inlet (gas combustion) */
  int           *itype;    /* type of inlet/outlet (compressible model) */
  double        *prein;    /* inlet pressure (compressible model) */
  double        *rhoin;    /* inlet density  (compressible model) */
  double        *tempin;   /* inlet temperature (compressible model) */
  double        *dh;       /* inlet hydraulic diameter */
  double        *xintur;   /* inlet turbulent intensity */
  int          **type_code;  /* type of boundary for each variable */
  double      ***distch;   /* ratio for each coal */
  double        *rough;    /* roughness size */
  bool          *head_loss_e;  /* formula for head loss (free inlet/outlet) */

  cs_meteo_t    *meteo;     /* inlet or outlet info for atmospheric flow */

} cs_gui_boundary_t;

/* xdef contexts associated to various cases
   ----------------------------------------- */

/*! Arguments passed by context pointer to cs_meg_* functions */

typedef struct {

  const  cs_zone_t    *zone;        /*<! Pointer to zone */

  const  char         *name;        /*<! Pointer to field or array name */
  const  char         *condition;   /*<! Pointer to condition name type */

  int                  dim;         /*<! Values dimension */

} cs_gui_boundary_meg_context_t;

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
  cs_real_t            c_tk;        /*<! imposed temperature  (compressible) */


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

static cs_gui_boundary_meg_context_t *
_add_boundary_meg_context(const  cs_zone_t   *zone,
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

/*----------------------------------------------------------------------------
 * Return a pointer to equation parameters based on a field or equation name.
 *
 * parameters:
 *   name <-- field or equation name
 *
 * return:
 *   pointer to matching child string
 *----------------------------------------------------------------------------*/

static cs_equation_param_t *
_get_equation_param(const char  *name)
{
  cs_equation_param_t *eqp = NULL;

  cs_field_t *f = cs_field_by_name_try(name);
  if (f != NULL)
    eqp = cs_field_get_equation_param(f);

  /* FIXME: else get by equation name */

  return eqp;
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
 * \brief Convert mass flow to velocity.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 *
 * Note: this function is designed to be called on a "per zone" basis, and
 *       requires MPI collegtives for reduction. In case of many zones,
 *       groupping reductions should be more efficient, but would add
 *       constraints for call grouping. As we generally have many more
 *       MPI_Allreduce calls per time step, this is not expected to have a
 *       significant impact on performance. In the rare cases where tens or
 *       hundreds of zones are considered, using an optimized user-defined
 *       function would be preferable.
 *
 * \param[in]      z             pointer to zone
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      mass_flow     prescribed mass flow
 * \param[in]      bc_pr         prescribed pressure for this zone, or huge
 * \param[in]      bc_tk         prescribed temperature for this zone, or huge
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_b_mass_flow_to_vel(const cs_zone_t  *z,
                    bool              dense_output,
                    cs_real_t         mass_flow,
                    cs_real_t         bc_pr,
                    cs_real_t         bc_tk,
                    cs_real_t        *retval)
{
  if (fabs(mass_flow) <= 0)    /* Nothing to do if mass flow is zero */
    return;

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *f_f_n = (const cs_real_3_t *)mq->b_f_face_normal;
  const cs_real_t *f_f_s = mq->b_f_face_surf;

  cs_real_t s[2] = {0., 0.};

  const cs_field_t *rho_b_f = CS_F_(rho_b);

  const cs_lnum_t n_elts = z->n_elts;
  const cs_lnum_t *elt_ids = z->elt_ids;

  /* Compressible prescribed inlet with given pressure and temperature */

  if (bc_pr < cs_math_infinite_r && bc_tk < cs_math_infinite_r) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t elt_id = elt_ids[i];
      cs_lnum_t j = dense_output ? i : elt_id;

      cs_real_t rho_b = cs_cf_thermo_b_rho_from_pt(elt_id, bc_pr, bc_tk);

      s[0] += f_f_s[elt_id];
      s[1] -= rho_b * cs_math_3_dot_product(retval + 3*j, f_f_n[elt_id]);
    }
  }

  /* Regular case */

  else if (rho_b_f != NULL) {
    const cs_real_t *rho_b = rho_b_f->val;
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t elt_id = elt_ids[i];
      cs_lnum_t j = dense_output ? i : elt_id;

      s[0] += f_f_s[elt_id];
      if (rho_b[elt_id] > 0)
        s[1] -= rho_b[elt_id] * cs_math_3_dot_product(retval + 3*j,
                                                      f_f_n[elt_id]);
      else
        s[1] -= cs_math_3_dot_product(retval + 3*j, f_f_n[elt_id]);
    }
  }

  else
    bft_error
      (__FILE__, __LINE__, 0,
       _("%s: no way to obtain boundary density for mass flow scaling\n"
         " for zone \"%s\"."),
       __func__, z->name);

  cs_parall_sum(2, CS_REAL_TYPE, s);

  /* Return to avoid division by zero if zone is empty */

  if (s[0] <= 0.) {
    bft_error
        (__FILE__, __LINE__, 0,
         _("%s: cannot normalize requested mass flow (%g) for zone \"%s\"\n"
           "  with zero surface)."),
         __func__, fabs(s[0]), z->name);

    return;
  }

  /* Now apply normalization */

  if (fabs(mass_flow) > 1e-30) {

    if (fabs(s[1]) < 1e-30)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: cannot normalize for zone \"%s\"\n"
           "  with zero or quasi-zero initial flow)."),
         __func__, z->name);

    cs_real_t sf = mass_flow / s[1];

    if (dense_output) {
      cs_lnum_t _n_elts = n_elts*3;
      for (cs_lnum_t i = 0; i < _n_elts; i++)
        retval[i] *= sf;
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t j = elt_ids[i];
        for (cs_lnum_t k = 0; k < 3; k++)
          retval[j*3 + k] *= sf;
      }
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute the velocity at boundary faces
 *        using a MEG generated norm and direction.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces * stride
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to cs_gui_boundary_vel_context_t
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_meg_vel_profile(cs_lnum_t         n_elts,
                     const cs_lnum_t  *elt_ids,
                     bool              dense_output,
                     void             *input,
                     cs_real_t        *retval)
{
  static const char *_vel_norm_meg_condition[] = {"norm_formula",
                                                  "flow1_formula",
                                                  "flow2_formula"};

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *f_n = (const cs_real_3_t *)mq->b_f_face_normal;
  const cs_real_t *f_s = mq->b_face_surf;

  cs_gui_boundary_vel_context_t  *c
    = (cs_gui_boundary_vel_context_t *)input;

  int  normalization = c->norm_type % 10;

  assert(normalization < 3);

  /* Prescribed mass flow, depending on options */
  cs_real_t mass_flow = 0;

  cs_real_t *v_dir = NULL;

  if (c->dir_type == 2)
    v_dir = cs_meg_boundary_function(c->zone, "direction", "formula");

  /* Local velocity norm */

  if (c->norm_type == 10) {

    cs_real_t *v_loc
      = cs_meg_boundary_function(c->zone,
                                 "velocity",
                                 _vel_norm_meg_condition[normalization]);

    switch(c->dir_type) {
    case 0:
      {
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
          cs_lnum_t j = dense_output ? i : elt_id;

          for (cs_lnum_t k = 0; k < 3; k++)
            retval[j*3 + k] = -f_n[elt_id][k]/f_s[elt_id] * v_loc[i];
        }
      }
      break;
    case 1:
      {
        cs_real_t dir_v[3];
        cs_math_3_normalize(c->dir_v, dir_v);
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
          cs_lnum_t j = dense_output ? i : elt_id;

          for (cs_lnum_t k = 0; k < 3; k++)
            retval[j*3 + k] = dir_v[k] * v_loc[i];
        }
      }
      break;
    case 2:
      {
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
          cs_lnum_t j = dense_output ? i : elt_id;

          /* Note: cs_meg_boundary_function output is not interleaved */
          for (cs_lnum_t k = 0; k < 3; k++)
            retval[j*3 + k] = v_dir[k*n_elts + i] * v_loc[i];
        }
      }
      break;
    }

    BFT_FREE(v_loc);

  }

  /* Global velocity norm */

  else {

    cs_real_t v = c->v;

    if (c->norm_type > 10) {
      cs_real_t *flow
        = cs_meg_boundary_function(c->zone,
                                   "velocity",
                                   _vel_norm_meg_condition[normalization]);
      v = flow[0];
      BFT_FREE(flow);
    }

    if (normalization == 1) { /* For mass flow, rescaling is done later,
                                 and some models assume an initial norm of 1. */
      mass_flow = v;
      v = 1.0;
    }

    switch(c->dir_type) {
    case 0:
      {
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
          cs_lnum_t j = dense_output ? i : elt_id;

          for (cs_lnum_t k = 0; k < 3; k++)
            retval[j*3 + k] = -f_n[elt_id][k]/f_s[elt_id] * v;
        }
      }
      break;
    case 1:
      {
        cs_real_t dir_v[3];
        cs_math_3_normalize(c->dir_v, dir_v);
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
          cs_lnum_t j = dense_output ? i : elt_id;

          for (cs_lnum_t k = 0; k < 3; k++)
            retval[j*3 + k] = dir_v[k] * v;
        }
      }
      break;
    case 2:
      {
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t elt_id = (elt_ids == NULL) ? i : elt_ids[i];
          cs_lnum_t j = dense_output ? i : elt_id;

          /* Note: cs_meg_boundary_function output is not interleaved */
          for (cs_lnum_t k = 0; k < 3; k++)
            retval[j*3 + k] = v_dir[k*n_elts + i] * v;
        }
      }
      break;
    }

  }

  BFT_FREE(v_dir);

  /* Normalization */

  if (normalization == 2) {
    cs_real_t sf = 1. / c->zone->f_measure;
    if (dense_output || elt_ids == NULL) {
      cs_lnum_t _n_elts = n_elts*3;
      for (cs_lnum_t i = 0; i < _n_elts; i++)
        retval[i] *= sf;
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t j = elt_ids[i];
        for (cs_lnum_t k = 0; k < 3; k++)
          retval[j*3 + k] *= sf;
      }
    }
  }

  else if (normalization == 1) {

    /* In all cases, update total mass flow if given by a MEG function
       (not be needed by all models, but provided in a consistent manner) */

    if (c->norm_type == 11)
      boundaries->qimp[c->zone->id - 1] = mass_flow;

    /* For most models, apply scaling to convert from mass flow
       to velocity; For combustion models, scaling is done later. */

    bool scale_velocity = true;

    for (int i = CS_COMBUSTION_3PT; i <= CS_COMBUSTION_FUEL; i++) {
      if (cs_glob_physical_model_flag[i] >= 0)
        scale_velocity = false;
    }

    if (scale_velocity)
      _b_mass_flow_to_vel(c->zone,
                          dense_output,
                          mass_flow,
                          c->c_pr,
                          c->c_tk,
                          retval);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign cs_dof_func_t function for Dirichlet velocity condition.
 *
 * Builds and adds matching context to list.
 *
 * \param[in]  tn_vp   tree node associated with BC velocity/pressure
 * \param[in]  z       pointer to associated zone
 *
 * \return associated iqimp(izone) value
 */
/*----------------------------------------------------------------------------*/

static void
_set_vel_profile(cs_tree_node_t    *tn_vp,
                 const  cs_zone_t  *z)
{
  cs_equation_param_t *eqp = _get_equation_param("velocity");

  if (cs_equation_find_bc(eqp, z->name) != NULL)   /* Ignore if already set */
    return;                                        /* (priority) */

  /* Add new context in list */

  BFT_REALLOC(_b_contexts,
              _n_b_contexts+1,
              void *);

  cs_gui_boundary_vel_context_t  *c = NULL;
  BFT_MALLOC(c, 1, cs_gui_boundary_vel_context_t);

  _b_contexts[_n_b_contexts] = c;
  _n_b_contexts += 1;

  /* Initialize context values */

  c->zone = z;
  c->norm_type = 0;
  c->dir_type = 0;
  c->v = 0;

  for (int i = 0; i < 3; i++)
    c->dir_v[i] = 0;

  c->c_pr = cs_math_infinite_r;
  c->c_tk = cs_math_infinite_r;

  const char *choice_v = cs_gui_node_get_tag(tn_vp, "choice");
  const char *choice_d = cs_gui_node_get_tag(tn_vp, "direction");

  /* Set norm from tree values */

  if (cs_gui_strcmp(choice_v, "norm")) {
    c->norm_type = 0;
    cs_gui_node_get_child_real(tn_vp, choice_v, &(c->v));
  }
  else if (cs_gui_strcmp(choice_v, "flow1")) {
    c->norm_type = 1;
    cs_gui_node_get_child_real(tn_vp, choice_v, &(c->v));
  }
  else if (cs_gui_strcmp(choice_v, "flow2")) {
    c->norm_type = 2;
    cs_gui_node_get_child_real(tn_vp, choice_v, &(c->v));
  }
  else if (cs_gui_strcmp(choice_v, "norm_formula")) {
    if (cs_tree_node_get_child_value_str(tn_vp, choice_v) != NULL)
      c->norm_type = 10;
  }
  else if (cs_gui_strcmp(choice_v, "flow1_formula")) {
    if (cs_tree_node_get_child_value_str(tn_vp, choice_v) != NULL)
      c->norm_type = 11;
  }
  else if (cs_gui_strcmp(choice_v, "flow2_formula")) {
    if (cs_tree_node_get_child_value_str(tn_vp, choice_v) != NULL)
      c->norm_type = 12;
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
       so only direction counts */

    else if (cs_gui_strcmp(choice_c, "subsonic_inlet_PH")) {
      c->norm_type = 0;
      c->v = 1.;
    }
  }

  /* Set direction from tree values */

  else if (cs_gui_strcmp(choice_d, "normal"))
    c->dir_type = 0;
  else if (cs_gui_strcmp(choice_d, "coordinates")) {
    c->dir_type = 1;
    cs_gui_node_get_child_real(tn_vp, "direction_x", c->dir_v);
    cs_gui_node_get_child_real(tn_vp, "direction_y", c->dir_v+1);
    cs_gui_node_get_child_real(tn_vp, "direction_z", c->dir_v+2);
  }
  else if (cs_gui_strcmp(choice_d, "formula")) {
    if (cs_tree_node_get_child_value_str(tn_vp, "direction_formula") != NULL)
      c->dir_type = 2;
  }
  else if (choice_d != NULL)
    bft_error
      (__FILE__, __LINE__, 0,
       _("%s: unexpected \"direction\" type \"%s\" for zone \"%s\"."),
       __func__, choice_d, z->name);

  /* Now associate cs_dof_funct_t function for velocity */

  cs_equation_add_bc_by_dof_func(eqp,
                                 CS_PARAM_BC_DIRICHLET,
                                 z->name,
                                 cs_flag_boundary_face,  // location flag
                                 _dof_meg_vel_profile,
                                 c);

  /* Also update legacy boundary condition structures */

  if (c->norm_type%10 == 1) {
    int izone = z->id - 1;
    boundaries->iqimp[izone] = 1;
    boundaries->qimp[izone] = c->v;
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

  cs_real_t *t_loc = cs_meg_boundary_function(c->zone,
                                              "temperature",
                                              c->condition);

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

  cs_real_t *v_loc = cs_meg_boundary_function(c->zone,
                                              c->name,
                                              c->condition);

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
  cs_equation_param_t *eqp = _get_equation_param(f_name);

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

  if (cs_gui_strcmp(choice, "hydraulic_diameter"))
    boundaries->icalke[izone] = 1;
  else if(cs_gui_strcmp(choice, "turbulent_intensity"))
    boundaries->icalke[izone] = 2;
  else if(cs_gui_strcmp(choice, "formula")) {
    boundaries->icalke[izone] = 0;
    return;
  }
  else
    return;

  cs_gui_node_get_child_real(tn_t, "hydraulic_diameter",
                             &boundaries->dh[izone]);

  if (cs_gui_strcmp(choice, "turbulent_intensity")) {
    const cs_real_t *v
      = cs_tree_node_get_child_values_real(tn_t, "turbulent_intensity");
    if (v != NULL)
      boundaries->xintur[izone] = v[0] * 0.01;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute a boundary profiles
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

static void
_dof_meg_profile(cs_lnum_t         n_elts,
                 const cs_lnum_t  *elt_ids,
                 bool              dense_output,
                 void             *input,
                 cs_real_t        *retval)
{
  cs_gui_boundary_meg_context_t  *c
    = (cs_gui_boundary_meg_context_t *)input;

  cs_real_t *v_loc = cs_meg_boundary_function(c->zone,
                                              c->name,
                                              c->condition);

  const cs_lnum_t dim = c->dim;

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

  if (dim > 3)
    bft_error(__FILE__, __LINE__, 0,
              _("In %s, variable dimension > 3 (%s) not handled yet\n"
                "for exchange coefficient boundary coefficients."),
                __func__, c->name);

  cs_real_t *v_loc = cs_meg_boundary_function(c->zone,
                                              c->name,
                                              c->condition);

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
                                     CS_PARAM_BC_DIRICHLET,
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
          = _add_boundary_meg_context(z, f->name, choice, f->dim);

        cs_equation_add_bc_by_dof_func(eqp,
                                       CS_PARAM_BC_DIRICHLET,
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
                                   CS_PARAM_BC_DIRICHLET,
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
                                   CS_PARAM_BC_DIRICHLET,
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

  cs_equation_param_t *eqp = _get_equation_param(f->name);

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
                                       CS_PARAM_BC_DIRICHLET,
                                       z->name,
                                       cs_flag_boundary_face,
                                       _dof_const_t2h,
                                       c);
      }

      else
        cs_equation_add_bc_by_value(eqp,
                                    CS_PARAM_BC_DIRICHLET,
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
        = _add_boundary_meg_context(z, f->name, choice, dim);

      if (f == CS_F_(h) && cs_gui_strcmp(cnv, "temperature"))
        cs_equation_add_bc_by_dof_func(eqp,
                                       CS_PARAM_BC_DIRICHLET,
                                       z->name,
                                       cs_flag_boundary_face,
                                       _dof_meg_t2h,
                                       c);

      else {
        cs_equation_add_bc_by_dof_func(eqp,
                                       CS_PARAM_BC_DIRICHLET,
                                       z->name,
                                       cs_flag_boundary_face,
                                       _dof_meg_profile,
                                       c);
      }

    }

  }

  else if (! strcmp(choice, "neumann_formula")) {

    const char *s = cs_tree_node_get_child_value_str(tn_s, choice);
    if (s != NULL) {
      cs_gui_boundary_meg_context_t  *c
        = _add_boundary_meg_context(z, f->name, choice, dim);

      cs_equation_add_bc_by_dof_func(eqp,
                                     CS_PARAM_BC_NEUMANN,
                                     z->name,
                                     cs_flag_boundary_face,
                                     _dof_meg_profile,
                                     c);
    }

  }

  else if (! strcmp(choice, "exchange_coefficient_formula")) {

    const char *s = cs_tree_node_get_child_value_str(tn_s, choice);
    if (s != NULL) {
      cs_gui_boundary_meg_context_t  *c
        = _add_boundary_meg_context(z, f->name, choice, dim);

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
 *   tn_vp    <-- tree node associated with velocity and pressure
 *   izone    <-- associated zone id
 *----------------------------------------------------------------------------*/

static void
_inlet_coal(cs_tree_node_t  *tn_vp,
            int              izone)
{
  const int n_coals = cs_glob_combustion_model->coal.n_coals;
  const int *nclpch = cs_glob_combustion_model->coal.n_classes_per_coal;

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
      boundaries->qimpcp[izone][icoal] = v[0];

    /* temperature of coal */
    v = cs_tree_node_get_child_values_real(tn0, "temperature");
    if (v != NULL)
      boundaries->timpcp[izone][icoal] = v[0];

    /* loop on number of class by coal for ratio (%) stored in distch */

    for (int iclass = 0; iclass < nclpch[icoal]; iclass++) {

      char classname[32];
      snprintf(classname, 31, "class%2.2i", iclass+1);

      cs_tree_node_t *tn1
        = cs_tree_get_node_with_tag(tn0, "ratio", "name", classname);
      v = cs_tree_node_get_values_real(tn1);
      if (v != NULL)
        boundaries->distch[izone][icoal][iclass] = v[0];

    }

  }

  /* if there is no coal, it is an inlet only for oxydant */
  if (_n_coals == 0) {
    boundaries->ientat[izone] = 1;
    boundaries->ientcp[izone] = 0;
  }
  else {
    boundaries->ientat[izone] = 0;
    boundaries->ientcp[izone] = 1;
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

    cs_equation_param_t *eqp = _get_equation_param("total_energy");
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

    cs_equation_param_t *eqp = _get_equation_param("total_energy");
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
    cs_equation_param_t *eqp = _get_equation_param(name);

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
    eqp = _get_equation_param("pressure_head"); /* CDO version */

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
        = _add_boundary_meg_context(z, "hydraulic_head", choice, 1);
      cs_equation_add_bc_by_dof_func(eqp,
                                     CS_PARAM_BC_DIRICHLET,
                                     z->name,
                                     cs_flag_boundary_face,
                                     _dof_meg_profile,
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

  cs_equation_param_t *eqp = _get_equation_param(name);

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
  int icharb, iclass;

  assert(boundaries == NULL);
  int n_fields = cs_field_n_fields();

  int n_zones = cs_tree_get_node_count(cs_glob_tree,
                                       "boundary_conditions/boundary");

  bool solid_fuels = false;
  if (   cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] > -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1)
    solid_fuels = true;
  bool gas_combustion = false;
  for (cs_physical_model_type_t m_type = CS_COMBUSTION_3PT;
       m_type <= CS_COMBUSTION_LW;
       m_type++) {
    if (cs_glob_physical_model_flag[m_type] > -1)
      gas_combustion = true;
  }
  if (   cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] > -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1)
    solid_fuels = true;

  BFT_MALLOC(boundaries, 1, cs_gui_boundary_t);

  boundaries->n_fields = n_fields;
  boundaries->n_zones = n_zones;

  BFT_MALLOC(boundaries->label,     n_zones,    const char *);
  BFT_MALLOC(boundaries->nature,    n_zones,    const char *);
  BFT_MALLOC(boundaries->bc_num,    n_zones,    int);

  BFT_MALLOC(boundaries->iqimp,     n_zones,    int);

  boundaries->ientfu = NULL;
  boundaries->ientox = NULL;
  boundaries->ientgb = NULL;
  boundaries->ientgf = NULL;
  boundaries->ientat = NULL;
  boundaries->ientcp = NULL;

  BFT_MALLOC(boundaries->icalke,    n_zones,    int);
  BFT_MALLOC(boundaries->qimp,      n_zones,    double);

  boundaries->inmoxy = NULL;
  boundaries->timpat = NULL;
  boundaries->tkent  = NULL;
  boundaries->qimpcp = NULL;
  boundaries->timpcp = NULL;
  boundaries->fment  = NULL;
  boundaries->itype = NULL;
  boundaries->prein = NULL;
  boundaries->rhoin = NULL;
  boundaries->tempin = NULL;

  BFT_MALLOC(boundaries->dh,        n_zones,    double);
  BFT_MALLOC(boundaries->xintur,    n_zones,    double);
  BFT_MALLOC(boundaries->type_code, n_fields,   int *);

  boundaries->distch = NULL;

  BFT_MALLOC(boundaries->rough,     n_zones,    double);

  BFT_MALLOC(boundaries->head_loss_e, n_zones,  bool);

  boundaries->meteo = NULL;

  if (solid_fuels) {

    const cs_combustion_model_t *cm = cs_glob_combustion_model;

    BFT_MALLOC(boundaries->ientat, n_zones, int);
    BFT_MALLOC(boundaries->inmoxy, n_zones, int);
    BFT_MALLOC(boundaries->timpat, n_zones, double);
    BFT_MALLOC(boundaries->ientcp, n_zones, int);
    BFT_MALLOC(boundaries->qimpcp, n_zones, double *);
    BFT_MALLOC(boundaries->timpcp, n_zones, double *);
    BFT_MALLOC(boundaries->distch, n_zones, double **);

    for (int izone=0; izone < n_zones; izone++) {
      BFT_MALLOC(boundaries->qimpcp[izone], cm->coal.n_coals, double);
      BFT_MALLOC(boundaries->timpcp[izone], cm->coal.n_coals, double);
      BFT_MALLOC(boundaries->distch[izone], cm->coal.n_coals, double *);

      for (icharb = 0; icharb < cm->coal.n_coals; icharb++)
        BFT_MALLOC(boundaries->distch[izone][icharb],
                   cm->coal.n_classes_per_coal[icharb],
                   double);
    }
  }
  else if (gas_combustion) {
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
    boundaries->iqimp[izone]     = 0;
    boundaries->qimp[izone]      = 0;
    boundaries->icalke[izone]    = 0;
    boundaries->dh[izone]        = 0;
    boundaries->xintur[izone]    = 0;
    boundaries->rough[izone]     = -999;
    boundaries->head_loss_e[izone] = false;

    if (solid_fuels) {
      const cs_combustion_model_t *cm = cs_glob_combustion_model;

      boundaries->ientat[izone] = 0;
      boundaries->inmoxy[izone] = 1;
      boundaries->ientcp[izone] = 0;
      boundaries->timpat[izone] = 0;

      for (icharb = 0; icharb < cm->coal.n_coals; icharb++) {
        boundaries->qimpcp[izone][icharb] = 0;
        boundaries->timpcp[izone][icharb] = 0;

        for (iclass = 0; iclass < cm->coal.n_classes_per_coal[icharb]; iclass++)
          boundaries->distch[izone][icharb][iclass] = 0;
      }
    }

    else if (gas_combustion) {
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
      if (cs_glob_physical_model_flag[CS_GROUNDWATER] < 0) {
        if (read_inlet_data == false)
          _set_vel_profile(tn_vp, z);
      }
      else
        _boundary_darcy(tn, z);

      /* Inlet: data for coal combustion */
      if (solid_fuels) {
        cs_gui_node_get_child_real
          (tn_vp, "temperature", &boundaries->timpat[izone]);
        cs_gui_node_get_child_int
          (tn_vp, "oxydant", &boundaries->inmoxy[izone]);
        _inlet_coal(tn_vp, izone);
      }

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
 *   nozppm               <-- max number of boundary conditions zone
 *   izfppp               <-- zone number for each boundary face
 *----------------------------------------------------------------------------*/

static void
_init_zones(const cs_lnum_t   n_b_faces,
            const int        *nozppm,
            int              *izfppp)
{

  assert(boundaries != NULL);

  int n_zones = cs_tree_get_node_count(cs_glob_tree,
                                       "boundary_conditions/boundary");

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++)
    izfppp[ifac] = 0;

  for (int izone = 0; izone < n_zones; izone++) {

    int zone_nbr = boundaries->bc_num[izone];

    if (nozppm && zone_nbr > *nozppm)
      bft_error(__FILE__, __LINE__, 0,
                _("zone's label number %i is greater than %i,"
                  " the maximum allowed \n"), zone_nbr, *nozppm);

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

/*----------------------------------------------------------------------------
 * Define automatic turbulence values.
 *
 * parameters:
 *   rdcodcl <-> boundary condition values
 *----------------------------------------------------------------------------*/

static void
_standard_turbulence_bcs(cs_real_t  rcodcl[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const int var_key_id = cs_field_key_id("variable_id");

  const cs_real_t *b_rho = CS_F_(rho_b)->val;
  const cs_real_t *c_mu = CS_F_(mu)->val;

  cs_lnum_t ivar_u = cs_field_get_key_int(CS_F_(vel), var_key_id) -1;

  const cs_real_t  *_rcodcl_v[3];
  for (int coo_id = 0; coo_id < 3; coo_id++)
    _rcodcl_v[coo_id] = rcodcl + (ivar_u + coo_id)*n_b_faces;

  /* Inlet BC's */

  for (int izone = 0; izone < boundaries->n_zones; izone++) {

    if (boundaries->icalke[izone] > 0) {

      int zone_nbr = boundaries->bc_num[izone];
      const cs_zone_t *z = cs_boundary_zone_by_id(zone_nbr);

      const cs_real_t dh = boundaries->dh[izone];

      if (boundaries->icalke[izone] == 1) {

        for (cs_lnum_t i = 0; i < z->n_elts; i++) {
          cs_lnum_t face_id = z->elt_ids[i];
          cs_lnum_t cell_id = b_face_cells[face_id];

          cs_real_t vel[3] = {_rcodcl_v[0][face_id],
                              _rcodcl_v[1][face_id],
                              _rcodcl_v[2][face_id]};
          cs_real_t uref2 = fmax(cs_math_3_square_norm(vel),
                                 cs_math_epzero);

          cs_turbulence_bc_inlet_hyd_diam(face_id, uref2, dh,
                                          b_rho[face_id], c_mu[cell_id],
                                          rcodcl);
        }

      }
      else if (boundaries->icalke[izone] == 2) {

        const cs_real_t t_intensity = boundaries->xintur[izone];

        for (cs_lnum_t i = 0; i < z->n_elts; i++) {
          cs_lnum_t face_id = z->elt_ids[i];

          cs_real_t vel[3] = {_rcodcl_v[0][face_id],
                              _rcodcl_v[1][face_id],
                              _rcodcl_v[2][face_id]};
          cs_real_t uref2 = fmax(cs_math_3_square_norm(vel),
                                 cs_math_epzero);

          cs_turbulence_bc_inlet_turb_intensity(face_id, uref2,
                                                t_intensity, dh,
                                                rcodcl);
        }

      }

    } /* End for inlet */

  }

  /* Automatic wall condition for alpha */

  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {

    cs_lnum_t ivar_alp = cs_field_get_key_int(CS_F_(alp_bl), var_key_id) -1;

    cs_real_t  *_rcodcl_alp = rcodcl + ivar_alp*n_b_faces;

    for (int izone = 0; izone < boundaries->n_zones; izone++) {
      if (cs_gui_strcmp(boundaries->nature[izone], "wall")) {

        int zone_nbr = boundaries->bc_num[izone];
        const cs_zone_t *z = cs_boundary_zone_by_id(zone_nbr);

        for (cs_lnum_t i = 0; i < z->n_elts; i++) {
          cs_lnum_t face_id = z->elt_ids[i];
          _rcodcl_alp[face_id] = 0.;
        }

      }
    }

  }

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
 * integer          nozppm   <-- max number of boundary conditions zone
 * integer          nclpch   <-- number of simulated class per coals
 * integer          iqimp    <-- 1 if mass flow rate is applied
 * integer          ientat   <-- 1 for air temperature boundary conditions (coal)
 * integer          ientcp   <-- 1 for coal temperature boundary conditions (coal)
 * integer          inmoxy   <-- coal: number of oxydant for the current inlet
 * integer          ientox   <-- 1 for an air fow inlet (gas combustion)
 * integer          ientfu   <-- 1 for fuel flow inlet (gas combustion)
 * integer          ientgb   <-- 1 for burned gas inlet (gas combustion)
 * integer          ientgf   <-- 1 for unburned gas inlet (gas combustion)
 * integer          iprofm   <-- atmospheric flows: on/off for profile from data
 * integer          iautom   <-- atmospheric flows: auto inlet/outlet flag
 * integer          itypfb   <-- type of boundary for each face
 * integer          izfppp   <-- zone number for each boundary face
 * integer          icodcl   <-- boundary conditions array type
 * double precision cgdfbo   <-- boundary faces center of gravity
 * double precision qimp     <-- inlet flow rate
 * double precision qimpat   <-- inlet air flow rate (coal)
 * double precision qimpcp   <-- inlet coal flow rate (coal)
 * double precision timpat   <-- air temperature boundary conditions (coal)
 * double precision timpcp   <-- inlet coal temperature (coal)
 * double precision tkent    <-- inlet temperature (gas combustion)
 * double precision fment    <-- Mean Mixture Fraction at Inlet (gas combustion)
 * double precision distch   <-- ratio for each coal
 * integer          nvar     <-- dimension for rcodcl
 * double precision rcodcl   <-- boundary conditions array value
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclim, UICLIM)(const int  *nozppm,
                               int        *iqimp,
                               int        *icalke,
                               int        *ientat,
                               int        *ientcp,
                               int        *inmoxy,
                               int        *ientox,
                               int        *ientfu,
                               int        *ientgf,
                               int        *ientgb,
                               int        *iprofm,
                               int        *iautom,
                               int        *itypfb,
                               int        *izfppp,
                               int        *icodcl,
                               double     *qimp,
                               double     *qimpat,
                               double     *qimpcp,
                               double     *dh,
                               double     *xintur,
                               double     *timpat,
                               double     *timpcp,
                               double     *tkent,
                               double     *fment,
                               double     *distch,
                               int        *nvar,
                               double     *rcodcl)
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  const int var_key_id = cs_field_key_id("variable_id");

  const int ncharm = CS_COMBUSTION_MAX_COALS;

  bool solid_fuels = false;
  if (   cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] > -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1)
    solid_fuels = true;
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
    _init_zones(n_b_faces, nozppm, izfppp);
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

      /* Update the zone's arrays (iqimp, dh, xintur, icalke, qimp,...)
         because they are re-initialized at each time step
         in PRECLI and PPPRCL routines */

      /* data by zone */
      iqimp[zone_nbr-1]  = boundaries->iqimp[izone];
      dh[zone_nbr-1]     = boundaries->dh[izone];
      xintur[zone_nbr-1] = boundaries->xintur[izone];
      icalke[zone_nbr-1] = boundaries->icalke[izone];

      if (solid_fuels) {
        const cs_combustion_model_t *cm = cs_glob_combustion_model;

        ientat[zone_nbr-1] = boundaries->ientat[izone];
        inmoxy[zone_nbr-1] = boundaries->inmoxy[izone];
        ientcp[zone_nbr-1] = boundaries->ientcp[izone];
        timpat[zone_nbr-1] = boundaries->timpat[izone];

        for (int icharb = 0; icharb < cm->coal.n_coals; icharb++) {
          int ich = icharb *(*nozppm)+zone_nbr-1;
          qimpcp[ich] = boundaries->qimpcp[izone][icharb];
          timpcp[ich] = boundaries->timpcp[izone][icharb];

          for (int iclass = 0;
               iclass < cm->coal.n_classes_per_coal[icharb];
               iclass++) {
            int icl = iclass*(*nozppm) * ncharm + ich;
            distch[icl] = boundaries->distch[izone][icharb][iclass];
          }
        }
      }
      else if (gas_combustion) {
        ientfu[zone_nbr-1] = boundaries->ientfu[izone];
        ientox[zone_nbr-1] = boundaries->ientox[izone];
        ientgb[zone_nbr-1] = boundaries->ientgb[izone];
        ientgf[zone_nbr-1] = boundaries->ientgf[izone];
        tkent[zone_nbr-1]  = boundaries->tkent[izone];
        fment[zone_nbr-1]  = boundaries->fment[izone];
      }
      else if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {

        if (  boundaries->itype[izone] == CS_ESICF
            ||boundaries->itype[izone] == CS_EPHCF) {
          const cs_field_t  *fp = cs_field_by_name_try("pressure");
          int ivarp = cs_field_get_key_int(fp, var_key_id) -1;

          for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
            cs_lnum_t face_id = bz->elt_ids[elt_id];
            rcodcl[ivarp * n_b_faces + face_id] = boundaries->prein[izone];
          }
        }

        if (boundaries->itype[izone] == CS_ESICF) {
          cs_field_t *b_rho = cs_field_by_name_try("boundary_density");
          const cs_field_t  *ft = cs_field_by_name_try("temperature");
          int ivart = cs_field_get_key_int(ft, var_key_id) -1;

          for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
            cs_lnum_t face_id = bz->elt_ids[elt_id];
            rcodcl[ivart * n_b_faces + face_id] = boundaries->tempin[izone];
            b_rho->val[face_id] = boundaries->rhoin[izone];
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
        izfppp[face_id] = zone_nbr;
        itypfb[face_id] = inlet_type;
      }

      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
        iprofm[zone_nbr-1] = boundaries->meteo[izone].read_data;
        if (boundaries->meteo[izone].automatic) {
          for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
            cs_lnum_t face_id = bz->elt_ids[elt_id];
            iautom[face_id] = 1;
          }
        }
      }

      /* turbulent inlet, with formula */
      if (boundaries->icalke[zone_nbr-1] == 0) {

        tn_bc = _get_zone_bc_node(tn_bc, izone);

        cs_tree_node_t *tn_t = cs_tree_node_get_child(tn_bc, "turbulence");
        const char *formula = cs_tree_node_get_child_value_str(tn_t, "formula");

        if (formula != NULL) {

          const char *model = cs_gui_get_thermophysical_model("turbulence");
          if (model == NULL)
            return;

          if (   cs_gui_strcmp(model, "k-epsilon")
              || cs_gui_strcmp(model, "k-epsilon-PL")) {

            cs_real_t *new_vals = cs_meg_boundary_function(bz,
                                                           "turbulence_ke",
                                                           "formula");

            cs_field_t *c_k   = cs_field_by_name("k");
            cs_field_t *c_eps = cs_field_by_name("epsilon");
            int ivark = cs_field_get_key_int(c_k, var_key_id) -1;
            int ivare = cs_field_get_key_int(c_eps, var_key_id) -1;

            for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
              cs_lnum_t face_id = bz->elt_ids[elt_id];
              rcodcl[ivark * n_b_faces + face_id]
                = new_vals[0 * bz->n_elts + elt_id];
              rcodcl[ivare * n_b_faces + face_id]
                = new_vals[1 * bz->n_elts + elt_id];
            }
            BFT_FREE(new_vals);
          }
          else if (  cs_gui_strcmp(model, "Rij-epsilon")
                   ||cs_gui_strcmp(model, "Rij-SSG")) {

            cs_real_t *new_vals = cs_meg_boundary_function(bz,
                                                           "turbulence_rije",
                                                           "formula");
            cs_field_t *cfld_rij = cs_field_by_name("rij");
            cs_field_t *c_eps = cs_field_by_name("epsilon");
            int ivarrij = cs_field_get_key_int(cfld_rij, var_key_id) - 1;
            int ivare   = cs_field_get_key_int(c_eps, var_key_id) -1;

            for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
              cs_lnum_t face_id = bz->elt_ids[elt_id];

              /* Values are stored for rij components then epsilon */
              for (int ii = 0; ii < 6; ii++)
                rcodcl[(ivarrij + ii) * n_b_faces + face_id]
                  = new_vals[bz->n_elts * ii + elt_id];

              rcodcl[ivare * n_b_faces + face_id]
                = new_vals[bz->n_elts * 6 + elt_id];
            }
            BFT_FREE(new_vals);
          }
          else if (cs_gui_strcmp(model, "Rij-EBRSM")) {

            cs_real_t *new_vals = cs_meg_boundary_function(bz,
                                                           "turbulence_rij_ebrsm",
                                                           "formula");

            cs_field_t *cfld_rij = cs_field_by_name("rij");
            cs_field_t *c_eps = cs_field_by_name("epsilon");
            cs_field_t *c_a   = cs_field_by_name("alpha");

            int ivarrij = cs_field_get_key_int(cfld_rij, var_key_id) - 1;
            int ivare   = cs_field_get_key_int(c_eps, var_key_id) -1;
            int ivara   = cs_field_get_key_int(c_a, var_key_id) -1;

            for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
              cs_lnum_t face_id = bz->elt_ids[elt_id];

              /* Values are stored for rij components then epsilon and alpha*/
              for (int ii = 0; ii < 6; ii++)
                rcodcl[(ivarrij + ii) * n_b_faces + face_id]
                  = new_vals[bz->n_elts * ii + elt_id];

              rcodcl[ivare * n_b_faces + face_id]
                = new_vals[bz->n_elts * 6 + elt_id];

              rcodcl[ivara   * n_b_faces + face_id]
                = new_vals[bz->n_elts * 7 + elt_id];
            }
            BFT_FREE(new_vals);
          }
          else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
            cs_real_t *new_vals = cs_meg_boundary_function(bz,
                                                           "turbulence_v2f",
                                                           "formula");

            cs_field_t *c_k   = cs_field_by_name("k");
            cs_field_t *c_eps = cs_field_by_name("epsilon");
            cs_field_t *c_phi = cs_field_by_name("phi");
            cs_field_t *c_a   = cs_field_by_name("alpha");
            int ivark = cs_field_get_key_int(c_k,   var_key_id) -1;
            int ivare = cs_field_get_key_int(c_eps, var_key_id) -1;
            int ivarp = cs_field_get_key_int(c_phi, var_key_id) -1;
            int ivara = cs_field_get_key_int(c_a,   var_key_id) -1;

            for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
              cs_lnum_t face_id = bz->elt_ids[elt_id];
              rcodcl[ivark * n_b_faces + face_id]
                = new_vals[0 * bz->n_elts + elt_id];
              rcodcl[ivare * n_b_faces + face_id]
                = new_vals[1 * bz->n_elts + elt_id];
              rcodcl[ivarp * n_b_faces + face_id]
                = new_vals[2 * bz->n_elts + elt_id];
              rcodcl[ivara * n_b_faces + face_id]
                = new_vals[3 * bz->n_elts + elt_id];
            }
            BFT_FREE(new_vals);
          }
          else if (cs_gui_strcmp(model, "k-omega-SST")) {
            cs_real_t *new_vals = cs_meg_boundary_function(bz,
                                                           "turbulence_kw",
                                                           "formula");

            cs_field_t *c_k = cs_field_by_name("k");
            cs_field_t *c_o = cs_field_by_name("omega");
            int ivark = cs_field_get_key_int(c_k,   var_key_id) -1;
            int ivaro = cs_field_get_key_int(c_o,   var_key_id) -1;

            for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
              cs_lnum_t face_id = bz->elt_ids[elt_id];
              rcodcl[ivark * n_b_faces + face_id]
                = new_vals[0 * bz->n_elts + elt_id];
              rcodcl[ivaro * n_b_faces + face_id]
                = new_vals[1 * bz->n_elts + elt_id];
            }
            BFT_FREE(new_vals);
          }
          else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
            cs_real_t *new_vals = cs_meg_boundary_function(bz,
                                                           "turbulence_spalart",
                                                           "formula");

            cs_field_t *c_nu = cs_field_by_name("nu_tilda");
            int ivarnu = cs_field_get_key_int(c_nu, var_key_id) -1;

            for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
              cs_lnum_t face_id = bz->elt_ids[elt_id];
              rcodcl[ivarnu * n_b_faces + face_id] = new_vals[elt_id];
            }
            BFT_FREE(new_vals);
          }
          else
            bft_error(__FILE__, __LINE__, 0,
                      _("Invalid turbulence model: %s.\n"), model);
        }
      }

#if _XML_DEBUG_
      if (solid_fuels) {
        const cs_combustion_model_t *cm = cs_glob_combustion_model;

        bft_printf("-----ientat=%i, ientcp=%i, timpat=%12.5e \n",
            ientat[zone_nbr-1], ientcp[zone_nbr-1], timpat[zone_nbr-1]);

        for (int icharb = 0; icharb < cm->coal.n_coals; icharb++) {
          bft_printf("-----coal=%i, qimpcp=%12.5e, timpcp=%12.5e \n",
                     icharb+1, qimpcp[icharb *(*nozppm)+zone_nbr-1],
                     timpcp[icharb *(*nozppm)+zone_nbr-1]);

          for (int iclass = 0; iclass < cm->coal.n_coals; iclass++)
            bft_printf("-----coal=%i, class=%i, distch=%f \n",
                       icharb+1, iclass+1,
                       distch[  iclass * (*nozppm) * ncharm
                              + icharb * (*nozppm) +zone_nbr-1]);
        }
      }
      else if (gas_combustion) {
        bft_printf("-----iqimp=%i \n",
                   iqimp[zone_nbr-1]);
        bft_printf("-----ientox=%i, ientfu=%i, ientgf=%i, ientgb=%i \n",
                   ientox[zone_nbr-1], ientfu[zone_nbr-1],
                   ientgf[zone_nbr-1], ientgb[zone_nbr-1]);
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
                   iqimp[zone_nbr-1], qimp[zone_nbr-1]);
      }

      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
        bft_printf("-----iprofm=%i, automatic=%i \n",
                   iprofm[zone_nbr-1], boundaries->meteo[izone].automatic);
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
        izfppp[face_id] = zone_nbr;
        itypfb[face_id] = iwall;
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "outlet")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        izfppp[face_id] = zone_nbr;
        if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1)
          itypfb[face_id] = boundaries->itype[izone];
        else
          itypfb[face_id] = CS_OUTLET;
      }

      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
        iprofm[zone_nbr-1] = boundaries->meteo[izone].read_data;
        if (boundaries->meteo[izone].automatic) {
          for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
            cs_lnum_t face_id = bz->elt_ids[elt_id];
            iautom[face_id] = 1;
          }
        }
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "imposed_p_outlet")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_OUTLET;
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "symmetry")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_SYMMETRY;
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "free_inlet_outlet")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_FREE_INLET;
      }

      if (boundaries->head_loss_e[izone]) {
        cs_real_t *new_vals =
          cs_meg_boundary_function(bz, "head_loss", "formula");

        const cs_field_t  *fp = cs_field_by_name_try("pressure");
        int ivarp = cs_field_get_key_int(fp, var_key_id) -1;

        for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
          cs_lnum_t face_id = bz->elt_ids[elt_id];
          rcodcl[1 * n_b_faces * (*nvar) + ivarp * n_b_faces + face_id]
            = new_vals[elt_id];
        }
        BFT_FREE(new_vals);
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "free_surface")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_FREE_SURFACE;
      }
    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "groundwater")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        izfppp[face_id] = zone_nbr;
        itypfb[face_id] = CS_INDEF;
      }

      /* set velocity to 0 */
      const cs_field_t  *fp2 = cs_field_by_name_try("velocity");
      if (fp2 != NULL) {
        int ivar2 = cs_field_get_key_int(fp2, var_key_id) -1;

        for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
          cs_lnum_t face_id = bz->elt_ids[elt_id];
          for (cs_lnum_t i = 0; i < 3; i++) {
            icodcl[(ivar2 + i) * n_b_faces + face_id] = 3;
            rcodcl[(ivar2 + i) * n_b_faces + face_id] = 0.;
          }
        }
      }

    }

    else if (cs_gui_strcmp(boundaries->nature[izone], "undefined")) {
      for (cs_lnum_t elt_id = 0; elt_id < bz->n_elts; elt_id++) {
        cs_lnum_t face_id = bz->elt_ids[elt_id];
        izfppp[face_id] = zone_nbr;
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
        cs_lnum_t ivar = cs_field_get_key_int(f, var_key_id) -1;
        if (f->type & CS_FIELD_VARIABLE) {
          bft_printf("------%s: icodcl=%i, "
                     "rcodcl(1)=%12.5e, rcodcl(2)=%12.5e, rcodcl(3)=%12.5e\n",
                     f->name,
                     icodcl[ivar *n_b_faces +face_id ],
                     rcodcl[0 * n_b_faces * (*nvar) +ivar * n_b_faces +face_id],
                     rcodcl[1 * n_b_faces * (*nvar) +ivar * n_b_faces +face_id],
                     rcodcl[2 * n_b_faces * (*nvar) +ivar * n_b_faces +face_id]);
        }
      }
    }
#endif

  } /*  for (izone=0; izone < boundaries->n_zones; izone++) */

  /* Define boundary conditions based on cs_equation_param_t structures */

  cs_boundary_conditions_compute(*nvar,
                                 itypfb,
                                 icodcl,
                                 rcodcl);

  /* Final adjustments */

  /* Automatic turbulence values
     TODO: spcecify which models do not need this instead of those who do.
     In most cases, even if a model changes this, doing it in all cases
     should be safe. */

  bool std_turbulence_bcs = false;
  if (   cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] == 0
      || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 0
      || cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0
      || cs_glob_physical_model_flag[CS_GAS_MIX] >= 0)
    std_turbulence_bcs = true;

  if (std_turbulence_bcs)
    _standard_turbulence_bcs(rcodcl);

  /* Mass flow (TODO: merge qimp and qimpat) */

  for (int izone = 0; izone < boundaries->n_zones; izone++) {
    qimp[izone] = boundaries->qimp[izone];
    if (solid_fuels)
      qimpat[izone] = boundaries->qimp[izone];
  }
}

/*----------------------------------------------------------------------------
 * Boundary conditions input verification
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICLVE
 * *****************
 *
 * integer          nozppm  <-- max number of boundary conditions zone
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclve, UICLVE)(const int  *nozppm)
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
        && (cs_glob_ale != 0)) {
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

    if (nozppm && zone_nbr > *nozppm)
      bft_error(__FILE__, __LINE__, 0,
                _("zone's label number %i is greater than %i,"
                  " the maximum allowed \n"), zone_nbr, *nozppm);

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
  int izone;
  int n_zones;
  int icharb;

  /* clean memory for global private structure boundaries */

  if (boundaries != NULL) {

    n_zones = boundaries->n_zones;

    for (int f_id = 0; f_id < boundaries->n_fields; f_id++) {
      const cs_field_t  *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        if (boundaries->type_code != NULL)
          BFT_FREE(boundaries->type_code[f->id]);
      }
    }

    if (   cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] > -1
        || cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1) {
      const int n_coals = cs_glob_combustion_model->coal.n_coals;
      for (izone = 0; izone < n_zones; izone++) {
        BFT_FREE(boundaries->qimpcp[izone]);
        BFT_FREE(boundaries->timpcp[izone]);
        for (icharb = 0; icharb < n_coals; icharb++)
          BFT_FREE(boundaries->distch[izone][icharb]);
        BFT_FREE(boundaries->distch[izone]);
      }
      BFT_FREE(boundaries->ientat);
      BFT_FREE(boundaries->ientcp);
      BFT_FREE(boundaries->inmoxy);
      BFT_FREE(boundaries->timpat);
      BFT_FREE(boundaries->qimpcp);
      BFT_FREE(boundaries->timpcp);
      BFT_FREE(boundaries->distch);
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

    BFT_FREE(boundaries->iqimp);
    BFT_FREE(boundaries->icalke);
    BFT_FREE(boundaries->qimp);
    BFT_FREE(boundaries->dh);
    BFT_FREE(boundaries->xintur);
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

END_C_DECLS
