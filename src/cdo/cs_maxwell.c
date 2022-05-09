/*============================================================================
 * Handle Maxwell module with CDO schemes
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

#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

#include "cs_hodge.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_maxwell.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_maxwell.h

  \brief Structure and functions handling the Maxwell module dedicated to
         the resolution of electro-magnetic equations
*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_MAXWELL_DBG     0

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Set of parameters related to the Maxwell module */

struct _maxwell_t {

  cs_flag_t   model;      /* Modelling for the Maxwell module */
  cs_flag_t   options;    /* Flag dedicated to general options to handle
                           * the Maxwell module */
  cs_flag_t   post_flag;  /* Flag dedicated to the post-processing
                           * of the Maxwell module */

  /* Properties associated to this module */

  cs_real_t       e_perm_ref;      /* Reference value of the electric
                                    * permeability */
  cs_property_t  *e_permeability;  /* Electric permeability oftenly denoted
                                    * by epsilon */

  cs_real_t       m_perm_ref;      /* Reference value of the magnetic
                                    * permittivity */
  cs_property_t  *m_permittivity;  /* Magnetic permittivity oftenly denoted
                                    * by mu */

  cs_property_t  *conductivity;    /* Conductivity in Ohm's law oftenly denoted
                                    * by sigma */

  /* Fields associated to this module */

  cs_field_t   *scal_pot;        /* Scalar potential at vertices called
                                    "electric_potential" */

  cs_field_t   *vect_pot;        /* Vector potential */

  cs_field_t   *e_field;         /* E: Electric field */
  cs_real_t    *e_field_array;   /* E: Electric field along edges */

  cs_field_t   *d_field;         /* D: Electric induction field
                                    electric flux density (As/m^2) */
  cs_real_t    *d_field_array;   /* D: Electric induction field across dual
                                    faces */

  cs_field_t   *h_field;         /* H = Magnetic field (faces) */
  cs_real_t    *h_field_array;   /* H along dual edges */

  cs_field_t   *b_field;         /* B = Magnetic induction field */
  cs_real_t    *b_field_array;   /* B across faces */

  cs_field_t   *j_field;         /* J = density flux field */
  cs_real_t    *j_field_array;   /* J across dual faces */

  cs_field_t   *joule_effect;    /* Joule effect power (W.m^-3) */
};

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_maxwell[] =
  " Stop execution.\n"
  " The structure related to the Maxwell module is empty.\n"
  " Please check your settings.\n";

static cs_maxwell_t  *cs_maxwell_structure = NULL;

/* Vacuum magnetic permeability constant (H/m) */

static const cs_real_t  cs_maxwell_vacuum_m_permeability = 1.25663706143591e-6;

/* Vacuum permittivity constant (F/m) */

static const cs_real_t  cs_maxwell_vacuum_e_permittivity = 8.85418782e-12;

/*============================================================================
 * Private static inline function prototypes
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cellwise constant vector-valued approximation of fields
 *
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      c2e          pointer to a cs_adjacency_t structure
 * \param[in]      ep_values    array of edge values
 * \param[in]      df_values    array of dual face values
 * \param[in, out] c_ep_values  array vector-valued cell arrays
 * \param[in, out] c_fd_values  array vector-valued cell arrays
 */
/*----------------------------------------------------------------------------*/

static void
_build_edge_based_vector_fields(const cs_cdo_quantities_t   *quant,
                                const cs_adjacency_t        *c2e,
                                const cs_real_t             *ep_values,
                                const cs_real_t             *fd_values,
                                cs_real_t                   *c_ep_values,
                                cs_real_t                   *c_fd_values)
{
  assert(ep_values != NULL && fd_values != NULL);
  assert(c_ep_values != NULL && c_fd_values != NULL);

  memset(c_ep_values, 0, 3*quant->n_cells*sizeof(cs_real_t));
  memset(c_fd_values, 0, 3*quant->n_cells*sizeof(cs_real_t));

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    cs_real_t  *c_ep = c_ep_values + 3*c_id;
    cs_real_t  *c_fd = c_fd_values + 3*c_id;

    for (cs_lnum_t j = c2e->idx[c_id]; j <c2e->idx[c_id+1]; j++) {

      const cs_lnum_t  e_id = c2e->ids[j];
      const cs_real_t  e_val = ep_values[e_id];
      const cs_real_t  *e_vect = quant->edge_vector + 3*e_id;
      const cs_real_t  *dface = quant->dface_normal + 3*j;
      for (int k = 0; k < 3; k++) {
        c_fd[k] += fd_values[j] * e_vect[k];
        c_ep[k] += e_val * dface[k];
      }

    }

    const double  invvol = 1/quant->cell_vol[c_id];
    for (int k = 0; k < 3; k++) {
      c_fd[k] *= invvol;
      c_ep[k] *= invvol;
    }

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cellwise constant vector-valued approximation of fields
 *         for fields associated to (1) primal faces and seen as a flux
 *         (fp_values) and to (2) dual edges and seen as a circulation
 *         (ed_values)
 *
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      c2f          pointer to a cs_adjacency_t structure
 * \param[in]      fp_values    array of (primal) face values
 * \param[in]      ed_values    array of dual edge values
 * \param[in, out] c_fp_values  array vector-valued cell arrays
 * \param[in, out] c_ed_values  array vector-valued cell arrays
 */
/*----------------------------------------------------------------------------*/

static void
_build_face_based_vector_fields(const cs_cdo_quantities_t   *quant,
                                const cs_adjacency_t        *c2f,
                                const cs_real_t             *fp_values,
                                const cs_real_t             *ed_values,
                                cs_real_t                   *c_fp_values,
                                cs_real_t                   *c_ed_values)
{
  assert(fp_values != NULL && ed_values != NULL);
  assert(c_fp_values != NULL && c_ed_values != NULL);

  memset(c_fp_values, 0, 3*quant->n_cells*sizeof(cs_real_t));
  memset(c_ed_values, 0, 3*quant->n_cells*sizeof(cs_real_t));

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    cs_real_t  *c_fp = c_fp_values + 3*c_id;
    cs_real_t  *c_ed = c_ed_values + 3*c_id;

    for (cs_lnum_t j = c2f->idx[c_id]; j <c2f->idx[c_id+1]; j++) {

      const cs_lnum_t  f_id = c2f->ids[j];
      const cs_nvec3_t  pfq = cs_quant_set_face_nvec(f_id, quant);
      const cs_real_t  ed_coef = ed_values[j] * pfq.meas;
      const cs_real_t  *ed_vect = quant->dedge_vector + 3*j;
      const cs_real_t  f_val = fp_values[f_id];

      for (int k = 0; k < 3; k++) {
        c_ed[k] += ed_coef * pfq.unitv[k];
        c_fp[k] += f_val * ed_vect[k];
      }

    }

    const double  invvol = 1/quant->cell_vol[c_id];
    for (int k = 0; k < 3; k++) {
      c_ed[k] *= invvol;
      c_fp[k] *= invvol;
    }

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a structure dedicated to manage the Maxwell module
 *
 * \return a pointer to a new allocated cs_maxwell_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_maxwell_t *
_maxwell_create(void)
{
  cs_maxwell_t  *mxl = NULL;

  BFT_MALLOC(mxl, 1, cs_maxwell_t);

  /* Default initialization */

  mxl->model = 0;
  mxl->options = 0;
  mxl->post_flag = 0;

  /* Properties */

  mxl->e_perm_ref = cs_maxwell_vacuum_e_permittivity;
  mxl->e_permeability = NULL;

  mxl->m_perm_ref = cs_maxwell_vacuum_m_permeability;
  mxl->m_permittivity = NULL;

  mxl->conductivity = NULL;

  /* Fields and related arrays */

  mxl->scal_pot = NULL;

  mxl->vect_pot = NULL;

  mxl->e_field = NULL;
  mxl->e_field_array = NULL;    /* different location that e_field */

  mxl->d_field = NULL;
  mxl->d_field_array = NULL;    /* different location that d_field */

  mxl->h_field = NULL;
  mxl->h_field_array = NULL;    /* different location that h_field */

  mxl->b_field = NULL;
  mxl->b_field_array = NULL;    /* different location that b_field */

  mxl->j_field = NULL;
  mxl->j_field_array = NULL;    /* different location that j_field */

  /* Additional quantities (source terms for instance) */

  mxl->joule_effect = NULL;

  return mxl;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if the computation of Maxwell equations is activated
 */
/*----------------------------------------------------------------------------*/

bool
cs_maxwell_is_activated(void)
{
  if (cs_maxwell_structure == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the future computation of the Maxwell equations
 *
 * \param[in]      model         type of modelling
 * \param[in]      options       flag to handle optional parameters
 *
 * \return a pointer to a new allocated Maxwell structure
 */
/*----------------------------------------------------------------------------*/

cs_maxwell_t *
cs_maxwell_activate(cs_flag_t     model,
                    cs_flag_t     options)
{
  if (model < 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid modelling. Model = %d\n", __func__, model);

  /* Allocate an empty structure */

  cs_maxwell_t  *mxl = _maxwell_create();

  /* Set members of the structure according to the given settings */

  mxl->model = model;
  mxl->options = options;

  if (model & CS_MAXWELL_MODEL_ELECTROSTATIC) {

    cs_equation_t  *e_static = cs_equation_add(CS_MAXWELL_ESTATIC_EQNAME,
                                               "electric_potential",
                                               CS_EQUATION_TYPE_MAXWELL,
                                               1,
                                               CS_PARAM_BC_HMG_NEUMANN);
    cs_equation_param_t  *eqp = cs_equation_get_param(e_static);

    mxl->e_permeability = cs_property_add("electric_permeability",
                                          CS_PROPERTY_ISO);

    /* By default, set the reference permeability */

    cs_property_def_iso_by_value(mxl->e_permeability, NULL, mxl->e_perm_ref);

    cs_equation_add_diffusion(eqp, mxl->e_permeability);

    /* Should be symmetric */

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vb");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "bubble");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_COEF, "frac23");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_FAMILY, "cs");
    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "cg");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL_EPS, "1e-6");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL_RESNORM_TYPE, "filtered");

  }

  if (model & CS_MAXWELL_MODEL_MAGNETOSTATIC) {

    cs_equation_t  *m_static = cs_equation_add(CS_MAXWELL_MSTATIC_EQNAME,
                                               "magnetic_potential",
                                               CS_EQUATION_TYPE_MAXWELL,
                                               3,
                                               CS_PARAM_BC_HMG_DIRICHLET);

    cs_equation_param_t  *eqp = cs_equation_get_param(m_static);

    mxl->m_permittivity = cs_property_add("magnetic_permittivity",
                                          CS_PROPERTY_ISO);

    /* By default, set the reference permeability */

    cs_property_def_iso_by_value(mxl->m_permittivity, NULL, mxl->m_perm_ref);

    cs_equation_add_curlcurl(eqp, mxl->m_permittivity,
                             1); /* Inverse of the property is requested */

    /* Should be symmetric */

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_eb");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "cost");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_COEF, "dga");
    cs_equation_param_set(eqp, CS_EQKEY_SOLVER_FAMILY, "cs");
    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "cg");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL_EPS, "1e-6");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL_RESNORM_TYPE, "filtered");

  }

  cs_maxwell_structure = mxl;

  return mxl;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the Maxwell module
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_maxwell_t *
cs_maxwell_destroy_all(void)
{
  if (cs_maxwell_structure == NULL)
    return NULL;

  cs_maxwell_t  *mxl = cs_maxwell_structure;

  /* The lifecycle of properties and fields is not managed by the current
     structure.
     Free only arrays which are owned by this structure */

  BFT_FREE(mxl->e_field_array);
  BFT_FREE(mxl->d_field_array);
  BFT_FREE(mxl->h_field_array);
  BFT_FREE(mxl->b_field_array);
  BFT_FREE(mxl->j_field_array);

  BFT_FREE(mxl);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup equations/properties related to the Maxwell module
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_init_setup(void)
{
  cs_maxwell_t  *mxl = cs_maxwell_structure;

  if (mxl == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_maxwell));

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_CDO;
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");

  if (mxl->model & CS_MAXWELL_MODEL_ELECTROSTATIC) {

    mxl->e_field = cs_field_create(CS_MAXWELL_EFIELD_NAME,
                                   field_mask,
                                   CS_MESH_LOCATION_CELLS,
                                   3,
                                   true);

    cs_field_set_key_int(mxl->e_field, log_key, 1);
    cs_field_set_key_int(mxl->e_field, post_key, 1);

    mxl->d_field = cs_field_create(CS_MAXWELL_DFIELD_NAME,
                                   field_mask,
                                   CS_MESH_LOCATION_CELLS,
                                   3,
                                   true);

    cs_field_set_key_int(mxl->d_field, log_key, 1);
    cs_field_set_key_int(mxl->d_field, post_key, 1);

    /* Add the variable field */

    cs_equation_t  *eq = cs_equation_by_name(CS_MAXWELL_ESTATIC_EQNAME);

    cs_equation_predefined_create_field(-1, eq); /* automatic */

  }

  if (mxl->model & CS_MAXWELL_MODEL_MAGNETOSTATIC) {

    mxl->b_field = cs_field_create(CS_MAXWELL_BFIELD_NAME,
                                   field_mask,
                                   CS_MESH_LOCATION_CELLS,
                                   3,
                                   true);

    cs_field_set_key_int(mxl->b_field, log_key, 1);
    cs_field_set_key_int(mxl->b_field, post_key, 1);

    mxl->h_field = cs_field_create(CS_MAXWELL_MFIELD_NAME,
                                   field_mask,
                                   CS_MESH_LOCATION_CELLS,
                                   3,
                                   true);

    cs_field_set_key_int(mxl->h_field, log_key, 1);
    cs_field_set_key_int(mxl->h_field, post_key, 1);

    /* Add the variable field */

    cs_equation_t  *eq = cs_equation_by_name(CS_MAXWELL_MSTATIC_EQNAME);

    cs_equation_predefined_create_field(-1, eq); /* automatic */

  }

  /* Optional settings */

  if (mxl->options & CS_MAXWELL_JOULE_EFFECT) {

    mxl->joule_effect = cs_field_create(CS_MAXWELL_JEFFECT_NAME,
                                        field_mask,
                                        CS_MESH_LOCATION_CELLS,
                                        1,
                                        true);
    cs_field_set_key_int(mxl->joule_effect, log_key, 1);
    cs_field_set_key_int(mxl->joule_effect, post_key, 1);

  }

  /* Add default post-processing related to the Maxwell module */

  cs_post_add_time_mesh_dep_output(cs_maxwell_extra_post, mxl);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup stage for equations related to the Maxwell module
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_finalize_setup(const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant)
{
  cs_maxwell_t  *mxl = cs_maxwell_structure;

  if (mxl == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_maxwell));

  if (mxl->model & CS_MAXWELL_MODEL_ELECTROSTATIC) {

    cs_equation_t  *es_eq = cs_equation_by_name(CS_MAXWELL_ESTATIC_EQNAME);
    assert(cs_equation_get_space_scheme(es_eq) == CS_SPACE_SCHEME_CDOVB);

    mxl->scal_pot = cs_equation_get_field(es_eq);

    /* Electric field array along edges */

    BFT_MALLOC(mxl->e_field_array, quant->n_edges, cs_real_t);
    memset(mxl->e_field_array, 0, quant->n_edges*sizeof(cs_real_t));

    /* Electric induction (flux density) across dual faces */

    const cs_adjacency_t  *c2e = connect->c2e;
    const cs_lnum_t  array_size = c2e->idx[quant->n_cells];
    BFT_MALLOC(mxl->d_field_array, array_size, cs_real_t);
    memset(mxl->d_field_array, 0, array_size*sizeof(cs_real_t));

  }

  if (mxl->model & CS_MAXWELL_MODEL_MAGNETOSTATIC) {

    cs_equation_t  *ms_eq = cs_equation_by_name(CS_MAXWELL_MSTATIC_EQNAME);
    assert(cs_equation_get_space_scheme(ms_eq) == CS_SPACE_SCHEME_CDOEB);

    mxl->vect_pot = cs_equation_get_field(ms_eq);

    /* Magnetic field array along dual edges */

    const cs_adjacency_t  *c2f = connect->c2f;
    const cs_lnum_t  array_size = c2f->idx[quant->n_cells];
    BFT_MALLOC(mxl->h_field_array, array_size, cs_real_t);
    memset(mxl->h_field_array, 0, array_size*sizeof(cs_real_t));

    /* Magnetic induction (flux density) across primal faces */

    BFT_MALLOC(mxl->b_field_array, quant->n_faces, cs_real_t);
    memset(mxl->b_field_array, 0, quant->n_faces*sizeof(cs_real_t));

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log a summary of the Maxwell module
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_log_setup(void)
{
  cs_maxwell_t  *mxl = cs_maxwell_structure;

  if (mxl == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the Maxwell module\n");
  cs_log_printf(CS_LOG_SETUP, "%s\n", cs_sep_h1);

  cs_log_printf(CS_LOG_SETUP, "  * Maxwell | Model:");
  if (mxl->model & CS_MAXWELL_MODEL_ELECTROSTATIC)
    cs_log_printf(CS_LOG_SETUP, "  Electro-static");
  if (mxl->model & CS_MAXWELL_MODEL_MAGNETOSTATIC)
    cs_log_printf(CS_LOG_SETUP, "+  Magneto-static");
  cs_log_printf(CS_LOG_SETUP, "\n");

  if (mxl->options & CS_MAXWELL_JOULE_EFFECT)
    cs_log_printf(CS_LOG_SETUP, "  * Maxwell | Joule effect\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve if needed the steady-state equations related to the Maxwell
 *         module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_compute_steady_state(const cs_mesh_t              *mesh,
                                const cs_time_step_t         *time_step,
                                const cs_cdo_connect_t       *connect,
                                const cs_cdo_quantities_t    *quant)
{
  cs_maxwell_t  *mxl = cs_maxwell_structure;

  if (mxl == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_maxwell));

  if (mxl->model & CS_MAXWELL_MODEL_ELECTROSTATIC) {

    cs_equation_t  *es_eq = cs_equation_by_name(CS_MAXWELL_ESTATIC_EQNAME);

    assert(es_eq != NULL);
    assert(cs_equation_uses_new_mechanism(es_eq));
    cs_equation_solve_steady_state(mesh, es_eq);

  }

  if (mxl->model & CS_MAXWELL_MODEL_MAGNETOSTATIC) {

    cs_equation_t  *ms_eq = cs_equation_by_name(CS_MAXWELL_MSTATIC_EQNAME);

    assert(ms_eq != NULL);
    assert(cs_equation_uses_new_mechanism(ms_eq));
    cs_equation_solve_steady_state(mesh, ms_eq);

  }

  /* Update fields and properties which are related to solved variables */

  cs_maxwell_update(mesh, connect, quant, time_step,
                    true); /* operate current to previous ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve equations related to the Maxwell module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_compute(const cs_mesh_t              *mesh,
                   const cs_time_step_t         *time_step,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *quant)
{
  cs_maxwell_t  *mxl = cs_maxwell_structure;

  if (mxl == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_maxwell));

  /* Add equations to be solved at each time step */

  /* Update fields and properties which are related to solved variables */

  cs_maxwell_update(mesh, connect, quant, time_step,
                    true); /* operate current to previous ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update/initialize the Maxwell module according to the settings
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_update(const cs_mesh_t             *mesh,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *quant,
                  const cs_time_step_t        *ts,
                  bool                         cur2prev)
{
  CS_UNUSED(mesh);

  cs_maxwell_t  *mxl = cs_maxwell_structure;

  if (mxl == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_maxwell));

  if (mxl->model & CS_MAXWELL_MODEL_ELECTROSTATIC) {

    cs_equation_t  *es_eq = cs_equation_by_name(CS_MAXWELL_ESTATIC_EQNAME);

    /* Retrieve the scalar electric potential */

    const cs_real_t  *pot = cs_equation_get_vertex_values(es_eq, false);

    /* Compute the electric field: E = -grad(scal_pot) */

    const cs_adjacency_t  *e2v = connect->e2v;
    for (cs_lnum_t e = 0; e < quant->n_edges; e++) {

      const cs_lnum_t  *v_ids = e2v->ids + 2*e;

      /* E = -grad(scal_pot) */
      mxl->e_field_array[e] = e2v->sgn[2*e]*(pot[v_ids[0]] - pot[v_ids[1]]);

    } /* Loop on edges */

    cs_equation_compute_diffusive_flux(es_eq,
                                       cs_flag_dual_face_byc,
                                       ts->t_cur,
                                       mxl->d_field_array);

    /* Update related vector-valued fields at cell centers */

    if (cur2prev) {
      cs_field_current_to_previous(mxl->e_field);
      cs_field_current_to_previous(mxl->d_field);
    }

    _build_edge_based_vector_fields(quant,
                                    connect->c2e,
                                    mxl->e_field_array, /* in */
                                    mxl->d_field_array, /* in */
                                    mxl->e_field->val,
                                    mxl->d_field->val);

  } /* Electrostatic updates */

  if (mxl->model & CS_MAXWELL_MODEL_MAGNETOSTATIC) {

    cs_equation_t  *ms_eq = cs_equation_by_name(CS_MAXWELL_MSTATIC_EQNAME);
    cs_equation_param_t  *ms_eqp = cs_equation_get_param(ms_eq);

    /* Retrieve the scalar electric potential */

    const cs_real_t  *pot = cs_equation_get_edge_values(ms_eq, false);

    /* Compute the magnetic induction field: B = curl(vect_pot) */

    cs_cdo_connect_discrete_curl(connect, pot, &(mxl->b_field_array));

    /* Compute the magnetic field using Hodge operator */

    cs_hodge_circulation_from_flux(connect, quant, ts->t_cur,
                                   ms_eqp->curlcurl_hodgep,
                                   ms_eqp->curlcurl_property,
                                   mxl->b_field_array,
                                   mxl->h_field_array);

    /* Update related vector-valued fields at cell centers */

    if (cur2prev) {
      cs_field_current_to_previous(mxl->b_field);
      cs_field_current_to_previous(mxl->h_field);
    }

    _build_face_based_vector_fields(quant,
                                    connect->c2f,
                                    mxl->b_field_array, /* in */
                                    mxl->h_field_array, /* in */
                                    mxl->b_field->val,
                                    mxl->h_field->val);

  } /* Magnetostatic updates */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the Maxwell module
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant     pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_extra_op(const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *quant)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_maxwell_t  *mxl = cs_maxwell_structure;

  if (mxl == NULL)
    return;

  /* TODO */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the Maxwell module.
 *         Prototype of this function is fixed since it is a function pointer
 *         defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to an optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_extra_post(void                      *input,
                      int                        mesh_id,
                      int                        cat_id,
                      int                        ent_flag[5],
                      cs_lnum_t                  n_cells,
                      cs_lnum_t                  n_i_faces,
                      cs_lnum_t                  n_b_faces,
                      const cs_lnum_t            cell_ids[],
                      const cs_lnum_t            i_face_ids[],
                      const cs_lnum_t            b_face_ids[],
                      const cs_time_step_t      *time_step)
{
  CS_UNUSED(mesh_id);
  CS_UNUSED(cat_id);
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);
  CS_UNUSED(time_step);

  if (input == NULL)
    return;

  /* TODO */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
