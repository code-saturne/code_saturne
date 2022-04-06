/*============================================================================
 * Manage the enforcement of values at interior degrees of freedom (DoFs) and
 * associated helper functions
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
#include <string.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_bc.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_enforcement.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_enforcement.c

  \brief Structure and functions handling the way to enforce interior degrees
         of freedom.
*/

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_ENFORCEMENT_DBG       0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and define a cs_enforcement_param_t structure
 *
 * \param[in] sel_type   type of elements which have been selected
 * \param[in] type       way to set values for the selected elements
 * \param[in] stride     number of values to enforce by element
 * \param[in] n_elts     number of selected elements locally
 * \param[in] elt_ids    list of element ids
 * \param[in] values     array of values to enforce
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_enforcement_param_create(cs_enforcement_selection_t    sel_type,
                            cs_enforcement_type_t         type,
                            int                           stride,
                            cs_lnum_t                     n_elts,
                            const cs_lnum_t              *elt_ids,
                            const cs_real_t              *values)
{
  cs_enforcement_param_t  *efp = NULL;

  BFT_MALLOC(efp, 1, cs_enforcement_param_t);

  efp->selection_type = sel_type;
  efp->type = type;
  efp->stride = stride;
  efp->n_elts = n_elts;

  if (n_elts > 0 && values == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: No value for the enforcement\n", __func__);

  BFT_MALLOC(efp->elt_ids, n_elts, cs_lnum_t);
  memcpy(efp->elt_ids, elt_ids, n_elts*sizeof(cs_lnum_t));

  switch (type) {

  case CS_ENFORCEMENT_BY_CONSTANT:
    BFT_MALLOC(efp->values, efp->stride, cs_real_t);
    for (int k = 0; k < stride; k++)
      efp->values[k] = values[k];
    break;

  case CS_ENFORCEMENT_BY_DOF_VALUES:
    BFT_MALLOC(efp->values, stride*n_elts, cs_real_t);
    memcpy(efp->values, values, stride*n_elts*sizeof(cs_real_t));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Undefined way to enforce values for interior DoFs\n",
              __func__);

  }

  return efp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset an existing cs_enforcement_param_t structure
 *
 * \param[in, out] efp        pointer to a cs_enforcement_param_t structure
 * \param[in]      sel_type   type of elements which have been selected
 * \param[in]      type       way to set values for the selected elements
 * \param[in]      stride     number of values to enforce by element
 * \param[in]      n_elts     number of selected elements locally
 * \param[in]      elt_ids    list of element ids
 * \param[in]      values     array of values to enforce
 */
/*----------------------------------------------------------------------------*/

void
cs_enforcement_param_reset(cs_enforcement_param_t       *efp,
                           cs_enforcement_selection_t    sel_type,
                           cs_enforcement_type_t         type,
                           int                           stride,
                           cs_lnum_t                     n_elts,
                           const cs_lnum_t              *elt_ids,
                           const cs_real_t              *values)
{
  if (efp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Enforcement param not allocated.\n",
              __func__);
  assert(efp->stride == stride);

  efp->selection_type = sel_type;
  efp->type = type;
  efp->n_elts = n_elts;

  if (n_elts > 0 && values == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: No value for the enforcement\n", __func__);

  BFT_REALLOC(efp->elt_ids, n_elts, cs_lnum_t);
  memcpy(efp->elt_ids, elt_ids, n_elts*sizeof(cs_lnum_t));

  switch (type) {

  case CS_ENFORCEMENT_BY_CONSTANT:
    assert(efp->values != NULL);
    for (int k = 0; k < stride; k++)
      efp->values[k] = values[k];
    break;

  case CS_ENFORCEMENT_BY_DOF_VALUES:
    BFT_REALLOC(efp->values, stride*n_elts, cs_real_t);
    memcpy(efp->values, values, stride*n_elts*sizeof(cs_real_t));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Undefined way to enforce values for interior DoFs\n",
              __func__);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy a cs_enforcement_param_t structure
 *
 * \param[in] ref    reference structure to copy
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_enforcement_param_copy(const cs_enforcement_param_t   *ref)
{
  cs_enforcement_param_t  *dst = NULL;

  if (ref == NULL)
    return dst;

  dst =  cs_enforcement_param_create(ref->selection_type,
                                     ref->type,
                                     ref->stride,
                                     ref->n_elts,
                                     ref->elt_ids,
                                     ref->values);

  return dst;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_enforcement_param_t structure
 *
 * \param[in, out] p_efp    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_enforcement_param_free(cs_enforcement_param_t   **p_efp)
{
  if (p_efp == NULL)
    return;

  cs_enforcement_param_t  *efp = *p_efp;

  if (efp == NULL)
    return;

  BFT_FREE(efp->elt_ids);
  BFT_FREE(efp->values);

  BFT_FREE(efp);
  *p_efp = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log a cs_enforcement_param_t structure
 *
 * \param[in] eqname  name of the related equation
 * \param[in] efp     pointer to a  cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_enforcement_param_log(const char                     *eqname,
                         const cs_enforcement_param_t   *efp)
{
  if (efp == NULL)
    return;

  if (efp->type == CS_ENFORCEMENT_BY_CONSTANT) {

    switch (efp->selection_type) {

    case CS_ENFORCEMENT_SELECTION_CELLS:
      cs_log_printf(CS_LOG_SETUP, "  * %s |   Cell selection | Constant value:",
                    eqname);
      break;

    case CS_ENFORCEMENT_SELECTION_FACES:
      cs_log_printf(CS_LOG_SETUP, "  * %s |   Face selection | Constant value:",
                    eqname);
      break;

    case CS_ENFORCEMENT_SELECTION_EDGES:
      cs_log_printf(CS_LOG_SETUP, "  * %s |   Edge selection | Constant value:",
                    eqname);
      break;

    case CS_ENFORCEMENT_SELECTION_VERTICES:
      cs_log_printf(CS_LOG_SETUP, "  * %s | Vertex selection | Constant value:",
                    eqname);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of selection.",
                __func__);

    } /* Switch on selection type */

    for (int i = 0; i < efp->stride; i++)
      cs_log_printf(CS_LOG_SETUP, " % -6.3e", efp->values[i]);
    cs_log_printf(CS_LOG_SETUP, "\n");

  }
  else if (efp->type == CS_ENFORCEMENT_BY_DOF_VALUES) {

    switch (efp->selection_type) {

    case CS_ENFORCEMENT_SELECTION_CELLS:
      cs_log_printf(CS_LOG_SETUP, "  * %s |   Cell selection | By DoF values\n",
                    eqname);
      break;

    case CS_ENFORCEMENT_SELECTION_FACES:
      cs_log_printf(CS_LOG_SETUP, "  * %s |   Face selection | By DoF values\n",
                    eqname);
      break;

    case CS_ENFORCEMENT_SELECTION_EDGES:
      cs_log_printf(CS_LOG_SETUP, "  * %s |   Edge selection | By DoF values\n",
                    eqname);
      break;

    case CS_ENFORCEMENT_SELECTION_VERTICES:
      cs_log_printf(CS_LOG_SETUP, "  * %s | Vertex selection | By DoF values\n",
                    eqname);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of selection.",
                __func__);

    } /* Switch on selection type */

  }
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of definition",
              __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a cs_enforcement_t structure for vertex-based scheme
 *
 * \param[in] connect    pointer to a cs_cdo_connect_t
 * \param[in] n_params   number of enforcement parameters
 * \param[in] efp_array  array of parameter structures defining the enforcement
 *
 * \return an array with the values to enforce
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_enforcement_define_at_vertices(const cs_cdo_connect_t     *connect,
                                  int                         n_params,
                                  cs_enforcement_param_t    **efp_array)
{
  if (n_params == 0)
    return NULL;

  cs_lnum_t  n_vertices = connect->n_vertices;
  cs_real_t  *values = NULL;
  int  stride = (efp_array[0])->stride;

  BFT_MALLOC(values, stride*n_vertices, cs_real_t);
  for (cs_lnum_t i = 0; i < stride*n_vertices; i++)
    values[i] = FLT_MAX; /* By default, max value */

  /* Define the enforcement value for each vertex related to an enforcement */

  for (int param_id = 0; param_id < n_params; param_id++) {

    cs_enforcement_param_t  *efp = efp_array[param_id];
    assert(stride == efp->stride);

    switch (efp->selection_type) {

    case CS_ENFORCEMENT_SELECTION_VERTICES:

      switch (efp->type) {

      case CS_ENFORCEMENT_BY_CONSTANT:
        for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

          cs_lnum_t  vtx_id = efp->elt_ids[i];
          for (int k = 0; k < stride; k++)
            values[stride*vtx_id + k] = efp->values[k];

        }
        break;

      case CS_ENFORCEMENT_BY_DOF_VALUES:
        for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

          cs_lnum_t  vtx_id = efp->elt_ids[i];
          for (int k = 0; k < stride; k++)
            values[stride*vtx_id + k] = efp->values[stride*i+k];

        }
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of definition.\n",
                  __func__);

      } /* End of switch on type */
      break;

    case CS_ENFORCEMENT_SELECTION_CELLS:
      {
        const cs_adjacency_t  *c2v = connect->c2v;

        switch (efp->type) {

        case CS_ENFORCEMENT_BY_CONSTANT:
          for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

            const cs_lnum_t  c_id = efp->elt_ids[i];
            for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
              for (int k = 0; k < stride; k++)
                values[stride*c2v->ids[j] + k] = efp->values[k];

          }
          break;

        case CS_ENFORCEMENT_BY_DOF_VALUES:
          for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

            const cs_lnum_t  c_id = efp->elt_ids[i];
            for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
              for (int k = 0; k < stride; k++)
                values[stride*c2v->ids[j] + k] = efp->values[stride*c_id + k];

          }
          break;

        default:
          bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of definition.\n",
                    __func__);

        } /* End of switch on type */

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of selection",
                __func__);

    } /* End of switch on selection */

  } /* Loop on parameters */

  /* Parallel synchronization: If there is a conflict between two definitions
     or if a DoF at the parallel interface is defined and its corresponding one
     is not defined, one takes the min. values (since one initializes with
     FLT_MAX). */

  if (connect->vtx_ifs != NULL)
    cs_interface_set_min(connect->vtx_ifs,
                         n_vertices,  /* array size */
                         stride,      /* array stride */
                         true,        /* interlace */
                         CS_REAL_TYPE,
                         values);

  return values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a cs_enforcement_t structure for face-based scheme
 *
 * \param[in] connect    pointer to a cs_cdo_connect_t
 * \param[in] n_params   number of enforcement parameters
 * \param[in] efp_array  array of parameter structures defining the enforcement
 *
 * \return an array with the values to enforce
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_enforcement_define_at_faces(const cs_cdo_connect_t     *connect,
                               int                         n_params,
                               cs_enforcement_param_t    **efp_array)
{
  if (n_params == 0)
    return NULL;

  cs_lnum_t  n_faces = connect->n_faces[0];
  cs_real_t  *values = NULL;
  int  stride = (efp_array[0])->stride;

  BFT_MALLOC(values, stride*n_faces, cs_real_t);
  for (cs_lnum_t i = 0; i < stride*n_faces; i++)
    values[i] = FLT_MAX; /* By default, max value */

  /* Define the enforcement value for each vertex related to an enforcement */

  for (int param_id = 0; param_id < n_params; param_id++) {

    cs_enforcement_param_t  *efp = efp_array[param_id];
    assert(stride == efp->stride);

    switch (efp->selection_type) {

    case CS_ENFORCEMENT_SELECTION_FACES:

      switch (efp->type) {

      case CS_ENFORCEMENT_BY_CONSTANT:
        for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

          cs_lnum_t  f_id = efp->elt_ids[i];
          for (int k = 0; k < stride; k++)
            values[stride*f_id + k] = efp->values[k];

        }
        break;

      case CS_ENFORCEMENT_BY_DOF_VALUES:
        for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

          cs_lnum_t  f_id = efp->elt_ids[i];
          for (int k = 0; k < stride; k++)
            values[stride*f_id + k] = efp->values[stride*i+k];

        }
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of definition.\n",
                  __func__);

      } /* End of switch on type */
      break;

    case CS_ENFORCEMENT_SELECTION_CELLS:
      {
        const cs_adjacency_t  *c2f = connect->c2f;

        switch (efp->type) {

        case CS_ENFORCEMENT_BY_CONSTANT:
          for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

            const cs_lnum_t  c_id = efp->elt_ids[i];
            for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++)
              for (int k = 0; k < stride; k++)
                values[stride*c2f->ids[j] + k] = efp->values[k];

          }
          break;

        case CS_ENFORCEMENT_BY_DOF_VALUES:
          for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

            const cs_lnum_t  c_id = efp->elt_ids[i];
            for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++)
              for (int k = 0; k < stride; k++)
                values[stride*c2f->ids[j] + k] = efp->values[stride*c_id + k];

          }
          break;

        default:
          bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of definition.\n",
                    __func__);

        } /* End of switch on type */

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of selection",
                __func__);

    } /* End of switch on selection */

  } /* Loop on parameters */

  /* Parallel synchronization: If there is a conflict between two definitions
     or if a DoF at the parallel interface is defined and its corresponding one
     is not defined, one takes the min. values (since one initializes with
     FLT_MAX). */

  if (connect->face_ifs != NULL)
    cs_interface_set_min(connect->face_ifs,
                         n_faces,     /* array size */
                         stride,      /* array stride */
                         true,        /* interlace */
                         CS_REAL_TYPE,
                         values);

  return values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a cs_enforcement_t structure for edge-based scheme
 *
 * \param[in] connect    pointer to a cs_cdo_connect_t
 * \param[in] n_params   number of enforcement parameters
 * \param[in] efp_array  array of parameter structures defining the enforcement
 *
 * \return an array with the values to enforce
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_enforcement_define_at_edges(const cs_cdo_connect_t     *connect,
                               int                         n_params,
                               cs_enforcement_param_t    **efp_array)
{
  if (n_params == 0)
    return NULL;

  cs_lnum_t  n_edges = connect->n_edges;
  cs_real_t  *values = NULL;
  int  stride = (efp_array[0])->stride;

  BFT_MALLOC(values, stride*n_edges, cs_real_t);
  for (cs_lnum_t i = 0; i < stride*n_edges; i++)
    values[i] = FLT_MAX; /* By default, max value */

  /* Define the enforcement value for each vertex related to an enforcement */

  for (int param_id = 0; param_id < n_params; param_id++) {

    cs_enforcement_param_t  *efp = efp_array[param_id];
    assert(stride == efp->stride);

    switch (efp->selection_type) {

    case CS_ENFORCEMENT_SELECTION_EDGES:

      switch (efp->type) {

      case CS_ENFORCEMENT_BY_CONSTANT:
        for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

          cs_lnum_t  f_id = efp->elt_ids[i];
          for (int k = 0; k < stride; k++)
            values[stride*f_id + k] = efp->values[k];

        }
        break;

      case CS_ENFORCEMENT_BY_DOF_VALUES:
        for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

          cs_lnum_t  f_id = efp->elt_ids[i];
          for (int k = 0; k < stride; k++)
            values[stride*f_id + k] = efp->values[stride*i+k];

        }
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of definition.\n",
                  __func__);

      } /* End of switch on type */
      break;

    case CS_ENFORCEMENT_SELECTION_CELLS:
      {
        const cs_adjacency_t  *c2e = connect->c2e;

        switch (efp->type) {

        case CS_ENFORCEMENT_BY_CONSTANT:
          for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

            const cs_lnum_t  c_id = efp->elt_ids[i];
            for (cs_lnum_t j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++)
              for (int k = 0; k < stride; k++)
                values[stride*c2e->ids[j] + k] = efp->values[k];

          }
          break;

        case CS_ENFORCEMENT_BY_DOF_VALUES:
          for (cs_lnum_t i = 0; i < efp->n_elts; i++) {

            const cs_lnum_t  c_id = efp->elt_ids[i];
            for (cs_lnum_t j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++)
              for (int k = 0; k < stride; k++)
                values[stride*c2e->ids[j] + k] = efp->values[stride*c_id + k];

          }
          break;

        default:
          bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of definition.\n",
                    __func__);

        } /* End of switch on type */

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid type of selection",
                __func__);

    } /* End of switch on selection */

  } /* Loop on parameters */

  /* Parallel synchronization: If there is a conflict between two definitions
     or if a DoF at the parallel interface is defined and its corresponding one
     is not defined, one takes the min. values (since one initializes with
     FLT_MAX). */

  if (connect->edge_ifs != NULL)
    cs_interface_set_min(connect->edge_ifs,
                         n_edges,     /* array size */
                         stride,      /* array stride */
                         true,        /* interlace */
                         CS_REAL_TYPE,
                         values);

  return values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the cell-wise value to enforce
 *
 * \param[in]      forced_values     values to enforce or FLT_MAX
 * \param[in, out] csys              pointer to a cs_cell_sys_t structure
 * \param[in, out] cw_forced_values  local values to enforce
 *
 * \return true if at least one DoF has to be enforced
 */
/*----------------------------------------------------------------------------*/

bool
cs_enforcement_dofs_cw(const cs_real_t      *forced_values,
                       cs_cell_sys_t        *csys,
                       cs_real_t            *cw_forced_values)
{
  assert(forced_values != NULL && cw_forced_values != NULL);

  /* cw_forced_values is initialized by the calling function */

  /* Initialization */

  bool  has_enforcement = false;

  for (int i = 0; i < csys->n_dofs; i++) {

    csys->dof_is_forced[i] = false;    /* Not enforced by default */

    /* In case of a Dirichlet BC or a circulation BC, this BC is applied and
       the enforcement is ignored. A dirichlet BC is included inside the
       circulation condition */

    if (!cs_cdo_bc_is_circulation(csys->dof_flag[i])) {

      cs_real_t  _val = forced_values[csys->dof_ids[i]];
      if (_val < FLT_MAX) {

        cw_forced_values[i] = _val;
        has_enforcement = true;
        csys->dof_is_forced[i] = true;

      }

    } /* DoF is not associated to a Dirichlet BCs */

  } /* Loop on cell-wise DoFs */

  return has_enforcement;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
