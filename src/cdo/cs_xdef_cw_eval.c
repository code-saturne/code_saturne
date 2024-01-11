/*============================================================================
 * Manage the (generic) cellwise evaluation of extended definitions
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_defs.h"
#include "cs_field.h"
#include "cs_mesh_location.h"
#include "cs_reco.h"
#include "cs_reco_cw.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_xdef_cw_eval.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definition (unset at the end of file)
 *============================================================================*/

/* Redefined the name of functions from cs_math to get shorter names */

#define _dp3  cs_math_3_dot_product

/*=============================================================================
 * Local variables
 *============================================================================*/

static const char _err_empty_array[] =
  " %s: Array storing the evaluation should be allocated before the call"
  " to this function.";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Integrate an analytic function over a face
 *
 * \param[in]      cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      t_eval   time at which the function is evaluated
 * \param[in]      f        face id in the local cell numbering
 * \param[in]      ana      analytic function to integrate
 * \param[in]      input    pointer to an input structure
 * \param[in]      qfunc    quadrature function to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_f_int_by_analytic(const cs_cell_mesh_t            *cm,
                                  double                           t_eval,
                                  short int                        f,
                                  cs_analytic_func_t              *ana,
                                  void                            *input,
                                  cs_quadrature_tria_integral_t   *qfunc,
                                  cs_real_t                       *eval)
{
  const cs_quant_t  pfq = cm->face[f];
  const int  start = cm->f2e_idx[f];
  const int  end = cm->f2e_idx[f+1];
  const short int n_vf = end - start;  /* #vertices (=#edges) */
  const short int *f2e_ids = cm->f2e_ids + start;

  switch (n_vf) {
  case CS_TRIANGLE_CASE:
    {
      short int  v0, v1, v2;
      cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

      qfunc(t_eval, cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2, pfq.meas,
            ana, input, eval);
    }
    break;

  default:
    {
      const double *tef = cm->tef + start;
      for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

        const short int  *e2v = cm->e2v_ids + 2*f2e_ids[e];
        const double  *xv0 = cm->xv + 3*e2v[0];
        const double  *xv1 = cm->xv + 3*e2v[1];

        qfunc(t_eval, xv0, xv1, pfq.center, tef[e], ana, input, eval);

      }
    }
  } /* Switch */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Integrate an analytic function over a cell
 *
 * \param[in]      cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      t_eval   time at which the function is evaluated
 * \param[in]      ana      analytic function to integrate
 * \param[in]      input    pointer to an input structure
 * \param[in]      qfunc    quadrature function to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_c_int_by_analytic(const cs_cell_mesh_t            *cm,
                                  double                           t_eval,
                                  cs_analytic_func_t              *ana,
                                  void                            *input,
                                  cs_quadrature_tetra_integral_t  *qfunc,
                                  cs_real_t                       *eval)
{
  switch (cm->type) {

  case FVM_CELL_TETRA:
    {
      assert(cm->n_fc == 4 && cm->n_vc == 4);
      qfunc(t_eval, cm->xv, cm->xv+3, cm->xv+6, cm->xv+9, cm->vol_c,
            ana, input, eval);
    }
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:
    {
      for (short int f = 0; f < cm->n_fc; ++f) {

        const cs_quant_t  pfq = cm->face[f];
        const double  hf_coef = cs_math_1ov3 * cm->hfc[f];
        const int  start = cm->f2e_idx[f];
        const int  end = cm->f2e_idx[f+1];
        const short int  n_vf = end - start; /* #vertices (=#edges) */
        const short int  *f2e_ids = cm->f2e_ids + start;

        assert(n_vf > 2);
        switch(n_vf){

        case CS_TRIANGLE_CASE: /* Optimized version, no subdivision */
          {
            short int  v0, v1, v2;
            cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids,
                                             &v0, &v1, &v2);

            const double  *xv0 = cm->xv + 3*v0;
            const double  *xv1 = cm->xv + 3*v1;
            const double  *xv2 = cm->xv + 3*v2;

            qfunc(t_eval, xv0, xv1, xv2, cm->xc, hf_coef * pfq.meas,
                  ana, input, eval);
          }
          break;

        default:
          {
            const double  *tef = cm->tef + start;

            for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

              const short int  *e2v = cm->e2v_ids + 2*f2e_ids[e];
              const double  *xv0 = cm->xv + 3*e2v[0];
              const double  *xv1 = cm->xv + 3*e2v[1];

              qfunc(t_eval, xv0, xv1, pfq.center, cm->xc, hf_coef * tef[e],
                    ana, input, eval);
            }
          }
          break;

        } /* End of switch */
      } /* End of loop on faces */

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,  _(" Unknown cell-type.\n"));
    break;

  } /* End of switch on the cell-type */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Routine to integrate an analytic function over a cell and its faces
 *
 * \param[in]  cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in]  t_eval   physical time at which one evaluates the term
 * \param[in]  ana      analytic function to integrate
 * \param[in]  input    pointer to an input structure
 * \param[in]  dim      dimension of the function
 * \param[in]  q_tet    quadrature function to use on tetrahedra
 * \param[in]  q_tri    quadrature function to use on triangles
 * \param[out] c_int    result of the evaluation on the cell
 * \param[out] f_int    result of the evaluation on the faces
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_fc_int_by_analytic(const cs_cell_mesh_t             *cm,
                                   cs_real_t                         t_eval,
                                   cs_analytic_func_t               *ana,
                                   void                             *input,
                                   const short int                   dim,
                                   cs_quadrature_tetra_integral_t   *q_tet,
                                   cs_quadrature_tria_integral_t    *q_tri,
                                   cs_real_t                        *c_int,
                                   cs_real_t                        *f_int)
{
  short int  v0, v1, v2;

  const short int  nf = cm->n_fc;

  /* Switch on the cell-type: optimised version for tetra */

  switch (cm->type) {

  case FVM_CELL_TETRA:
    {
      assert(cm->n_fc == 4 && cm->n_vc == 4);
      q_tet(t_eval, cm->xv, cm->xv+3, cm->xv+6, cm->xv+9, cm->vol_c,
            ana, input, c_int);

      for (short int f = 0; f < nf; ++f) {

        const cs_quant_t  pfq = cm->face[f];
        const short int  *f2e_ids = cm->f2e_ids + cm->f2e_idx[f];

        cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);
        q_tri(t_eval, cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2, pfq.meas,
              ana, input, f_int + dim*f);
      }
    }
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:
  {
    for (short int f = 0; f < nf; ++f) {

      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = cs_math_1ov3 * cm->hfc[f];
      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int  n_vf = end - start; /* #vertices (=#edges) */
      const short int  *f2e_ids = cm->f2e_ids + start;

      assert(n_vf > 2);
      switch(n_vf){

      case CS_TRIANGLE_CASE: /* triangle (optimized version, no subdivision) */
        {
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          const double  *xv0 = cm->xv + 3*v0;
          const double  *xv1 = cm->xv + 3*v1;
          const double  *xv2 = cm->xv + 3*v2;

          q_tet(t_eval, xv0, xv1, xv2, cm->xc, hf_coef * pfq.meas,
                ana, input, c_int);
          q_tri(t_eval, xv0, xv1, xv2, pfq.meas, ana, input, f_int + dim*f);
        }
        break;

      default:
        {
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            const short int  *e2v = cm->e2v_ids + 2*f2e_ids[e];
            const double  *xv0 = cm->xv + 3*e2v[0];
            const double  *xv1 = cm->xv + 3*e2v[1];

            q_tet(t_eval, xv0, xv1, pfq.center, cm->xc, hf_coef*tef[e],
                  ana, input, c_int);
            q_tri(t_eval, xv0, xv1, pfq.center, tef[e],
                  ana, input, f_int + dim*f);

          }
        }
        break;

      } /* End of switch */

    } /* End of loop on faces */

  }
  break;

  default:
    bft_error(__FILE__, __LINE__, 0,  _(" Unknown cell-type.\n"));
    break;

  } /* End of switch on the cell-type */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for evaluating the average on a face of a scalar
 *        function defined through a descriptor (\ref cs_xdef_t structure) by a
 *        cellwise process (usage of a \ref cs_cell_mesh_t structure)
 *
 * \param[in]      cm          pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      f           local face id
 * \param[in]      time_eval   physical time at which one evaluates the term
 * \param[in]      context     pointer to a context structure
 * \param[in]      qtype       level of quadrature to use
 * \param[in, out] eval        result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_scalar_face_avg_by_analytic(const cs_cell_mesh_t   *cm,
                                            short int               f,
                                            cs_real_t               time_eval,
                                            void                   *context,
                                            cs_quadrature_type_t    qtype,
                                            cs_real_t              *eval)
{
  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(context != NULL);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  cs_quadrature_tria_integral_t
    *qfunc = cs_quadrature_get_tria_integral(1, qtype);
  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;

  cs_xdef_cw_eval_f_int_by_analytic(cm, time_eval, f,
                                    ac->func, ac->input, qfunc, eval);

  /* Average */

  eval[0] /= cm->face[f].meas;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a vector
 *         function defined through a descriptor (\ref cs_xdef_t structure) by
 *         a cellwise process (usage of a \ref cs_cell_mesh_t structure)
 *
 * \param[in]      cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      f        local face id
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      context  pointer to a context structure
 * \param[in]      qtype    level of quadrature to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_vector_face_avg_by_analytic(const cs_cell_mesh_t    *cm,
                                            short int                f,
                                            cs_real_t                t_eval,
                                            void                    *context,
                                            cs_quadrature_type_t     qtype,
                                            cs_real_t               *eval)
{
  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(context != NULL);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  cs_quadrature_tria_integral_t
    *qfunc = cs_quadrature_get_tria_integral(3, qtype);
  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;

  cs_xdef_cw_eval_f_int_by_analytic(cm, t_eval, f,
                                    ac->func, ac->input, qfunc, eval);

  /* Average */

  const double _oversurf = 1./cm->face[f].meas;
  for (short int xyz = 0; xyz < 3; xyz++)
    eval[xyz] *= _oversurf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the average on a face of a tensor
 *         function defined through a descriptor (\ref cs_xdef_t structure) by
 *         a cellwise process (usage of a \ref cs_cell_mesh_t structure)
 *
 * \param[in]      cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      f        local face id
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      context  pointer to a context structure
 * \param[in]      qtype    level of quadrature to use
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_tensor_face_avg_by_analytic(const cs_cell_mesh_t    *cm,
                                            short int                f,
                                            cs_real_t                t_eval,
                                            void                    *context,
                                            cs_quadrature_type_t     qtype,
                                            cs_real_t               *eval)
{
  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(context != NULL);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  cs_quadrature_tria_integral_t
    *qfunc = cs_quadrature_get_tria_integral(9, qtype);;
  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;

  cs_xdef_cw_eval_f_int_by_analytic(cm, t_eval, f,
                                    ac->func, ac->input, qfunc, eval);

  /* Average */

  const double _oversurf = 1./cm->face[f].meas;
  for (short int i = 0; i < 9; i++)
    eval[i] *= _oversurf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (\ref cs_xdef_t structure) by a cellwise process (usage
 *         of a \ref cs_cell_mesh_t structure).
 *         This evaluation hinges on the computation of integrals
 *
 * \param[in]      cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      qtype    quadrature type
 * \param[in]      context  pointer to a context structure
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_scalar_avg_by_analytic(const cs_cell_mesh_t     *cm,
                                       cs_real_t                 t_eval,
                                       void                     *context,
                                       cs_quadrature_type_t      qtype,
                                       cs_real_t                *eval)
{
  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(context != NULL);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  cs_quadrature_tetra_integral_t
    *qfunc = cs_quadrature_get_tetra_integral(1, qtype);
  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;

  cs_xdef_cw_eval_c_int_by_analytic(cm, t_eval,
                                    ac->func, ac->input, qfunc,
                                    eval);

  /* Average */

  eval[0] /= cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (\ref cs_xdef_t structure) by a cellwise process (usage
 *         of a \ref cs_cell_mesh_t structure).
 *         This evaluation hinges on the computation of integrals
 *         Vector-valued case.
 *
 * \param[in]      cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      t_eval    physical time at which one evaluates the term
 * \param[in]      qtype     quadrature type
 * \param[in]      context   pointer to a context structure
 * \param[in, out] eval      result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_vector_avg_by_analytic(const cs_cell_mesh_t     *cm,
                                       cs_real_t                 t_eval,
                                       void                     *context,
                                       cs_quadrature_type_t      qtype,
                                       cs_real_t                *eval)
{
  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(context != NULL);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  cs_quadrature_tetra_integral_t
    *qfunc = cs_quadrature_get_tetra_integral(3, qtype);
  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;

  cs_xdef_cw_eval_c_int_by_analytic(cm, t_eval,
                                    ac->func, ac->input,
                                    qfunc, eval);

  /* Average */

  const double _overvol = 1./cm->vol_c;
  for (short int xyz = 0; xyz < 3; xyz++)
    eval[xyz] *= _overvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined through a
 *         descriptor (\ref cs_xdef_t structure) by a cellwise process (usage
 *         of a \ref cs_cell_mesh_t structure).
 *         This evaluation hinges on the computation of integrals
 *         Tensor-valued case.
 *
 * \param[in]      cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      t_eval    physical time at which one evaluates the term
 * \param[in]      qtype     quadrature type
 * \param[in]      context   pointer to a context structure
 * \param[in, out] eval      result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_tensor_avg_by_analytic(const cs_cell_mesh_t     *cm,
                                       cs_real_t                 t_eval,
                                       void                     *context,
                                       cs_quadrature_type_t      qtype,
                                       cs_real_t                *eval)
{
  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(context != NULL);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  cs_quadrature_tetra_integral_t
    *qfunc = cs_quadrature_get_tetra_integral(9, qtype);
  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;

  cs_xdef_cw_eval_c_int_by_analytic(cm, t_eval,
                                    ac->func, ac->input, qfunc,
                                    eval);

  /* Average */

  const double _overvol = 1./cm->vol_c;
  for (short int xyz = 0; xyz < 9; xyz++)
    eval[xyz] *= _overvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a  quantity by a cellwise process using a definition by
 *         time function
 *
 * \param[in]  cm          pointer to a \ref cs_cell_mesh_t structure
 * \param[in]  time_eval   physical time at which one evaluates the term
 * \param[in]  context     pointer to a context structure
 * \param[out] eval        result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_by_time_func(const cs_cell_mesh_t     *cm,
                             cs_real_t                 time_eval,
                             void                     *context,
                             cs_real_t                *eval)
{
  CS_UNUSED(cm);
  cs_xdef_time_func_context_t *tfc = (cs_xdef_time_func_context_t *)context;

  /* Evaluate the quantity */
  tfc->func(time_eval, tfc->input, eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity at the cell center defined using an analytic
 *         function by a cellwise process (usage of a \ref cs_cell_mesh_t
 *         structure)
 *
 * \param[in]      cm          pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      time_eval   physical time at which one evaluates the term
 * \param[in]      context     pointer to a context structure
 * \param[in, out] eval        result of the evaluation at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_by_analytic(const cs_cell_mesh_t       *cm,
                            cs_real_t                   time_eval,
                            void                       *context,
                            cs_real_t                  *eval)
{
  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;

  /* Evaluate the function for this time at the cell center */

  ac->func(time_eval,
           1, NULL, cm->xc,
           true, /* compacted output ? */
           ac->input,
           eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity at cells defined by an array.
 *         Array is assumed to be interlaced.
 *         Variation using a \ref cs_cell_mesh_t structure
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in, out] eval       result of the evaluation at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_by_array(const cs_cell_mesh_t      *cm,
                         cs_real_t                  time_eval,
                         void                      *context,
                         cs_real_t                 *eval)
{
  CS_UNUSED(time_eval);

  cs_xdef_array_context_t  *cx = (cs_xdef_array_context_t *)context;

  const int  stride = cx->stride;

  /* Array is assumed to be interlaced */

  if (cs_flag_test(cx->value_location, cs_flag_primal_cell)) {
    /*                                 =================== */

    if (cx->full_length)
      for (int k = 0; k < stride; k++)
        eval[k] = cx->values[stride*cm->c_id + k];

    else {

      assert(cx->full2subset != NULL);
      cs_lnum_t  compact_id = cx->full2subset[cm->c_id];
      assert(compact_id > -1);

      for (int k = 0; k < stride; k++)
        eval[k] = cx->values[stride*compact_id + k];

    }

  }
  else if (cs_flag_test(cx->value_location, cs_flag_primal_vtx)) {
    /*                                      ================== */

    assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));
    assert(cx->full_length == true);

    /* Reconstruct (or interpolate) value at the current cell center */

    if (stride == 1)
      eval[0] = cs_reco_cw_scalar_v2c(cm, cx->values);
    else
      cs_reco_cw_stride_v2c(stride, cm, cx->values, eval);

  }
  else if (cs_flag_test(cx->value_location, cs_flag_dual_cell_byc)) {
    /*                                      ===================== */

    const cs_adjacency_t  *adj = cx->adjacency;
    assert(adj != NULL);
    assert(cx->full_length == true);

    /* Reconstruct (or interpolate) value at the current cell center */

    if (stride == 1)
      eval[0] = cs_reco_cw_scalar_vbyc2c(cm, cx->values + adj->idx[cm->c_id]);
    else
      cs_reco_cw_stride_vbyc2c(stride, cm, cx->values + adj->idx[cm->c_id],
                               eval);

  }
  else if (cs_flag_test(cx->value_location, cs_flag_primal_edge_byc)) {
    /*                                      ======================= */

    const cs_adjacency_t  *adj = cx->adjacency;
    assert(adj != NULL);
    assert(cx->full_length == true);
    assert(stride == 1);

    /* Reconstruct (or interpolate) value at the current cell center */

    eval[0] = cs_reco_cw_scalar_ebyc2c(cm, cx->values + adj->idx[cm->c_id]);

  }
  else if (cs_flag_test(cx->value_location, cs_flag_dual_face_byc)) {
    /*                                      ===================== */

    const cs_adjacency_t  *adj = cx->adjacency;
    assert(adj != NULL);
    assert(cx->full_length == true);

    /* Reconstruct (or interpolate) value at the current cell center */

    cs_reco_dfbyc_in_cell(cm,
                          cx->values + adj->idx[cm->c_id],
                          eval);

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid support for the input array", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a quantity inside a cell defined using a field
 *         Variation using a \ref cs_cell_mesh_t structure
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in, out] eval       value of the property at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_by_field(const cs_cell_mesh_t      *cm,
                         cs_real_t                  time_eval,
                         void                      *context,
                         cs_real_t                 *eval)
{
  CS_UNUSED(time_eval);

  cs_field_t  *field = (cs_field_t *)context;
  assert(field != NULL);
  cs_real_t  *values = field->val;
  const int  c_ml_id = cs_mesh_location_get_id_by_name(N_("cells"));
  const int  v_ml_id = cs_mesh_location_get_id_by_name(N_("vertices"));

  if (field->location_id == c_ml_id) {

    for (int k = 0; k < field->dim; k++)
      eval[k] = values[field->dim*cm->c_id + k];

  }
  else if (field->location_id == v_ml_id) {

    assert(field->dim == 1);
    assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

    /* Reconstruct (or interpolate) value at the current cell center */

    for (short int v = 0; v < cm->n_vc; v++)
      eval[0] += cm->wvc[v] * values[cm->v_ids[v]];


  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid support for the input array", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined by analytic
 *         function at a precise location (x, y, z) inside a cell
 *         Use of a \ref cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      n_points   number of points where to compute the evaluation
 * \param[in]      xyz        where to compute the evaluation
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in, out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_at_xyz_by_analytic(const cs_cell_mesh_t       *cm,
                                   cs_lnum_t                   n_points,
                                   const cs_real_t            *xyz,
                                   cs_real_t                   time_eval,
                                   void                       *context,
                                   cs_real_t                  *eval)
{
  CS_UNUSED(cm);

  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;

  /* Evaluate the function for this time at the given coordinates */

  ac->func(time_eval,
           n_points, NULL, xyz, true, /* compacted output ? */
           ac->input,
           eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined by analytic
 *         function at a precise location inside a cell
 *         Use of a \ref cs_cell_mesh_t structure.
 *         Vector-valued case.
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      n_points   number of points where to compute the evaluation
 * \param[in]      xyz        where to compute the evaluation
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in, out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_vector_at_xyz_by_array(const cs_cell_mesh_t     *cm,
                                       cs_lnum_t                 n_points,
                                       const cs_real_t          *xyz,
                                       cs_real_t                 time_eval,
                                       void                     *context,
                                       cs_real_t                *eval)
{
  CS_UNUSED(xyz);
  CS_UNUSED(time_eval);

  cs_xdef_array_context_t  *cx = (cs_xdef_array_context_t *)context;

  const int  stride = cx->stride;

  /* Array is assumed to be interlaced */

  if (cs_flag_test(cx->value_location, cs_flag_primal_cell)) {

    assert(stride == 3);
    assert(stride == n_points);

    cs_real_3_t  cell_vector;
    const cs_real_t  *_val = cx->values + stride*cm->c_id;
    for (int k = 0; k < stride; k++)
      cell_vector[k] = _val[k];

    for (int i = 0; i < n_points; i++) {
      eval[3*i    ] = cell_vector[0];
      eval[3*i + 1] = cell_vector[1];
      eval[2*i + 2] = cell_vector[2];
    }

  }
  else if (cs_flag_test(cx->value_location, cs_flag_primal_vtx)) {

    assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));
    assert(stride == 3);

    /* Reconstruct (or interpolate) value at the current cell center */

    for (int k = 0; k < stride; k++)
      for (short int v = 0; v < cm->n_vc; v++)
        eval[k] += cm->wvc[v] * cx->values[stride*cm->v_ids[v] + k];

  }
  else if (cs_flag_test(cx->value_location, cs_flag_dual_face_byc)) {

    const cs_adjacency_t  *adj = cx->adjacency;
    assert(adj != NULL);

    /* Reconstruct (or interpolate) value at the current cell center */

    cs_real_3_t  cell_vector;
    cs_reco_dfbyc_in_cell(cm,
                          cx->values + adj->idx[cm->c_id],
                          cell_vector);

    for (int i = 0; i < n_points; i++) {
      eval[3*i    ] = cell_vector[0];
      eval[3*i + 1] = cell_vector[1];
      eval[3*i + 2] = cell_vector[2];
    }

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid support for the input array", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating a quantity defined by a field
 *         at a precise location inside a cell
 *         Use of a \ref cs_cell_mesh_t structure.
 *         Vector-valued case.
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      n_points   number of points where to compute the evaluation
 * \param[in]      xyz        where to compute the evaluation
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in, out] eval       result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_vector_at_xyz_by_field(const cs_cell_mesh_t    *cm,
                                       cs_lnum_t                n_points,
                                       const cs_real_t         *xyz,
                                       cs_real_t                time_eval,
                                       void                    *context,
                                       cs_real_t               *eval)
{
  CS_UNUSED(xyz);
  CS_UNUSED(time_eval);

  cs_field_t  *field = (cs_field_t *)context;
  const cs_real_t  *values = field->val;

  assert(field != NULL);
  assert(field->dim == 3);

  const int  c_ml_id = cs_mesh_location_get_id_by_name(N_("cells"));
  const int  v_ml_id = cs_mesh_location_get_id_by_name(N_("vertices"));

  /* Array is assumed to be interlaced */
  if (field->location_id == c_ml_id) {

    cs_real_3_t  cell_vector;
    for (int k = 0; k < 3; k++)
      cell_vector[k] = values[3*cm->c_id + k];

    for (int i = 0; i < n_points; i++) { /* No interpolation */
      eval[3*i    ] = cell_vector[0];
      eval[3*i + 1] = cell_vector[1];
      eval[3*i + 2] = cell_vector[2];
    }

  }
  else if (field->location_id == v_ml_id) {

    assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

    /* Reconstruct (or interpolate) value at the current cell center */

    for (int k = 0; k < 3; k++)
      for (short int v = 0; v < cm->n_vc; v++)
        eval[k] += cm->wvc[v] * values[3*cm->v_ids[v] + k];

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid support for the input array", __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux defined by
 *         scalar-valued quantities. The normal flux is then added to each
 *         portion of face related to a vertex.
 *         Use of a \ref cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in, out] eval       array of Neumann fluxes at DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_flux_v_by_scalar_val(const cs_cell_mesh_t     *cm,
                                     short int                 f,
                                     cs_real_t                 time_eval,
                                     void                     *context,
                                     cs_real_t                *eval)
{
  CS_UNUSED(time_eval);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_EV | CS_FLAG_COMP_FE));

  const cs_real_t  *flux = (cs_real_t *)context;

  /* Scaled by 1/2 since the part for a vertex v from the triangle tef is
     weighted by 1/2 */

  const cs_real_t  half_full_flux = 0.5 * flux[0];

  if (cs_eflag_test(cm->flag, CS_FLAG_COMP_FEQ)) {

    /* Loop on face edges */

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const double  _flx = cm->tef[i] * half_full_flux;
      const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];

      eval[e2v[0]] += _flx;
      eval[e2v[1]] += _flx;

    }

  }
  else {

    /* Loop on face edges */

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const double  tef = cs_compute_area_from_quant(cm->edge[e],
                                                     cm->face[f].center);
      const double  _flx = tef * half_full_flux;

      eval[cm->e2v_ids[2*e  ]] += _flx;
      eval[cm->e2v_ids[2*e+1]] += _flx;

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux defined by
 *         vector-valued quantities. The normal flux is then added to each
 *         portion of face related to a vertex.
 *         Use of a \ref cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in, out] eval       array of Neumann fluxes at DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_flux_v_by_vector_val(const cs_cell_mesh_t     *cm,
                                     short int                 f,
                                     cs_real_t                 time_eval,
                                     void                     *context,
                                     cs_real_t                *eval)
{
  CS_UNUSED(time_eval);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_EV | CS_FLAG_COMP_FE));

  const cs_real_t  *flux = (cs_real_t *)context;

  /* Scaled by 1/2 since the part for a vertex v from the triangle tef is
     weighted by 1/2 */

  const cs_real_t  half_full_flux = 0.5 * _dp3(flux, cm->face[f].unitv);

  if (cs_eflag_test(cm->flag, CS_FLAG_COMP_FEQ)) {

    /* Loop on face edges */

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const double  _flx = cm->tef[i] * half_full_flux;
      const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];

      eval[e2v[0]] += _flx;
      eval[e2v[1]] += _flx;

    }

  }
  else {

    /* Loop on face edges */

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const double  tef = cs_compute_area_from_quant(cm->edge[e],
                                                     cm->face[f].center);
      const double  _flx = tef * half_full_flux;

      eval[cm->e2v_ids[2*e  ]] += _flx;
      eval[cm->e2v_ids[2*e+1]] += _flx;

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux defined by a
 *         scalar-valued quantities and relying on an analytic function. The
 *         normal flux is then added to each portion of face related to a
 *         vertex.
 *         Use of a \ref cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in]      qtype      level of quadrature to use
 * \param[in, out] eval       array of Neumann fluxes at DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_flux_v_by_scalar_analytic(const cs_cell_mesh_t      *cm,
                                          short int                  f,
                                          cs_real_t                  time_eval,
                                          void                      *context,
                                          cs_quadrature_type_t       qtype,
                                          cs_real_t                 *eval)
{
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PFQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_FE));

  const cs_xdef_analytic_context_t *ac = (cs_xdef_analytic_context_t *)context;
  const cs_quant_t  fq = cm->face[f];

  switch (qtype) {

  case CS_QUADRATURE_NONE:
  case CS_QUADRATURE_BARY:
    {
      cs_real_t  flux_xf = 0;

      /* Evaluate the function for this time at the face center */

      ac->func(time_eval, 1, NULL, fq.center, true, /* compacted output ? */
               ac->input,
               &flux_xf);

      /* Plug into the evaluation by value now */

      cs_xdef_cw_eval_flux_v_by_scalar_val(cm, f, time_eval, &flux_xf, eval);
    }
    break;

  case CS_QUADRATURE_BARY_SUBDIV:
    {
      assert(cs_flag_test(cm->flag, CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ));
      cs_real_t  _val[2];
      cs_real_3_t  _xyz[2];

      if (cs_flag_test(cm->flag, CS_FLAG_COMP_FEQ)) {

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];

          for (int k = 0; k < 3; k++) {
            const double xef = cm->edge[e].center[k] + fq.center[k];
            _xyz[0][k] = cs_math_1ov3 * (xef + cm->xv[3*v1+k]);
            _xyz[1][k] = cs_math_1ov3 * (xef + cm->xv[3*v2+k]);
          }

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 2, NULL,
                   (const cs_real_t *)_xyz, true, /* compacted output ? */
                   ac->input,
                   (cs_real_t *)_val);

          eval[v1] += 0.5*cm->tef[i] * _val[0];
          eval[v2] += 0.5*cm->tef[i] * _val[1];

        }
      }
      else {

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];

          for (int k = 0; k < 3; k++) {
            const double xef = cm->edge[e].center[k] + fq.center[k];
            _xyz[0][k] = cs_math_1ov3 * (xef + cm->xv[3*v1+k]);
            _xyz[1][k] = cs_math_1ov3 * (xef + cm->xv[3*v2+k]);
          }

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 2, NULL,
                   (const cs_real_t *)_xyz, true, /* compacted output ? */
                   ac->input,
                   (cs_real_t *)_val);

          const double tef = cs_compute_area_from_quant(cm->edge[e], fq.center);

          eval[v1] += 0.5 * tef * _val[0];
          eval[v2] += 0.5 * tef * _val[1];

        }

      }

    }
    break; /* BARY_SUBDIV */

  case CS_QUADRATURE_HIGHER:
    {
      assert(cs_flag_test(cm->flag, CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ));

      /* Two triangles s_{vef} related to a vertex and four values by triangle
       * --> 2*4 = 8 Gauss points
       * The flux returns by the analytic function is a vector. So the size
       * of _val is 24 = 8*3
       */

      cs_real_t _val[8], w[8];
      cs_real_3_t  gpts[8];

      if (cs_flag_test(cm->flag, CS_FLAG_COMP_FEQ)) { /* tef is pre-computed */

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];
          const cs_real_t  svef = 0.5 * cm->tef[i];

          /* Two triangles composing the portion of face related to a vertex
             Evaluate the field at the four quadrature points */

          cs_quadrature_tria_4pts(cm->edge[e].center, fq.center, cm->xv + 3*v1,
                                  svef,
                                  gpts, w);

          cs_quadrature_tria_4pts(cm->edge[e].center, fq.center, cm->xv + 3*v2,
                                  svef,
                                  gpts + 4, w + 4);

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 8, NULL,
                   (const cs_real_t *)gpts, true, /* compacted output ? */
                   ac->input,
                   _val);

          cs_real_t  add0 = 0, add1 = 0;
          for (int p = 0; p < 4; p++)
            add0 += w[p] * _val[p], add1 += w[4 + p] * _val[4+p];

          eval[v1] += add0;
          eval[v2] += add1;

        }

      }
      else {

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];
          const double svef = 0.5 * cs_compute_area_from_quant(cm->edge[e],
                                                               fq.center);

          /* Two triangles composing the portion of face related to a vertex
             Evaluate the field at the four quadrature points */

          cs_quadrature_tria_4pts(cm->edge[e].center, fq.center, cm->xv + 3*v1,
                                  svef,
                                  gpts, w);

          cs_quadrature_tria_4pts(cm->edge[e].center, fq.center, cm->xv + 3*v2,
                                  svef,
                                  gpts + 4, w + 4);

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 8, NULL,
                   (const cs_real_t *)gpts, true, /* compacted output ? */
                   ac->input,
                   _val);

          cs_real_t  add0 = 0, add1 = 0;
          for (int p = 0; p < 4; p++)
            add0 += w[p] * _val[p], add1 += w[4+p] * _val[4+p];

          eval[v1] += add0;
          eval[v2] += add1;

        }

      } /* Is tef already computed ? */
    }
    break;

  case CS_QUADRATURE_HIGHEST:
    {
      assert(cs_flag_test(cm->flag, CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ));

      /* Two triangles s_{vef} related to a vertex and seven values by triangle
       * --> 2*7 = 14 Gauss points
       * The flux returns by the analytic function is a scalar. So the size
       * of _val is 14
       */

      cs_real_t _val[14], w[14];
      cs_real_3_t  gpts[14];

      if (cs_flag_test(cm->flag, CS_FLAG_COMP_FEQ)) { /* tef is pre-computed */

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];
          const cs_real_t  svef = 0.5 * cm->tef[i];

          /* Two triangles composing the portion of face related to a vertex
             Evaluate the field at the seven quadrature points */

          cs_quadrature_tria_7pts(cm->edge[e].center, fq.center, cm->xv + 3*v1,
                                  svef,
                                  gpts, w);

          cs_quadrature_tria_7pts(cm->edge[e].center, fq.center, cm->xv + 3*v2,
                                  svef,
                                  gpts + 7, w + 7);

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 14, NULL,
                   (const cs_real_t *)gpts, true, /* compacted output ? */
                   ac->input,
                   _val);

          cs_real_t  add0 = 0, add1 = 0;
          for (int p = 0; p < 7; p++)
            add0 += w[p] * _val[p], add1 += w[7+p] * _val[7+p];

          eval[v1] += add0;
          eval[v2] += add1;

        }

      }
      else {

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];
          const double svef = 0.5 * cs_compute_area_from_quant(cm->edge[e],
                                                               fq.center);

          /* Two triangles composing the portion of face related to a vertex
             Evaluate the field at the seven quadrature points */

          cs_quadrature_tria_7pts(cm->edge[e].center, fq.center, cm->xv + 3*v1,
                                  svef,
                                  gpts, w);

          cs_quadrature_tria_7pts(cm->edge[e].center, fq.center, cm->xv + 3*v2,
                                  svef,
                                  gpts + 7, w + 7);

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 14, NULL,
                   (const cs_real_t *)gpts, true, /* compacted output ? */
                   ac->input,
                   _val);

          cs_real_t  add0 = 0, add1 = 0;
          for (int p = 0; p < 7; p++)
            add0 += w[p] * _val[p], add1 += w[7+p] * _val[7+p];

          eval[v1] += add0;
          eval[v2] += add1;

        }

      } /* Is tef already computed ? */

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of quadrature.", __func__);
    break;

  }  /* switch type of quadrature */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux defined by a
 *         vector-valued quantities and relying on an analytic function. The
 *         normal flux is then added to each portion of face related to a
 *         vertex.
 *         Use of a \ref cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in]      qtype      level of quadrature to use
 * \param[in, out] eval       array of Neumann fluxes at DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_flux_v_by_vector_analytic(const cs_cell_mesh_t      *cm,
                                          short int                  f,
                                          cs_real_t                  time_eval,
                                          void                      *context,
                                          cs_quadrature_type_t       qtype,
                                          cs_real_t                 *eval)
{
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PFQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_FE));

  const cs_xdef_analytic_context_t *ac = (cs_xdef_analytic_context_t *)context;
  const cs_quant_t  fq = cm->face[f];

  switch (qtype) {

  case CS_QUADRATURE_NONE:
  case CS_QUADRATURE_BARY:
    {
      cs_real_3_t  flux_xf = {0, 0, 0};

      /* Evaluate the function for this time at the face center */

      ac->func(time_eval, 1, NULL, fq.center, true, /* compacted output ? */
               ac->input,
               flux_xf);

      /* Plug into the evaluation by value now */

      cs_xdef_cw_eval_flux_v_by_vector_val(cm, f, time_eval, flux_xf, eval);
    }
    break;

  case CS_QUADRATURE_BARY_SUBDIV:
    {
      assert(cs_flag_test(cm->flag, CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ));
      cs_real_3_t  _val[2], _xyz[2];

      if (cs_flag_test(cm->flag, CS_FLAG_COMP_FEQ)) {

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];

          for (int k = 0; k < 3; k++) {
            const double xef = cm->edge[e].center[k] + fq.center[k];
            _xyz[0][k] = cs_math_1ov3 * (xef + cm->xv[3*v1+k]);
            _xyz[1][k] = cs_math_1ov3 * (xef + cm->xv[3*v2+k]);
          }

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 2, NULL,
                   (const cs_real_t *)_xyz, true, /* compacted output ? */
                   ac->input,
                   (cs_real_t *)_val);

          eval[v1] += 0.5*cm->tef[i] * _dp3(_val[0], fq.unitv);
          eval[v2] += 0.5*cm->tef[i] * _dp3(_val[1], fq.unitv);

        }

      }
      else {

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];

          for (int k = 0; k < 3; k++) {
            const double xef = cm->edge[e].center[k] + fq.center[k];
            _xyz[0][k] = cs_math_1ov3 * (xef + cm->xv[3*v1+k]);
            _xyz[1][k] = cs_math_1ov3 * (xef + cm->xv[3*v2+k]);
          }

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 2, NULL,
                   (const cs_real_t *)_xyz, true, /* compacted output ? */
                   ac->input,
                   (cs_real_t *)_val);

          const double tef = cs_compute_area_from_quant(cm->edge[e], fq.center);

          eval[v1] += 0.5 * tef * _dp3(_val[0], fq.unitv);
          eval[v2] += 0.5 * tef * _dp3(_val[1], fq.unitv);

        }

      }

    }
    break; /* BARY_SUBDIV */

  case CS_QUADRATURE_HIGHER:
    {
      assert(cs_flag_test(cm->flag, CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ));

      /* Two triangles s_{vef} related to a vertex and four values by triangle
       * --> 2*4 = 8 Gauss points
       * The flux returns by the analytic function is a vector. So the size
       * of _val is 24 = 8*3
       */

      cs_real_t _val[24], w[8];
      cs_real_3_t  gpts[8];

      if (cs_flag_test(cm->flag, CS_FLAG_COMP_FEQ)) { /* tef is pre-computed */

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];
          const cs_real_t  svef = 0.5 * cm->tef[i];

          /* Two triangles composing the portion of face related to a vertex
             Evaluate the field at the four quadrature points */

          cs_quadrature_tria_4pts(cm->edge[e].center, fq.center, cm->xv + 3*v1,
                                  svef,
                                  gpts, w);

          cs_quadrature_tria_4pts(cm->edge[e].center, fq.center, cm->xv + 3*v2,
                                  svef,
                                  gpts + 4, w + 4);

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 8, NULL,
                   (const cs_real_t *)gpts, true, /* compacted output ? */
                   ac->input,
                   _val);

          cs_real_t  add0 = 0, add1 = 0;
          for (int p = 0; p < 4; p++) {
            add0 += w[    p] * _dp3(_val + 3*p,     fq.unitv);
            add1 += w[4 + p] * _dp3(_val + 3*(4+p), fq.unitv);
          }

          eval[v1] += add0;
          eval[v2] += add1;

        }

      }
      else {

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];
          const double svef = 0.5 * cs_compute_area_from_quant(cm->edge[e],
                                                               fq.center);

          /* Two triangles composing the portion of face related to a vertex
             Evaluate the field at the four quadrature points */

          cs_quadrature_tria_4pts(cm->edge[e].center, fq.center, cm->xv + 3*v1,
                                  svef,
                                  gpts, w);

          cs_quadrature_tria_4pts(cm->edge[e].center, fq.center, cm->xv + 3*v2,
                                  svef,
                                  gpts + 4, w + 4);

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 8, NULL,
                   (const cs_real_t *)gpts, true, /* compacted output ? */
                   ac->input,
                   _val);

          cs_real_t  add0 = 0, add1 = 0;
          for (int p = 0; p < 4; p++) {
            add0 += w[    p] * _dp3(_val + 3*p,     fq.unitv);
            add1 += w[4 + p] * _dp3(_val + 3*(4+p), fq.unitv);
          }

          eval[v1] += add0;
          eval[v2] += add1;

        }

      } /* Is tef already computed ? */

    }
    break;

  case CS_QUADRATURE_HIGHEST:
    {
      assert(cs_flag_test(cm->flag, CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ));

      /* Two triangles s_{vef} related to a vertex and seven values by triangle
       * --> 2*7 = 14 Gauss points
       * The flux returns by the analytic function is a vector. So the size
       * of _val is 42 = 14*3
       */

      cs_real_t _val[42], w[14];
      cs_real_3_t  gpts[14];

      if (cs_flag_test(cm->flag, CS_FLAG_COMP_FEQ)) { /* tef is pre-computed */

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];
          const cs_real_t  svef = 0.5 * cm->tef[i];

          /* Two triangles composing the portion of face related to a vertex
             Evaluate the field at the seven quadrature points */

          cs_quadrature_tria_7pts(cm->edge[e].center, fq.center, cm->xv + 3*v1,
                                  svef,
                                  gpts, w);

          cs_quadrature_tria_7pts(cm->edge[e].center, fq.center, cm->xv + 3*v2,
                                  svef,
                                  gpts + 7, w + 7);

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 14, NULL,
                   (const cs_real_t *)gpts, true, /* compacted output ? */
                   ac->input,
                   _val);

          cs_real_t  add0 = 0, add1 = 0;
          for (int p = 0; p < 7; p++) {
            add0 += w[    p] * _dp3(_val + 3*p,     fq.unitv);
            add1 += w[7 + p] * _dp3(_val + 3*(7+p), fq.unitv);
          }

          eval[v1] += add0;
          eval[v2] += add1;

        }

      }
      else {

        /* Loop on face edges */

        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int v1 = cm->e2v_ids[2*e];
          const short int v2 = cm->e2v_ids[2*e+1];
          const double svef = 0.5 * cs_compute_area_from_quant(cm->edge[e],
                                                               fq.center);

          /* Two triangles composing the portion of face related to a vertex
             Evaluate the field at the seven quadrature points */

          cs_quadrature_tria_7pts(cm->edge[e].center, fq.center, cm->xv + 3*v1,
                                  svef,
                                  gpts, w);

          cs_quadrature_tria_7pts(cm->edge[e].center, fq.center, cm->xv + 3*v2,
                                  svef,
                                  gpts + 7, w + 7);

          /* Evaluate the function for this time at the given coordinates */

          ac->func(time_eval, 14, NULL,
                   (const cs_real_t *)gpts, true, /* compacted output ? */
                   ac->input,
                   _val);

          cs_real_t  add0 = 0, add1 = 0;
          for (int p = 0; p < 7; p++) {
            add0 += w[    p] * _dp3(_val + 3*p,     fq.unitv);
            add1 += w[7 + p] * _dp3(_val + 3*(7+p), fq.unitv);
          }

          eval[v1] += add0;
          eval[v2] += add1;

        }

      } /* Is tef already computed ? */

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of quadrature.", __func__);
    break;

  }  /* switch type of quadrature */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal flux of a quantity
 *         defined by analytic function (scalar-valued) at the face f.
 *         Use of a \ref cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in]      qtype      level of quadrature to use
 * \param[in, out] eval       array of values at DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_flux_by_scalar_analytic(const cs_cell_mesh_t      *cm,
                                        short int                  f,
                                        cs_real_t                  time_eval,
                                        void                      *context,
                                        cs_quadrature_type_t       qtype,
                                        cs_real_t                 *eval)
{
  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;

  switch (qtype) {

  case CS_QUADRATURE_NONE:
  case CS_QUADRATURE_BARY:
    {
      double  flux = 0;

      ac->func(time_eval,
               1, NULL, cm->face[f].center,
               true, /* compacted output ? */
               ac->input,
               &flux);

      eval[f] = cm->face[f].meas * flux;
    }
    break;

  case CS_QUADRATURE_BARY_SUBDIV:
    {
      assert(cs_flag_test(cm->flag,
                          CS_FLAG_COMP_EV |CS_FLAG_COMP_FE |CS_FLAG_COMP_FEQ));

      const cs_quant_t  fq = cm->face[f];
      const cs_real_t  *xv = cm->xv;

      eval[f] = 0.; /* Reset value */

      /* Loop on face edges */

      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

        const short int  e = cm->f2e_ids[i];
        const short int v1 = cm->e2v_ids[2*e];
        const short int v2 = cm->e2v_ids[2*e+1];

        cs_real_3_t  _xyz;
        for (int k = 0; k < 3; k++)
          _xyz[k] = cs_math_1ov3 * (fq.center[k] + xv[3*v1+k] + xv[3*v2+k]);

        /* Evaluate the function for this time at the given coordinates */

        double  _eval = 0.;
        ac->func(time_eval, 1, NULL,
                 (const cs_real_t *)_xyz, true, /* compacted output ? */
                 ac->input,
                 &_eval);

        eval[f] += cm->tef[i] * _eval;

      } /* Loop on face edges */
    }
    break; /* BARY_SUBDIV */

  case CS_QUADRATURE_HIGHER:
    {
      assert(cs_flag_test(cm->flag,
                          CS_FLAG_COMP_EV |CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ));

      /* Four values by triangle --> 4 Gauss points
       * The flux returns by the analytic function is a vector. So the size of
       * _eval is 4 */

      cs_real_t  w[4], _eval[4];
      cs_real_3_t  gpts[4];

      const cs_quant_t  fq = cm->face[f];

      eval[f] = 0.; /* Reset value */

      /* Loop on face edges */

      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

        const short int  e = cm->f2e_ids[i];
        const short int v1 = cm->e2v_ids[2*e];
        const short int v2 = cm->e2v_ids[2*e+1];

        /* Evaluate the field at the three quadrature points */

        cs_quadrature_tria_4pts(fq.center, cm->xv + 3*v1, cm->xv + 3*v2,
                                cm->tef[i],
                                gpts, w);

        /* Evaluate the function for this time at the given coordinates */

        ac->func(time_eval, 4, NULL,
                 (const cs_real_t *)gpts, true,  /* compacted output ? */
                 ac->input,
                 _eval);

        cs_real_t  add = 0;
        for (int p = 0; p < 4; p++)
          add += w[p] * _eval[p];

        eval[f] += add;

      } /* Loop on face edges */
    }
    break;


  case CS_QUADRATURE_HIGHEST:
    {
      assert(cs_flag_test(cm->flag,
                          CS_FLAG_COMP_EV |CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ));

      /* Seven values by triangle --> 7 Gauss points
       * The flux returns by the analytic function is a vector. So the size of
       * _eval is 7
       */

      cs_real_t  w[7], _eval[7];
      cs_real_3_t  gpts[7];

      const cs_quant_t  fq = cm->face[f];

      eval[f] = 0.; /* Reset value */

      /* Loop on face edges */

      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

        const short int  e = cm->f2e_ids[i];
        const short int v1 = cm->e2v_ids[2*e];
        const short int v2 = cm->e2v_ids[2*e+1];

        /* Evaluate the field at the three quadrature points */

        cs_quadrature_tria_7pts(fq.center, cm->xv + 3*v1, cm->xv + 3*v2,
                                cm->tef[i],
                                gpts, w);

        /* Evaluate the function for this time at the given coordinates */

        ac->func(time_eval, 7, NULL,
                 (const cs_real_t *)gpts, true,  /* compacted output ? */
                 ac->input,
                 _eval);

        cs_real_t  add = 0;
        for (int p = 0; p < 7; p++)
          add += w[p] * _eval[p];

        eval[f] += add;

      } /* Loop on face edges */
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of quadrature.", __func__);
    break;

  } /* Switch on the quadrature type */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the (scalar-valued) normal flux of
 *         a quantity defined by an analytic function.
 *         Use of a \ref cs_cell_mesh_t structure.
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in]      qtype      level of quadrature to use
 * \param[in, out] eval       array of evaluations at DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_flux_by_vector_analytic(const cs_cell_mesh_t      *cm,
                                        short int                  f,
                                        cs_real_t                  time_eval,
                                        void                      *context,
                                        cs_quadrature_type_t       qtype,
                                        cs_real_t                 *eval)
{
  assert(cs_flag_test(cm->flag, CS_FLAG_COMP_PFQ));
  assert(eval != NULL);

  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;

  const cs_quant_t  pfq = cm->face[f];

  eval[f] = 0.; /* Reset value */

  switch (qtype) {

  case CS_QUADRATURE_NONE:
  case CS_QUADRATURE_BARY:
    {
      cs_real_3_t  flux_xf = {0, 0, 0};

      /* Evaluate the function for this time at the face center */

      ac->func(time_eval, 1, NULL, pfq.center,
               true, /* dense output ? */
               ac->input,
               flux_xf);

      eval[f] = pfq.meas * _dp3(pfq.unitv, flux_xf);
    }
    break;

  case CS_QUADRATURE_BARY_SUBDIV:
    {
      assert(cs_flag_test(cm->flag,
                          CS_FLAG_COMP_EV |CS_FLAG_COMP_FE |CS_FLAG_COMP_FEQ));

      cs_real_3_t  _val, _xyz;

      /* Loop on face edges */

      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

        const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];
        const cs_real_t  *xv0 = cm->xv + 3*e2v[0];
        const cs_real_t  *xv1 = cm->xv + 3*e2v[1];

        for (int k = 0; k < 3; k++)
          _xyz[k] = cs_math_1ov3 * (pfq.center[k] + xv0[k] + xv1[k]);

        /* Evaluate the function for this time at the given coordinates */

        ac->func(time_eval, 1, NULL,
                 (const cs_real_t *)_xyz, true, /* dense output ? */
                 ac->input,
                 (cs_real_t *)_val);

        eval[f] += cm->tef[i] * _dp3(_val, pfq.unitv);
      }
    }
    break; /* BARY_SUBDIV */

  case CS_QUADRATURE_HIGHER:
    {
      assert(cs_flag_test(cm->flag,
                          CS_FLAG_COMP_EV |CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ));

      /* Four values by triangle --> 4 Gauss points
       * The flux returns by the analytic function is a vector. So the
       * size of _val is 12=4*3 */

      cs_real_t  w[4], _val[12];
      cs_real_3_t  gpts[4];

      /* Loop on face edges */

      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

        const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];
        const cs_real_t  *xv0 = cm->xv + 3*e2v[0];
        const cs_real_t  *xv1 = cm->xv + 3*e2v[1];

        /* Evaluate the field at the three quadrature points */

        cs_quadrature_tria_4pts(pfq.center, xv0, xv1, cm->tef[i], gpts, w);

        /* Evaluate the function for this time at the given coordinates */

        ac->func(time_eval, 4, NULL,
                 (const cs_real_t *)gpts, true,  /* dense output ? */
                 ac->input,
                 _val);

        cs_real_t  add = 0;
        for (int p = 0; p < 4; p++)
          add += w[p] * _dp3(_val+3*p, pfq.unitv);

        eval[f] += add;

      } /* Loop on face edges */
    }
    break;

  case CS_QUADRATURE_HIGHEST:
    {
      assert(cs_flag_test(cm->flag,
                          CS_FLAG_COMP_EV |CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ));

      /* Seven values by triangle --> 7 Gauss points
       * The flux returns by the analytic function is a vector. So the
       * size of _val is 21=7*3
       */

      cs_real_t  w[7], _val[21];
      cs_real_3_t  gpts[7];

      /* Loop on face edges */

      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

        const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];
        const cs_real_t  *xv0 = cm->xv + 3*e2v[0];
        const cs_real_t  *xv1 = cm->xv + 3*e2v[1];

        /* Evaluate the field at the three quadrature points */

        cs_quadrature_tria_7pts(pfq.center, xv0, xv1, cm->tef[i], gpts, w);

        /* Evaluate the function for this time at the given coordinates */

        ac->func(time_eval, 7, NULL,
                 (const cs_real_t *)gpts, true,  /* dense output ? */
                 ac->input,
                 _val);

        cs_real_t  add = 0;
        for (int p = 0; p < 7; p++)
          add += w[p] * _dp3(_val+3*p, pfq.unitv);

        eval[f] += add;

      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of quadrature.", __func__);
    break;

  }  /* Switch type of quadrature */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the normal vector-valued flux of
 *         a quantity defined by analytic function for the face f
 *         Use of a \ref cs_cell_mesh_t structure.
 *         Case of a vector-valued quantities.
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      f          local face id
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      context    pointer to a context structure
 * \param[in]      qtype      level of quadrature to use
 * \param[in, out] eval       array of evaluations at DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_vector_flux_by_analytic(const cs_cell_mesh_t      *cm,
                                        short int                  f,
                                        cs_real_t                  time_eval,
                                        void                      *context,
                                        cs_quadrature_type_t       qtype,
                                        cs_real_t                 *eval)
{
  assert(cs_flag_test(cm->flag, CS_FLAG_COMP_PFQ));
  assert(eval != NULL);

  cs_xdef_analytic_context_t  *cx = (cs_xdef_analytic_context_t *)context;

  const cs_quant_t  pfq = cm->face[f];

  cs_real_t  *f_eval = eval + 3*f;
  f_eval[0] = f_eval[1] = f_eval[2] = 0.;

  switch (qtype) {

  case CS_QUADRATURE_NONE:
  case CS_QUADRATURE_BARY:
    {
      cs_real_3_t  flux_xf = {0, 0, 0};

      /* Evaluate the function for this time at the cell center */

      cx->func(time_eval, 1, NULL, pfq.center, true,  /* dense output ? */
               cx->input,
               (cs_real_t *)flux_xf);

      for (int k = 0; k < 3; k++)
        f_eval[k] = pfq.meas * flux_xf[k];
    }
    break;

  case CS_QUADRATURE_BARY_SUBDIV:
    {
      assert(cs_flag_test(cm->flag,
                          CS_FLAG_COMP_EV| CS_FLAG_COMP_FE| CS_FLAG_COMP_FEQ));

      cs_real_3_t  _xyz, _eval;

      /* Loop on face edges */

      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

        const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];
        const cs_real_t  *xv0 = cm->xv + 3*e2v[0];
        const cs_real_t  *xv1 = cm->xv + 3*e2v[1];

        for (int k = 0; k < 3; k++)
          _xyz[k] = cs_math_1ov3 * (pfq.center[k] + xv0[k] + xv1[k]);

        /* Evaluate the function for this time at the given coordinates */

        cx->func(time_eval, 1, NULL,
                 (const cs_real_t *)_xyz, true,  /* dense output ? */
                 cx->input,
                 (cs_real_t *)_eval);

        for (int k = 0; k < 3; k++)
          f_eval[k] += cm->tef[i] * _eval[k];
      }
    }
    break; /* BARY_SUBDIV */

  case CS_QUADRATURE_HIGHER:
    {
      assert(cs_flag_test(cm->flag,
                          CS_FLAG_COMP_EV| CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ));

      /* Four values by triangle --> 4 Gauss points
       * The flux returns by the analytic function is a 3x3 tensor. */

      cs_real_t  w[4], _eval[12];
      cs_real_3_t  gpts[4];

      /* Loop on face edges */

      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

        const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];
        const cs_real_t  *xv0 = cm->xv + 3*e2v[0];
        const cs_real_t  *xv1 = cm->xv + 3*e2v[1];

        /* Evaluate the field at the three quadrature points */

        cs_quadrature_tria_4pts(pfq.center, xv0, xv1, cm->tef[i], gpts, w);

        /* Evaluate the function for this time at the given coordinates */

        cx->func(time_eval, 4, NULL,
                 (const cs_real_t *)gpts, true,  /* dense output ? */
                 cx->input,
                 _eval);

        for (int p = 0; p < 4; p++)
          for (int k = 0; k < 3; k++)
            f_eval[k] += w[p] * _eval[3*p+k];
      }
    }
    break;

  case CS_QUADRATURE_HIGHEST:
    {
      assert(cs_flag_test(cm->flag,
                          CS_FLAG_COMP_EV| CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ));

      /* Seven values by triangle --> 7 Gauss points
       * The flux returns by the analytic function is a 3x3 tensor. */

      cs_real_t  w[7], _eval[21];
      cs_real_3_t  gpts[7];

      /* Loop on face edges */

      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

        const short int  *e2v = cm->e2v_ids + 2*cm->f2e_ids[i];
        const cs_real_t  *xv0 = cm->xv + 3*e2v[0];
        const cs_real_t  *xv1 = cm->xv + 3*e2v[1];

        /* Evaluate the field at the three quadrature points */

        cs_quadrature_tria_7pts(pfq.center, xv0, xv1, cm->tef[i], gpts, w);

        /* Evaluate the function for this time at the given coordinates */

        cx->func(time_eval, 7, NULL,
                 (const cs_real_t *)gpts, true,  /* dense output ? */
                 cx->input,
                 _eval);

        for (int p = 0; p < 7; p++)
          for (int k = 0; k < 3; k++)
            f_eval[k] += w[p] * _eval[3*p+k];
      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of quadrature.", __func__);
    break;

  }  /* Switch type of quadrature */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the reduction by averages of a
 *         analytic function by a cellwise process (usage of a
 *         \ref cs_cell_mesh_t structure).
 *         This evaluation hinges on the computation of integrals (faces first,
 *         then cell)
 *         Scalar-valued case.
 *
 * \param[in]      cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      qtype    quadrature type
 * \param[in]      context  pointer to a context structure
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_scal_avg_reduction_by_analytic(const cs_cell_mesh_t    *cm,
                                               cs_real_t                t_eval,
                                               void                    *context,
                                               cs_quadrature_type_t     qtype,
                                               cs_real_t               *eval)
{
  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(context != NULL);
  assert(cs_flag_test(cm->flag,
                      CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                      CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  const int dim = 1;
  const short int nf = cm->n_fc;

  cs_quadrature_tetra_integral_t
    *q_tet = cs_quadrature_get_tetra_integral(dim, qtype);
  cs_quadrature_tria_integral_t
    *q_tri = cs_quadrature_get_tria_integral(dim, qtype);
  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;
  cs_real_t *c_eval = eval + nf;

  cs_xdef_cw_eval_fc_int_by_analytic(cm, t_eval,
                                     ac->func, ac->input,
                                     dim,
                                     q_tet, q_tri,
                                     c_eval, eval);

  /* Compute the averages */

  for (short int f = 0; f < nf; f++)
    eval[f] /= cm->face[f].meas;
  eval[nf] /= cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer for evaluating the reduction by averages of a
 *         analytic function by a cellwise process (usage of a
 *         \ref cs_cell_mesh_t structure).
 *         This evaluation hinges on the computation of integrals (faces first,
 *         then cell)
 *         Vector-valued case.
 *
 * \param[in]      cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      qtype    quadrature type
 * \param[in]      context  pointer to a context structure
 * \param[in, out] eval     result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_cw_eval_vect_avg_reduction_by_analytic(const cs_cell_mesh_t    *cm,
                                               cs_real_t                t_eval,
                                               void                    *context,
                                               cs_quadrature_type_t     qtype,
                                               cs_real_t               *eval)
{
  if (eval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(context != NULL);
  assert(cs_flag_test(cm->flag,
                      CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                      CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  const int  dim = 3;
  const short int nf = cm->n_fc;

  cs_quadrature_tetra_integral_t
    *q_tet = cs_quadrature_get_tetra_integral(3, qtype);
  cs_quadrature_tria_integral_t
    *q_tri = cs_quadrature_get_tria_integral(3, qtype);
  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)context;
  cs_real_t *c_eval = eval + dim*nf;

  cs_xdef_cw_eval_fc_int_by_analytic(cm, t_eval,
                                     ac->func, ac->input,
                                     dim,
                                     q_tet, q_tri,
                                     c_eval, eval);

  /* Compute the averages */

  for (short int f = 0; f < nf; f++) {
    const cs_real_t _os = 1. / cm->face[f].meas;
    cs_real_t *f_eval = eval + dim*f;
    f_eval[0] *= _os, f_eval[1] *= _os, f_eval[2] *= _os;
  }

  const cs_real_t _ov = 1. / cm->vol_c;
  c_eval[0] *= _ov, c_eval[1] *= _ov, c_eval[2] *= _ov;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
