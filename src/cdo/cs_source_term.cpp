/*============================================================================
 * Functions and structures to deal with source term computations
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_mem.h"

#include "cdo/cs_basis_func.h"
#include "cdo/cs_evaluate.h"
#include "cdo/cs_hho_builder.h"
#include "cdo/cs_hodge.h"
#include "base/cs_log.h"
#include "cdo/cs_macfb_builder.h"
#include "base/cs_math.h"
#include "cdo/cs_scheme_geometry.h"
#include "base/cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cdo/cs_source_term.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_SOURCE_TERM_DBG 0

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_st[] =
  " Stop setting an empty cs_xdef_t structure.\n"
  " Please check your settings.\n";

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the mask associated to each cell from the mask related to
 *         the given source term structure
 *
 * \param[in]      st          pointer to a cs_xdef_t structure
 * \param[in]      st_id       id related to this source term
 * \param[in, out] cell_mask   mask related to each cell to be updated
 */
/*----------------------------------------------------------------------------*/

static void
_set_mask(const cs_xdef_t     *st,
          int                  st_id,
          cs_mask_t           *cell_mask)
{
  if (st == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  /* value of the mask for the source term */

  const cs_mask_t  mask = (1 << st_id);

  if (st->meta & CS_FLAG_FULL_LOC) /* All cells are selected */
#   pragma omp parallel for if (cs_cdo_quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < cs_cdo_quant->n_cells; i++) cell_mask[i] |= mask;

  else {

    /* Retrieve information from the volume zone structure */

    const cs_zone_t *z = cs_volume_zone_by_id(st->z_id);
    for (cs_lnum_t i = 0; i < z->n_elts; i++)
      cell_mask[z->elt_ids[i]] |= mask;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the reduction onto the cell polynomial space of a function
 *         defined by a constant value
 *
 * \param[in]       const_val  constant value
 * \param[in]       cbf        pointer to a structure for face basis functions
 * \param[in]       xv1        first vertex
 * \param[in]       xv2        second vertex
 * \param[in]       xv3        third vertex
 * \param[in]       xv4        fourth vertex
 * \param[in]       vol        volume of the tetrahedron
 * \param[in, out]  cb         pointer to a cs_cell_builder_structure_t
 * \param[in, out]  array      array storing values to compute
 */
/*----------------------------------------------------------------------------*/

static void
_hho_add_tetra_by_val(cs_real_t                        const_val,
                      const cs_basis_func_t           *cbf,
                      const cs_real_3_t                xv1,
                      const cs_real_3_t                xv2,
                      const cs_real_3_t                xv3,
                      const cs_real_3_t                xv4,
                      const double                     vol,
                      cs_cell_builder_t               *cb,
                      cs_real_t                        array[])
{
  cs_real_3_t  *gpts = cb->vectors;
  cs_real_t  *gw = cb->values;
  cs_real_t  *phi_eval = cb->values + 15;

  /* Compute Gauss points and related weights */

  cs_quadrature_tet_15pts(xv1, xv2, xv3, xv4, vol, gpts, gw);

  for (short int gp = 0; gp < 15; gp++) {

    cbf->eval_all_at_point(cbf, gpts[gp], phi_eval);

    const cs_real_t  w = gw[gp] * const_val;
    for (short int i = 0; i < cbf->size; i++)
      array[i] += w * phi_eval[i];

  }  /* End of loop on Gauss points */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the reduction onto the cell polynomial space of a function
 *         defined by an analytical expression depending on the location and
 *         the current time
 *
 * \param[in]       ac        pointer to an analytical definition
 * \param[in]       cbf       pointer to a structure for face basis functions
 * \param[in]       xv1       first vertex
 * \param[in]       xv2       second vertex
 * \param[in]       xv3       third vertex
 * \param[in]       xv4       fourth vertex
 * \param[in]       vol       volume of the tetrahedron
 * \param[in]       time_eval physical time at which one evaluates the term
 * \param[in, out]  cb        pointer to a cs_cell_builder_structure_t
 * \param[in, out]  array     array storing values to compute
 */
/*----------------------------------------------------------------------------*/

static void
_hho_add_tetra_by_ana(const cs_xdef_analytic_context_t   *ac,
                      const cs_basis_func_t              *cbf,
                      const cs_real_3_t                   xv1,
                      const cs_real_3_t                   xv2,
                      const cs_real_3_t                   xv3,
                      const cs_real_3_t                   xv4,
                      const double                        vol,
                      cs_real_t                           time_eval,
                      cs_cell_builder_t                  *cb,
                      cs_real_t                           array[])
{
  cs_real_3_t  *gpts = cb->vectors;
  cs_real_t  *gw = cb->values;
  cs_real_t  *ana_eval = cb->values + 15;
  cs_real_t  *phi_eval = cb->values + 30;

  /* Compute Gauss points and related weights */

  cs_quadrature_tet_15pts(xv1, xv2, xv3, xv4, vol, gpts, gw);

  /* Evaluate the analytical function at the Gauss points */

  ac->func(time_eval, 15, nullptr, (const cs_real_t *)gpts, true,
           ac->input, ana_eval);

  for (short int gp = 0; gp < 15; gp++) {

    cbf->eval_all_at_point(cbf, gpts[gp], phi_eval);

    const cs_real_t  w = gw[gp] * ana_eval[gp];
    for (short int i = 0; i < cbf->size; i++)
      array[i] += w * phi_eval[i];

  }  /* End of loop on Gauss points */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the reduction onto the cell polynomial space of a function
 *         defined by an analytical expression depending on the location and
 *         the current time. This is the vector case.
 *
 * \param[in]       ac         pointer to an analytical definition
 * \param[in]       cbf        pointer to a structure for cell basis functions
 * \param[in]       xv1        first vertex
 * \param[in]       xv2        second vertex
 * \param[in]       xv3        third vertex
 * \param[in]       xv4        third vertex
 * \param[in]       vol        volume of the tetrahedron
 * \param[in]       time_eval  physical time at which one evaluates the term
 * \param[in, out]  cb         pointer to a cs_cell_builder_structure_t
 * \param[in, out]  array      array storing values to compute
 */
/*----------------------------------------------------------------------------*/

static void
_hho_add_tetra_by_ana_vd(const cs_xdef_analytic_context_t  *ac,
                         const cs_basis_func_t             *cbf,
                         const cs_real_3_t                  xv1,
                         const cs_real_3_t                  xv2,
                         const cs_real_3_t                  xv3,
                         const cs_real_3_t                  xv4,
                         const double                       vol,
                         cs_real_t                          time_eval,
                         cs_cell_builder_t                 *cb,
                         cs_real_t                          array[])
{
  cs_real_3_t  *gpts = cb->vectors;

  /* cb->values is big enough to store:
     gw+ana_eval+phi_eval= 15+3*15+cell_basis
     = 64  < 78  + 12      for k=1
     = 70  < 465 + 30      for k=2
  */

  cs_real_t  *gw = cb->values;
  cs_real_t  *ana_eval = cb->values + 15;
  cs_real_t  *phi_eval = cb->values + 15 + 3*15;

  /* Compute Gauss points and related weights */

  cs_quadrature_tet_15pts(xv1, xv2, xv3, xv4, vol, gpts, gw);

  /* Evaluate the analytical function at the Gauss points */

  ac->func(time_eval, 15, nullptr, (const cs_real_t *)gpts, true,
           ac->input, ana_eval);

  for (short int gp = 0; gp < 15; gp++) {

    cbf->eval_all_at_point(cbf, gpts[gp], phi_eval);

    for (short int i = 0; i < cbf->size; i++) {

      const double  gcoef = gw[gp] * phi_eval[i];

      /* x-component */
      array[i              ] += gcoef * ana_eval[3*gp];
      /* y-component */
      array[i +   cbf->size] += gcoef * ana_eval[3*gp+1];
      /* z-component */
      array[i + 2*cbf->size] += gcoef * ana_eval[3*gp+2];

    }

  }  /* End of loop on Gauss points */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a function to compute a source term and update flag storing
 *        which cellwise quantities have to be defined for this kind of
 *        computation.
 *        Case of CDO vertex-based schemes
 *
 * \param[in]      st_def    pointer to the definition
 * \param[in, out] sys_flag  metadata about the algebraic system
 * \param[in, out] msh_flag  metadata about cellwise quantities
 *
 * \return a pointer to the function used for the source term computation
 */
/*----------------------------------------------------------------------------*/

static cs_source_term_cellwise_t *
_set_vb_function(const cs_xdef_t *st_def,
                 cs_flag_t       *sys_flag,
                 cs_eflag_t      *msh_flag)
{
  assert(st_def != nullptr);

  cs_source_term_cellwise_t *func = nullptr;

  switch (st_def->type) {

  case CS_XDEF_BY_VALUE:
    /* ---------------- */

    if (st_def->meta & CS_FLAG_DUAL) {

      *msh_flag |= CS_FLAG_COMP_PVQ;
      if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
        func = cs_source_term_dcvd_by_value;
      else
        func = cs_source_term_dcsd_by_value;

    }
    else {

      assert(st_def->meta & CS_FLAG_PRIMAL);
      *msh_flag |= CS_FLAG_COMP_PV;
      func = cs_source_term_pvsp_by_value;

    }
    break; /* definition by value */

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    /* ---------------------------- */

    if (st_def->meta & CS_FLAG_DUAL) {

      switch (st_def->qtype) {
      case CS_QUADRATURE_NONE:
        *msh_flag |= CS_FLAG_COMP_PVQ;
        func = cs_source_term_dcsd_none_by_analytic;
        break;

      case CS_QUADRATURE_BARY:
        *msh_flag |=
          CS_FLAG_COMP_PVQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PFC |
          CS_FLAG_COMP_FE  | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV  |
          CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PEC | CS_FLAG_COMP_DEQ;
        func = cs_source_term_dcsd_bary_by_analytic;
        break;

      case CS_QUADRATURE_BARY_SUBDIV:
        *msh_flag |=
          CS_FLAG_COMP_EV | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PFC |
          CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_DEQ;
        func = cs_source_term_dcsd_q1o1_by_analytic;
        break;

      case CS_QUADRATURE_HIGHER:
        *msh_flag |=
          CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PFC | CS_FLAG_COMP_FE  |
          CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV  | CS_FLAG_COMP_PVQ |
          CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DEQ;
        func = cs_source_term_dcsd_q10o2_by_analytic;
        break;

      case CS_QUADRATURE_HIGHEST:
        *msh_flag |=
          CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE  |
          CS_FLAG_COMP_EV  | CS_FLAG_COMP_PFC | CS_FLAG_COMP_FEQ |
          CS_FLAG_COMP_DEQ;
        func = cs_source_term_dcsd_q5o3_by_analytic;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid type of quadrature for computing a source term"
                  " with CDOVB schemes");
      } /* quad_type */

    }
    else {

      assert(st_def->meta & CS_FLAG_PRIMAL);
      *msh_flag |= CS_FLAG_COMP_PV;
      func = cs_source_term_pvsp_by_analytic;

    }
    break; /* definition by analytic */

  case CS_XDEF_BY_ARRAY:
    /* ---------------- */
    {
      cs_xdef_array_context_t *ctx
        = static_cast<cs_xdef_array_context_t *>(st_def->context);

      *msh_flag |= CS_FLAG_COMP_PVQ;

      if (cs_flag_test(ctx->value_location, cs_flag_primal_vtx)) {

        if (st_def->meta & CS_FLAG_DUAL) {

          if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
            func = cs_source_term_dcvd_by_pv_array;
          else
            func = cs_source_term_dcsd_by_pv_array;

        }
        else {
          assert(st_def->meta & CS_FLAG_PRIMAL);
          func = cs_source_term_pvsp_by_array;
        }

      }
      else if (cs_flag_test(ctx->value_location, cs_flag_dual_cell_byc))  {

        /* Should be before cs_flag_dual_cell */
        if (st_def->meta & CS_FLAG_DUAL)
          func = cs_source_term_dcsd_by_c2v_array;
        else {
          assert(st_def->meta & CS_FLAG_PRIMAL);
          func = cs_source_term_pvsp_by_c2v_array;
        }

      }
      else if (cs_flag_test(ctx->value_location, cs_flag_dual_cell)) {

        if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
          func = cs_source_term_dcvd_by_pv_array;
        else
          func = cs_source_term_dcsd_by_pv_array;

      }
      else if (cs_flag_test(ctx->value_location, cs_flag_primal_cell))
        func = cs_source_term_dcsd_by_pc_array;
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of definition for a source term in CDOVB",
                  __func__);
    }
    break;

  case CS_XDEF_BY_DOF_FUNCTION:
    /* ----------------------- */

    assert(st_def->meta & CS_FLAG_DUAL);
    *msh_flag |= CS_FLAG_COMP_PVQ;
    func = cs_source_term_dcsd_by_dof_func;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of definition for a source term in CDOVB",
              __func__);
    break;

  } /* switch one the type of definition */

  return func;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_init_sharing(const cs_cdo_quantities_t    *quant,
                            const cs_cdo_connect_t       *connect)
{
  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default flag related to a source term according to the
 *         numerical scheme chosen for discretizing an equation
 *
 * \param[in]   scheme      numerical scheme used for the discretization
 * \param[out]  state_flag  flag describing the status of the source term
 * \param[out]  meta_flag   additional flags associated to a source term
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_default_flag(cs_param_space_scheme_t    scheme,
                                cs_flag_t                 *state_flag,
                                cs_flag_t                 *meta_flag)
{
  switch (scheme) {
  case CS_SPACE_SCHEME_CDOVB:
    *state_flag = CS_FLAG_STATE_DENSITY;
    *meta_flag = cs_flag_dual_cell; /* Predefined mask */
    break;

  case CS_SPACE_SCHEME_CDOEB:
    *state_flag = CS_FLAG_STATE_FLUX;
    *meta_flag = cs_flag_dual_face; /* Predefined mask */
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_CDOCB:
  case CS_SPACE_SCHEME_MACFB:
    *state_flag = CS_FLAG_STATE_DENSITY;
    *meta_flag = cs_flag_primal_cell; /* Predefined mask */
    break;

  case CS_SPACE_SCHEME_CDOVCB:
  case CS_SPACE_SCHEME_HHO_P0:
  case CS_SPACE_SCHEME_HHO_P1:
  case CS_SPACE_SCHEME_HHO_P2:
    *state_flag = CS_FLAG_STATE_DENSITY;
    *meta_flag = CS_FLAG_PRIMAL;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid numerical scheme to set a source term.",
              __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set advanced parameters which are defined by default in a
 *         source term structure.
 *
 * \param[in, out]  st        pointer to a cs_xdef_t structure
 * \param[in]       flag      CS_FLAG_DUAL or CS_FLAG_PRIMAL
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_reduction(cs_xdef_t     *st,
                             cs_flag_t      flag)
{
  if (st == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  if (st->meta & flag)
    return; /* Nothing to do */

  cs_flag_t  save_meta = st->meta;

  st->meta = 0;

  /* Set unchanged parts of the existing flag */

  if (save_meta & CS_FLAG_SCALAR) st->meta |= CS_FLAG_SCALAR;
  if (save_meta & CS_FLAG_VECTOR) st->meta |= CS_FLAG_VECTOR;
  if (save_meta & CS_FLAG_TENSOR) st->meta |= CS_FLAG_TENSOR;
  if (save_meta & CS_FLAG_BORDER) st->meta |= CS_FLAG_BORDER;
  if (save_meta & CS_FLAG_BY_CELL) st->meta |= CS_FLAG_BY_CELL;
  if (save_meta & CS_FLAG_FULL_LOC) st->meta |= CS_FLAG_FULL_LOC;

  if (flag & CS_FLAG_DUAL) {
    assert(save_meta & CS_FLAG_PRIMAL);
    if (save_meta & CS_FLAG_VERTEX)
      st->meta |= CS_FLAG_DUAL | CS_FLAG_CELL;
    else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Stop modifying the source term flag.\n"
                " This case is not handled.", __func__);
  }
  else if (flag & CS_FLAG_PRIMAL) {
    assert(save_meta & CS_FLAG_DUAL);
    if (save_meta & CS_FLAG_CELL)
      st->meta |= CS_FLAG_PRIMAL | CS_FLAG_VERTEX;
    else
      bft_error(__FILE__, __LINE__, 0,
                " Stop modifying the source term flag.\n"
                " This case is not handled.");
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " Stop modifying the source term flag.\n"
              " This case is not handled.");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get metadata related to the given source term structure
 *
 * \param[in, out]  st          pointer to a cs_xdef_t structure
 *
 * \return the value of the flag related to this source term
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_source_term_get_flag(const cs_xdef_t  *st)
{
  if (st == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  return st->meta;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize data to build the source terms
 *
 * \param[in]      space_scheme    scheme used to discretize in space
 * \param[in]      n_source_terms  number of source terms
 * \param[in]      source_terms    pointer to the definitions of source terms
 * \param[in, out] compute_source  array of function pointers
 * \param[in, out] sys_flag        metadata about the algebraic system
 * \param[in, out] source_mask     pointer to an array storing in a compact way
 *                                 which source term is defined in a given cell
 *
 * \return a flag which indicates what to build in a cell mesh structure
 */
/*----------------------------------------------------------------------------*/

cs_eflag_t
cs_source_term_init(cs_param_space_scheme_t       space_scheme,
                    const int                     n_source_terms,
                    cs_xdef_t             *const *source_terms,
                    cs_source_term_cellwise_t    *compute_source[],
                    cs_flag_t                    *sys_flag,
                    cs_mask_t                    *source_mask[])
{
  if (n_source_terms > CS_N_MAX_SOURCE_TERMS)
    bft_error(__FILE__, __LINE__, 0,
              " Limitation to %d source terms has been reached!",
              CS_N_MAX_SOURCE_TERMS);

  cs_eflag_t  msh_flag = 0;
  *source_mask         = nullptr;
  for (short int i = 0; i < CS_N_MAX_SOURCE_TERMS; i++)
    compute_source[i] = nullptr;

  if (n_source_terms == 0)
    return msh_flag;

  bool  need_mask = false;

  for (int st_id = 0; st_id < n_source_terms; st_id++) {

    const cs_xdef_t  *st_def = source_terms[st_id];

    if (st_def->meta & CS_FLAG_PRIMAL) {
      if (space_scheme == CS_SPACE_SCHEME_CDOVB ||
          space_scheme == CS_SPACE_SCHEME_CDOVCB) {
        msh_flag |= CS_FLAG_COMP_PVQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ |
          CS_FLAG_COMP_EV  | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ;
        *sys_flag |= CS_FLAG_SYS_MASS_MATRIX | CS_FLAG_SYS_SOURCES_HLOC;
      }
      else if (space_scheme == CS_SPACE_SCHEME_CDOEB) {
        msh_flag |= CS_FLAG_COMP_DFQ | CS_FLAG_COMP_EV;
        *sys_flag |= CS_FLAG_SYS_MASS_MATRIX | CS_FLAG_SYS_SOURCES_HLOC;
      }
    }

    /* Not defined on the whole mesh */

    if ((st_def->meta & CS_FLAG_FULL_LOC) == 0)
      need_mask = true;

    /* Assign the right function according to the space discretization */

    switch (space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      /* --------------------- */
      msh_flag |= CS_FLAG_COMP_PV;
      compute_source[st_id] = _set_vb_function(st_def, sys_flag, &msh_flag);
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      /* ---------------------- */
      if (st_def->meta & CS_FLAG_DUAL) {

        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of definition for a source term in CDOVCB",
                  __func__);

        /* TODO:
           case CS_XDEF_BY_VALUE:
           cs_source_term_vcsd_by_value; --> case CS_QUADRATURE_BARY:

           case CS_XDEF_BY_ANALYTIC_FUNCTION:
           cs_source_term_vcsd_q1o1_by_analytic; --> case CS_QUADRATURE_BARY:
           cs_source_term_vcsd_q10o2_by_analytic; --> case CS_QUADRATURE_HIGHER:
           cs_source_term_vcsd_q5o3_by_analytic; --> case CS_QUADRATURE_HIGHEST:
        */

      }
      else {
        assert(st_def->meta & CS_FLAG_PRIMAL);

        switch (st_def->type) {

        case CS_XDEF_BY_VALUE:
          msh_flag |= CS_FLAG_COMP_PV;
          compute_source[st_id] = cs_source_term_vcsp_by_value;
          break;

        case CS_XDEF_BY_ANALYTIC_FUNCTION:
          msh_flag |= CS_FLAG_COMP_PV;
          compute_source[st_id] = cs_source_term_vcsp_by_analytic;
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    " Invalid type of definition for a source term in CDOVCB");
          break;

        } /* switch def_type */

      }
      break; /* CDOVCB */

    case CS_SPACE_SCHEME_CDOEB:
      /* --------------------- */
      switch (st_def->type) {

      case CS_XDEF_BY_VALUE:
        compute_source[st_id] = cs_source_term_dfsf_by_value;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of definition for a source term in CDOEB",
                  __func__);
        break;
      }
      break; /* CDOEB */

    case CS_SPACE_SCHEME_CDOCB:
      /* --------------------- */

      if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
        bft_error(__FILE__, __LINE__, 0, "%s: Invalid case", __func__);

      switch (st_def->type) {

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        msh_flag |= CS_FLAG_COMP_PV;

        /* Switch only to allow more precise quadratures */

        switch (st_def->qtype) {

        case CS_QUADRATURE_HIGHEST:
        case CS_QUADRATURE_HIGHER:
        case CS_QUADRATURE_BARY_SUBDIV:
          msh_flag |= CS_FLAG_COMP_EV | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FEQ |
            CS_FLAG_COMP_HFQ| CS_FLAG_COMP_PEQ | CS_FLAG_COMP_FE;
          compute_source[st_id] = cs_source_term_fcb_pcsd_by_analytic;
          break;

        default:
          /* CS_QUADRATURE_BARY or CS_QUADRATURE_NONE */
          compute_source[st_id] = cs_source_term_fcb_pcsd_bary_by_analytic;
          break;

        } /* Switch on the type of quadrature */
        break;

      case CS_XDEF_BY_ARRAY:
        compute_source[st_id] = cs_source_term_fcb_pcsd_by_array;
        break;

      case CS_XDEF_BY_DOF_FUNCTION:
        compute_source[st_id] = cs_source_term_fcb_pcsd_by_dof_func;
        break;

      case CS_XDEF_BY_VALUE:
        compute_source[st_id] = cs_source_term_fcb_pcsd_by_value;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of source term definition in CDO-CB",
                  __func__);
        break;

      }
      break;

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
      /* ---------------------- */
      switch (st_def->type) {

      case CS_XDEF_BY_VALUE:
        if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
          compute_source[st_id] = cs_source_term_fb_pcvd_by_value;
        else
          compute_source[st_id] = cs_source_term_fcb_pcsd_by_value;
        break;

      case CS_XDEF_BY_DOF_FUNCTION:
        if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
          compute_source[st_id] = cs_source_term_fb_pcvd_by_dof_func;
        else
          compute_source[st_id] = cs_source_term_fcb_pcsd_by_dof_func;
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        msh_flag |= CS_FLAG_COMP_PV;
        if ((*sys_flag) & CS_FLAG_SYS_VECTOR) {

          /* Switch only to allow more precise quadratures */

          if (st_def->qtype == CS_QUADRATURE_BARY)
            compute_source[st_id] = cs_source_term_fb_pcvd_bary_by_analytic;

          else {

            /* TODO: Are all these flags really necessary? Check in the */
            /* integration */

            msh_flag |= CS_FLAG_COMP_EV  |CS_FLAG_COMP_EF | CS_FLAG_COMP_PFQ |
                        CS_FLAG_COMP_FEQ |CS_FLAG_COMP_HFQ| CS_FLAG_COMP_PEQ |
                        CS_FLAG_COMP_FE;
            compute_source[st_id] = cs_source_term_fb_pcvd_by_analytic;
          }
        }
        else { /* Scalar-valued case */

          /* Switch only to allow more precise quadratures */

          switch (st_def->qtype) {

          case CS_QUADRATURE_HIGHEST:
          case CS_QUADRATURE_HIGHER:
          case CS_QUADRATURE_BARY_SUBDIV:
            /* TODO: Are all these flags really necessary? Check in the */
            /* integration */

            msh_flag |= CS_FLAG_COMP_EV | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FEQ |
              CS_FLAG_COMP_HFQ| CS_FLAG_COMP_PEQ | CS_FLAG_COMP_FE;
            compute_source[st_id] = cs_source_term_fcb_pcsd_by_analytic;
            break;

          default:
            /* CS_QUADRATURE_BARY, CS_QUADRATURE_NONE */
            compute_source[st_id] = cs_source_term_fcb_pcsd_bary_by_analytic;
            break;

          } /* Switch on the type of quadrature */
        }
        break;

      case CS_XDEF_BY_ARRAY:
        if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
          compute_source[st_id] = cs_source_term_fb_pcvd_by_array;
        else
          compute_source[st_id] = cs_source_term_fcb_pcsd_by_array;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of definition for a source term in CDOFB",
                  __func__);
        break;

      } /* switch def_type */
      break;

    case CS_SPACE_SCHEME_HHO_P1:
    case CS_SPACE_SCHEME_HHO_P2:
      /* ---------------------- */
      switch (st_def->type) {

      case CS_XDEF_BY_VALUE:
        if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: Invalid type of definition for a source term in HHO",
                    __func__);
        else
          compute_source[st_id] = cs_source_term_hhosd_by_value;
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
          compute_source[st_id] = cs_source_term_hhovd_by_analytic;
        else
          compute_source[st_id] = cs_source_term_hhosd_by_analytic;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid type of definition for a source term in HHO",
                  __func__);
        break;

      } /* switch def_type */
      break;

    case CS_SPACE_SCHEME_MACFB:
      /* ---------------------- */
      switch (st_def->type) {

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        msh_flag |= CS_FLAG_COMP_PV;
        if ((*sys_flag) & CS_FLAG_SYS_VECTOR) {

          /* TODO: Are all these flags really necessary? Check in the */
          /* integration */

          msh_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_HFQ;
          compute_source[st_id] = cs_source_term_macfb_pcvd_by_analytic;
        }
        else { /* Scalar-valued case */

          bft_error(__FILE__, __LINE__, 0,
                    "%s: Invalid type of definition for a scalar source term in MACFB",
                    __func__);
        }
        break;

      case CS_XDEF_BY_VALUE:
      case CS_XDEF_BY_DOF_FUNCTION:
      case CS_XDEF_BY_ARRAY:
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid type of definition for a source term in MACFB",
                  __func__);
        break;

      } /* switch def_type */
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid space scheme for setting the source term.",
                __func__);
      break;

    } /* Switch on space scheme */

  } /* Loop on source terms */

  if (need_mask) {

    const cs_lnum_t n_cells = cs_cdo_quant->n_cells;

    /* Initialize mask buffer */

    cs_mask_t *mask = nullptr;
    CS_MALLOC(mask, n_cells, cs_mask_t);
#pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (int i = 0; i < n_cells; i++)
      mask[i] = 0;

    for (int st_id = 0; st_id < n_source_terms; st_id++)
      _set_mask(source_terms[st_id], st_id, mask);

    *source_mask = mask;

  } /* Build a tag related to the source terms defined in each cell */

  return msh_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the local contributions of source terms in a cell
 *
 * \param[in]      n_source_terms  number of source terms
 * \param[in]      source_terms    pointer to the definitions of source terms
 * \param[in]      cm              pointer to a cs_cell_mesh_t structure
 * \param[in]      source_mask     array storing in a compact way which source
 *                                 term is defined in a given cell
 * \param[in]      compute_source  array of function pointers
 * \param[in]      time_eval       physical time at which one evaluates the term
 * \param[in, out] input           pointer to an element cast on-the-fly
 * \param[in, out] cb              pointer to a cs_cell_builder_t structure
 * \param[in, out] result          array storing the result of the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute_cellwise(const int                  n_source_terms,
                                cs_xdef_t *const          *source_terms,
                                const cs_cell_mesh_t      *cm,
                                const cs_mask_t           *source_mask,
                                cs_source_term_cellwise_t *compute_source[],
                                cs_real_t                  time_eval,
                                void                      *input,
                                cs_cell_builder_t         *cb,
                                cs_real_t                 *result)
{
  if (n_source_terms < 1)
    return;

  if (source_mask == nullptr) { /* Source terms are defined on the whole mesh */

    for (short int st_id = 0; st_id < n_source_terms; st_id++) {

      cs_source_term_cellwise_t *compute = compute_source[st_id];

      /* result is updated inside */

      compute(source_terms[st_id], cm, time_eval, cb, input, result);

    } /* Loop on source terms */
  }
  else { /* Some source terms are only defined on a selection of cells */

    for (short int st_id = 0; st_id < n_source_terms; st_id++) {

      const cs_mask_t st_mask = (1 << st_id);
      if (source_mask[cm->c_id] & st_mask) {

        cs_source_term_cellwise_t *compute = compute_source[st_id];

        /* result is updated inside */

        compute(source_terms[st_id], cm, time_eval, cb, input, result);

      } /* Compute the source term on this cell */

    } /* Loop on source terms */

  } /* Source terms are defined on the whole domain or not ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar potential defined at primal vertices by a constant
 *         value.
 *         A discrete Hodge operator has to be computed before this call and
 *         given as input parameter
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_pvsp_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_real_t                  time_eval,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values)
{
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  cs_hodge_t  *mass_hodge = (cs_hodge_t *)input;

  assert(values != nullptr && cm != nullptr && cb != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PV));
  assert(mass_hodge != nullptr);
  assert(mass_hodge->matrix != nullptr);

  const cs_real_t *s_values = (const cs_real_t *)source->context;
  const cs_real_t  pot_value = s_values[0];

  /* Retrieve the values of the potential at each cell vertices */

  double  *eval = cb->values;
  for (short int v = 0; v < cm->n_vc; v++)
    eval[v] = pot_value;

  /* Multiply these values by a cellwise Hodge operator previously computed */

  double  *hdg_eval = cb->values + cm->n_vc;
  cs_sdm_square_matvec(mass_hodge->matrix, eval, hdg_eval);

  for (short int v = 0; v < cm->n_vc; v++)
    values[v] += hdg_eval[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar potential defined at primal vertices by an
 *         analytical function.
 *         A discrete Hodge operator has to be computed before this call and
 *         given as input parameter
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_pvsp_by_analytic(const cs_xdef_t           *source,
                                const cs_cell_mesh_t      *cm,
                                cs_real_t                  time_eval,
                                cs_cell_builder_t         *cb,
                                void                      *input,
                                double                    *values)
{
  if (source == nullptr)
    return;

  cs_hodge_t  *mass_hodge = (cs_hodge_t *)input;

  assert(values != nullptr && cm != nullptr && cb != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PV));
  assert(mass_hodge != nullptr);
  assert(mass_hodge->matrix != nullptr);

  cs_xdef_analytic_context_t  *ac =
    (cs_xdef_analytic_context_t *)source->context;

  /* Retrieve the values of the potential at each cell vertices */

  double  *eval = cb->values;
  ac->func(time_eval, cm->n_vc, nullptr, cm->xv,
           true, /* compacted output ? */
           ac->input,
           eval);

  /* Multiply these values by a cellwise Hodge operator previously computed */

  double  *hdg_eval = cb->values + cm->n_vc;
  cs_sdm_square_matvec(mass_hodge->matrix, eval, hdg_eval);

  for (short int v = 0; v < cm->n_vc; v++)
    values[v] += hdg_eval[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the contribution for a cell related to a source term and
 *        add it to the given array of values.
 *        Case of a scalar potential defined at primal vertices by an array.
 *        A discrete Hodge operator has to be computed before this call and
 *        given as an input parameter
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_pvsp_by_array(const cs_xdef_t      *source,
                             const cs_cell_mesh_t *cm,
                             cs_real_t             time_eval,
                             cs_cell_builder_t    *cb,
                             void                 *input,
                             double               *values)
{
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  cs_hodge_t *mass_hodge = (cs_hodge_t *)input;

  assert(values != nullptr && cm != nullptr && cb != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PV));
  assert(mass_hodge != nullptr);
  assert(mass_hodge->matrix != nullptr);

  cs_xdef_array_context_t *ac
    = static_cast<cs_xdef_array_context_t *>(source->context);

  assert(ac->stride == 1);
  assert(ac->value_location == cs_flag_primal_vtx);

  /* Retrieve the values of the potential at each cell vertices */

  double *vals = cb->values;
  double *eval = cb->values + cm->n_vc;

  for (int v = 0; v < cm->n_vc; v++)
    vals[v] = ac->values[cm->v_ids[v]];

  /* Multiply these values by a cellwise Hodge operator previously computed */

  cs_sdm_square_matvec(mass_hodge->matrix, vals, eval);

  for (short int v = 0; v < cm->n_vc; v++)
    values[v] += eval[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the contribution for a cell related to a source term and
 *        add it to the given array of values.
 *        Case of a scalar potential defined at primal vertices by an array.
 *        A discrete Hodge operator has to be computed before this call and
 *        given as an input parameter
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_pvsp_by_c2v_array(const cs_xdef_t      *source,
                                 const cs_cell_mesh_t *cm,
                                 cs_real_t             time_eval,
                                 cs_cell_builder_t    *cb,
                                 void                 *input,
                                 double               *values)
{
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  cs_hodge_t  *mass_hodge = (cs_hodge_t *)input;

  assert(values != nullptr && cm != nullptr && cb != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PV));
  assert(mass_hodge != nullptr);
  assert(mass_hodge->matrix != nullptr);

  cs_xdef_array_context_t  *ac =
    (cs_xdef_array_context_t *)source->context;

  assert(ac->stride == 1);
  assert(ac->value_location == cs_flag_primal_vtx);
  assert(ac->adjacency != nullptr);

  double  *eval = cb->values;

  /* Multiply these values by a cellwise Hodge operator previously computed */

  cs_sdm_square_matvec(mass_hodge->matrix,
                       ac->values + ac->adjacency->idx[cm->c_id],
                       eval);

  for (short int v = 0; v < cm->n_vc; v++)
    values[v] += eval[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar density defined at dual cells by a value.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_real_t                  time_eval,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  const cs_real_t *s_values = (const cs_real_t *)source->context;
  const cs_real_t  density_value = s_values[0];

  for (int v = 0; v < cm->n_vc; v++)
    values[v] += density_value * cm->wvc[v] * cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a vector-valued density defined at dual cells by a value.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcvd_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_real_t                  time_eval,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  const cs_real_t *st_vect = (const cs_real_t *)source->context;

  for (int v = 0; v < cm->n_vc; v++)
    for (int k = 0; k < 3; k++)
      values[3*v+k] += st_vect[k] * cm->wvc[v] * cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the contribution for a cell related to a source term and
 *        add it to the given array of values.
 *        Case of a scalar density defined at dual cells by an array defined
 *        at (primal) vertices or dual cells
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_by_pv_array(const cs_xdef_t      *source,
                                const cs_cell_mesh_t *cm,
                                cs_real_t             time_eval,
                                cs_cell_builder_t    *cb,
                                void                 *input,
                                double               *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  const cs_xdef_array_context_t *ac
    = static_cast<const cs_xdef_array_context_t *>(source->context);

  assert(cs_flag_test(ac->value_location, cs_flag_primal_vtx) ||
         cs_flag_test(ac->value_location, cs_flag_dual_cell));

  for (int v = 0; v < cm->n_vc; v++)
    values[v] += ac->values[cm->v_ids[v]] * cm->wvc[v] * cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the contribution for a cell related to a source term and
 *        add it to the given array of values.
 *        Case of a vector-valued density defined at dual cells by an array
 *        defined at (primal) vertices or dual cells
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcvd_by_pv_array(const cs_xdef_t           *source,
                                const cs_cell_mesh_t      *cm,
                                cs_real_t                  time_eval,
                                cs_cell_builder_t         *cb,
                                void                      *input,
                                double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  const cs_xdef_array_context_t *ac
    = static_cast<const cs_xdef_array_context_t *>(source->context);

  assert(cs_flag_test(ac->value_location, cs_flag_primal_vtx));
  for (int v = 0; v < cm->n_vc; v++) {

    const double  vc_coef = cm->wvc[v] * cm->vol_c;
    const double  *ac_v_values = ac->values + 3*cm->v_ids[v];
    for (int k = 0; k < 3; k++)
      values[3*v+k] += ac_v_values[k] * vc_coef;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of a source term in a cell and add it to
 *         the given array of values.
 *         Case of a scalar density defined at dual cells by an array relying
 *         on the c2v connectivity to access the values
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_by_c2v_array(const cs_xdef_t      *source,
                                 const cs_cell_mesh_t *cm,
                                 cs_real_t             time_eval,
                                 cs_cell_builder_t    *cb,
                                 void                 *input,
                                 double               *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  const cs_xdef_array_context_t *ac
    = static_cast<const cs_xdef_array_context_t *>(source->context);
  const cs_adjacency_t  *adj = ac->adjacency;

  assert(cs_flag_test(ac->value_location, cs_flag_dual_cell_byc));
  assert(adj != nullptr);

  const cs_real_t  *_val = ac->values + adj->idx[cm->c_id];

  for (int v = 0; v < cm->n_vc; v++) {
    assert(cm->v_ids[v] == adj->ids[adj->idx[cm->c_id] + v]);
    values[v] += _val[v] * cm->wvc[v] * cm->vol_c;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the contribution for a cell related to a source term and
 *        add it to the given array of values.
 *        Case of a scalar density defined at dual cells by an array defined
 *        at (primal) cells
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_by_pc_array(const cs_xdef_t      *source,
                                const cs_cell_mesh_t *cm,
                                cs_real_t             time_eval,
                                cs_cell_builder_t    *cb,
                                void                 *input,
                                double               *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  const cs_xdef_array_context_t *ac
    = static_cast<const cs_xdef_array_context_t *>(source->context);

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));
  assert(cs_flag_test(ac->value_location, cs_flag_primal_cell));

  const cs_real_t  val_c = ac->values[cm->c_id] * cm->vol_c;
  for (int v = 0; v < cm->n_vc; v++)
    values[v] += val_c * cm->wvc[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar density defined at dual cells by a DoF function.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_by_dof_func(const cs_xdef_t           *source,
                                const cs_cell_mesh_t      *cm,
                                cs_real_t                  time_eval,
                                cs_cell_builder_t         *cb,
                                void                      *input,
                                double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(time_eval);
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  const cs_xdef_dof_context_t  *dc = (cs_xdef_dof_context_t *)source->context;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  /* Up to now this should be the only location allowed */

  assert(cs_flag_test(dc->dof_location, cs_flag_primal_cell));

  /* Call the DoF function to evaluate the function at xc */

  double  cell_eval;
  dc->func(1, &(cm->c_id), true,  /* dense output ? */
           dc->input,
           &cell_eval);

  cell_eval *= cm->vol_c;
  for (int v = 0; v < cm->n_vc; v++)
    values[v] += cell_eval * cm->wvc[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar density defined at dual cells by an analytic
 *         function without any quadrature rule.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_none_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_real_t                  time_eval,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  cs_xdef_analytic_context_t  *ac =
    (cs_xdef_analytic_context_t *)source->context;

  /* Call the DoF function to evaluate the function at xc */

  double  cell_eval = 0.;

  ac->func(time_eval, 1, nullptr, (const cs_real_t *)cm->xc,
           true, /* compacted output ? */
           ac->input,
           &cell_eval);

  cell_eval *= cm->vol_c;
  for (int v = 0; v < cm->n_vc; v++)
    values[v] += cell_eval * cm->wvc[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_bary_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_real_t                  time_eval,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values)
{
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PVQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PFC |
                       CS_FLAG_COMP_FE  | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV  |
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PEC));

  cs_xdef_analytic_context_t  *ac =
    (cs_xdef_analytic_context_t *)source->context;

  /* Compute the barycenter of each portion of dual cells */

  double  *vol_vc = cb->values;
  for (short int v = 0; v < cm->n_vc; v++)
    vol_vc[v] = cm->vol_c * cm->wvc[v];

  /* 1. cell and vertex contributions */

  cs_real_3_t  *xgv = cb->vectors;
  for (short int v = 0; v < cm->n_vc; v++) {
    xgv[v][0] = 0.25 * vol_vc[v] * (cm->xc[0] + cm->xv[3*v]);
    xgv[v][1] = 0.25 * vol_vc[v] * (cm->xc[1] + cm->xv[3*v+1]);
    xgv[v][2] = 0.25 * vol_vc[v] * (cm->xc[2] + cm->xv[3*v+2]);
  }

  /* 2. face and edge contribution */

  for (short int f = 0; f < cm->n_fc; f++) {

    const double  *xf = cm->face[f].center;
    const double  f_coef = 0.25 * cm->pvol_f[f] / cm->face[f].meas;

    for (short int ie = cm->f2e_idx[f]; ie < cm->f2e_idx[f+1]; ie++) {

      const short int e = cm->f2e_ids[ie];
      const double *xe = cm->edge[e].center;
      const double ef_coef = 0.5 * cm->tef[ie] * f_coef;

      double *xgv1 = xgv[cm->e2v_ids[2*e]];
      double *xgv2 = xgv[cm->e2v_ids[2*e+1]];

      for (int k = 0; k < 3; k++) {
        xgv1[k] += ef_coef * (xe[k] + xf[k]);
        xgv2[k] += ef_coef * (xe[k] + xf[k]);
      }

    } // Loop on face edges

  } // Loop on cell faces

  for (short int v = 0; v < cm->n_vc; v++) {
    const double  invvol = 1/vol_vc[v];
    for (int k = 0; k < 3; k++) xgv[v][k] *= invvol;
  }

  /* Call the analytic function to evaluate the function at xgv */

  double  *eval_xgv = vol_vc + cm->n_vc;
  ac->func(time_eval, cm->n_vc, nullptr, (const cs_real_t *)xgv,
           true, /* compacted output ? */
           ac->input,
           eval_xgv);

  for (short int v = 0; v < cm->n_vc; v++)
    values[v] = vol_vc[v] * eval_xgv[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use a the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q1o1_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_real_t                  time_eval,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PFC | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  cs_xdef_analytic_context_t  *ac =
    (cs_xdef_analytic_context_t *)source->context;

  for (short int f = 0; f < cm->n_fc; f++) {

    cs_real_3_t  xg[2], xfc;
    cs_real_t  eval_xg[2];

    const double  *xf = cm->face[f].center;
    const double  hf_coef = 0.5 * cm->pvol_f[f]/cm->face[f].meas;

    for (int k = 0; k < 3; k++) xfc[k] = 0.25*(xf[k] + cm->xc[k]);

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];
      const double  *xv1 = cm->xv + 3*v1, *xv2 = cm->xv + 3*v2;

      /* xg = 0.25(xv1 + xe + xf + xc) where xe = 0.5*(xv1 + xv2) */
      for (int k = 0; k < 3; k++)
        xg[0][k] = xfc[k] + 0.375*xv1[k] + 0.125*xv2[k];

      /* xg = 0.25(xv2 + xe + xf + xc) where xe = 0.5*(xv1 + xv2) */
      for (int k = 0; k < 3; k++)
        xg[1][k] = xfc[k] + 0.375*xv2[k] + 0.125*xv1[k];

      ac->func(time_eval, 2, nullptr, (const cs_real_t *)xg,
               true, /* compacted output ? */
               ac->input,
               eval_xg);

      const double  half_pef_vol = cm->tef[i]*hf_coef;
      values[v1] += half_pef_vol * eval_xg[0];
      values[v2] += half_pef_vol * eval_xg[1];

    } /* Loop on face edges */

  } /* Loop on cell faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use a ten-point quadrature rule to evaluate the integral.
 *         Exact for quadratic function.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q10o2_by_analytic(const cs_xdef_t           *source,
                                      const cs_cell_mesh_t      *cm,
                                      cs_real_t                  time_eval,
                                      cs_cell_builder_t         *cb,
                                      void                      *input,
                                      double                    *values)
{
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cb != nullptr);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PFC | CS_FLAG_COMP_FE  |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV  | CS_FLAG_COMP_PVQ |
                       CS_FLAG_COMP_PEQ));

  cs_xdef_analytic_context_t *ac =
    (cs_xdef_analytic_context_t *)source->context;

  /* Temporary buffers */

  double  *contrib = cb->values; /* size n_vc */

  /* 1) Compute the contributions seen by the whole portion of dual cell */

  /* Cell evaluation */

  double  eval_c;
  ac->func(time_eval,
           1,
           nullptr,
           cm->xc,
           true, /* compacted output ? */
           ac->input,
           &eval_c);

  /* Contributions related to vertices */

  double  *eval_v = cb->values + cm->n_vc; /* size n_vc */
  ac->func(time_eval,
           cm->n_vc,
           nullptr,
           cm->xv,
           true, /* compacted output ? */
           ac->input,
           eval_v);

  cs_real_3_t  *xvc = cb->vectors;
  for (short int v = 0; v < cm->n_vc; v++) {
    const double  *xv = cm->xv + 3*v;
    for (int k = 0; k < 3; k++) xvc[v][k] = 0.5*(cm->xc[k] + xv[k]);
  }

  double  *eval_vc = cb->values + 2*cm->n_vc; /* size n_vc */
  ac->func(time_eval,
           cm->n_vc,
           nullptr,
           (const cs_real_t *)xvc,
           true, /* compacted output ? */
           ac->input,
           eval_vc);

  for (short int v = 0; v < cm->n_vc; v++) {

    /* Set the initial values
       -1/20 on extremity points and 1/5 on midpoints */

    contrib[v] = cm->wvc[v] * cm->vol_c
                 * (-0.05 * (eval_c + eval_v[v]) + 0.2 * eval_vc[v]);

  } /* Loop on vertices */

  /* 2) Compute the contribution related to edge
     The portion of dual cell seen by each vertex is 1/2 |pec| */

  cs_real_3_t  *x_e = cb->vectors;            // size = n_ec (overwrite xvc)
  cs_real_3_t  *x_ec = cb->vectors + cm->n_ec; // size = n_ec (overwrite xvc)

  for (short int e = 0; e < cm->n_ec; e++) {
    for (int k = 0; k < 3; k++) {
      x_e[e][k] = cm->edge[e].center[k];
      x_ec[e][k] = 0.5*(cm->xc[k] + x_e[e][k]);
    }
  }

  /* Evaluate the analytic function at xe and xec */

  double  *eval_e = cb->values + cm->n_vc; // size=n_ec (overwrite eval_v)
  double  *eval_ec = eval_e + cm->n_ec;    // size=n_ec (overwrite eval_vc)
  ac->func(time_eval,
           2 * cm->n_ec,
           nullptr,
           (const cs_real_t *)x_e, // x_e (n_ec) then xec (n_ec)
           true, /* compacted output ? */
           ac->input,
           eval_e);

  /* Evaluation at x_ve (size = 2*n_ec) */

  cs_real_3_t  *x_ve = cb->vectors; // size=2*n_ec (overwrite xe and xec)
  for (short int e = 0; e < cm->n_ec; e++) {

    const double *xe = cm->edge[e].center;
    const double *xv1 = cm->xv + 3*(cm->e2v_ids[2*e]);
    const double *xv2 = cm->xv + 3*(cm->e2v_ids[2*e+1]);

    for (int k = 0; k < 3; k++) {
      x_ve[2*e  ][k] = 0.5*(xv1[k] + xe[k]);
      x_ve[2*e+1][k] = 0.5*(xv2[k] + xe[k]);
    }

  } /* Loop on edges */

  double  *eval_ve = eval_ec + cm->n_ec; // size = 2*n_ec
  ac->func(time_eval,
           2 * cm->n_ec,
           nullptr,
           (const cs_real_t *)x_ve,
           true, /* compacted output ? */
           ac->input,
           eval_ve);

  /* 3) Main loop on faces */

  double  *pvf_vol = eval_ve + 2*cm->n_ec;  /* size n_vc */

  for (short int f = 0; f < cm->n_fc; f++) {

    const double  *x_f = cm->face[f].center;
    const double  hf_coef = 0.5* cm->pvol_f[f]/cm->face[f].meas;

    /* Reset volume of the face related to a vertex */

    for (short int v = 0; v < cm->n_vc; v++) pvf_vol[v] = 0;

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];
      const double  half_pef_vol = hf_coef * cm->tef[i];

      pvf_vol[v1] += half_pef_vol;
      pvf_vol[v2] += half_pef_vol;

      cs_real_3_t  x_ef;
      cs_real_t  eval_ef;
      for (int k = 0; k < 3; k++) x_ef[k] = 0.5*(cm->edge[e].center[k] + x_f[k]);
      ac->func(time_eval,
               1,
               nullptr,
               x_ef,
               true, /* compacted output ? */
               ac->input,
               &eval_ef);

      /* 1/5 (EF + EC) -1/20 * (E) */

      const double  common_ef_contrib =
        0.2*(eval_ef + eval_ec[e]) -0.05*eval_e[e];

      contrib[v1] += half_pef_vol*(common_ef_contrib + 0.2*eval_ve[2*e]);
      contrib[v2] += half_pef_vol*(common_ef_contrib + 0.2*eval_ve[2*e+1]);

    } /* Loop on face edges */

    /* Contributions related to this face */

    cs_real_3_t  *x_vfc = cb->vectors;  // size=2+n_vc (overwrite x_ve)
    for (int k = 0; k < 3; k++) {
      x_vfc[0][k] = x_f[k];
      x_vfc[1][k] = 0.5*(x_f[k] + cm->xc[k]);  /* x_fc */
    }

    short int  n_vf = 0;
    for (short int v = 0; v < cm->n_vc; v++) {
      if (pvf_vol[v] > 0) {
        cb->ids[n_vf] = v;
        for (int k = 0; k < 3; k++)
          x_vfc[2+n_vf][k] = 0.5*(x_f[k] + cm->xv[3*v+k]);
        n_vf++;
      }
    }

    double  *eval_vfc = pvf_vol + cm->n_vc; /* size=n_vf + 2 */
    ac->func(time_eval,
             2 + n_vf,
             nullptr,
             (const cs_real_t *)x_vfc,
             true, /* compacted output ? */
             ac->input,
             eval_vfc);

    for (short int i = 0; i < n_vf; i++) {
      short int  v = cb->ids[i];
      const double  val_vfc = -0.05*eval_vfc[0] + 0.2*eval_vfc[1];
      contrib[v] += pvf_vol[v] * (val_vfc + 0.2*eval_vfc[2+i]);
    }

  } /* Loop on cell faces */

  /* Add the computed contributions to the return values */

  for (short int v = 0; v < cm->n_vc; v++)
    values[v] += contrib[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use a five-point quadrature rule to evaluate the integral.
 *         Exact for cubic function.
 *         This function may be expensive since many evaluations are needed.
 *         Please use it with care.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q5o3_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_real_t                  time_eval,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values)
{
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  double sum, weights[5], results[5];
  cs_real_3_t gauss_pts[5];

  assert(values != nullptr && cm != nullptr);
  assert(cb != nullptr);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_EV  | CS_FLAG_COMP_PFC | CS_FLAG_COMP_FEQ));

  cs_xdef_analytic_context_t *ac =
    (cs_xdef_analytic_context_t *)source->context;

  /* Temporary buffers */

  double *contrib = cb->values;
  memset(contrib, 0, cm->n_vc*sizeof(double));

  /* Main loop on faces */

  for (short int f = 0; f < cm->n_fc; f++) {

    const double *x_f = cm->face[f].center;
    const double hf_coef = 0.5* cm->pvol_f[f]/cm->face[f].meas;

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int e = cm->f2e_ids[i];
      const short int v1 = cm->e2v_ids[2*e];
      const short int v2 = cm->e2v_ids[2*e+1];
      const double half_pef_vol = hf_coef * cm->tef[i];

      /* Compute Gauss points and its weights */

      cs_quadrature_tet_5pts(cm->xv + 3*v1, cm->edge[e].center, x_f, cm->xc,
                             half_pef_vol,
                             gauss_pts, weights);

      ac->func(time_eval,
               5,
               nullptr,
               (const cs_real_t *)gauss_pts,
               true, /* compacted output ? */
               ac->input,
               results);

      sum = 0.;
      for (int p = 0; p < 5; p++) sum += results[p] * weights[p];
      contrib[v1] += sum;

      /* Compute Gauss points and its weights */

      cs_quadrature_tet_5pts(cm->xv + 3*v2, cm->edge[e].center, x_f, cm->xc,
                             half_pef_vol,
                             gauss_pts, weights);

      ac->func(time_eval,
               5,
               nullptr,
               (const cs_real_t *)gauss_pts,
               true, /* compacted output ? */
               ac->input,
               results);

      sum = 0.;
      for (int p = 0; p < 5; p++) sum += results[p] * weights[p];
      contrib[v2] += sum;

    } /* Loop on face edges */

  } /* Loop on cell faces */

  /* Add the computed contributions to the return values */
  for (short int v = 0; v < cm->n_vc; v++)
    values[v] += contrib[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar potential defined at primal vertices and cells
 *         by a constant value.
 *         A discrete Hodge operator has to be computed before this call and
 *         given as input parameter
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_vcsp_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_real_t                  time_eval,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values)
{
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  cs_hodge_t  *mass_hodge = (cs_hodge_t *)input;

  assert(values != nullptr && cm != nullptr && cb != nullptr);
  assert(mass_hodge != nullptr);
  assert(mass_hodge->matrix != nullptr);

  const cs_real_t *s_values = (const cs_real_t *)source->context;
  const cs_real_t  pot_value = s_values[0];

  /* Retrieve the values of the potential at each cell vertices */

  double  *eval = cb->values;
  for (short int v = 0; v < cm->n_vc; v++)
    eval[v] = pot_value;
  eval[cm->n_vc] = pot_value;

  /* Multiply these values by a cellwise Hodge operator previously computed */

  double  *hdg_eval = cb->values + cm->n_vc + 1;
  cs_sdm_square_matvec(mass_hodge->matrix, eval, hdg_eval);

  for (short int v = 0; v < cm->n_vc + 1; v++)
    values[v] += hdg_eval[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar potential defined at primal vertices and cells by
 *         an analytical function.
 *         A discrete Hodge operator has to be computed before this call and
 *         given as input parameter
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_vcsp_by_analytic(const cs_xdef_t           *source,
                                const cs_cell_mesh_t      *cm,
                                cs_real_t                  time_eval,
                                cs_cell_builder_t         *cb,
                                void                      *input,
                                double                    *values)
{
  if (source == nullptr)
    return;

  cs_hodge_t  *mass_hodge = (cs_hodge_t *)input;

  assert(values != nullptr && cm != nullptr && cb != nullptr);
  assert(mass_hodge != nullptr);
  assert(mass_hodge->matrix != nullptr);

  cs_xdef_analytic_context_t  *ac =
    (cs_xdef_analytic_context_t *)source->context;

  /* Retrieve the values of the potential at each cell vertices */

  double  *eval = cb->values;

  ac->func(time_eval,
           cm->n_vc,
           nullptr,
           cm->xv,
           true, /* compacted output ? */
           ac->input,
           eval);

  /* Retrieve the value at the cell center */

  ac->func(time_eval,
           1,
           nullptr,
           cm->xc,
           true, /* compacted output ? */
           ac->input,
           eval + cm->n_vc);

  /* Multiply these values by a cellwise Hodge operator previously computed */

  double  *hdg_eval = cb->values + cm->n_vc + 1;
  cs_sdm_square_matvec(mass_hodge->matrix, eval, hdg_eval);

  for (short int v = 0; v < cm->n_vc + 1; v++)
    values[v] += hdg_eval[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar density (sd) defined on primal cells by a value.
 *         Case of face-based/cell-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fcb_pcsd_by_value(const cs_xdef_t           *source,
                                 const cs_cell_mesh_t      *cm,
                                 cs_real_t                  time_eval,
                                 cs_cell_builder_t         *cb,
                                 void                      *input,
                                 double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);

  const cs_real_t *s_values = (const cs_real_t *)source->context;

  values[cm->n_fc] += s_values[0] * cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a density defined on primal cells by a DoF function.
 *         Case of scalar-valued face-based schemes.
 *         Case of face-based and cell-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fcb_pcsd_by_dof_func(const cs_xdef_t           *source,
                                    const cs_cell_mesh_t      *cm,
                                    cs_real_t                  time_eval,
                                    cs_cell_builder_t         *cb,
                                    void                      *input,
                                    double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(time_eval);
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  cs_xdef_dof_context_t  *cx = (cs_xdef_dof_context_t *)source->context;

  assert(values != nullptr && cm != nullptr);

  /* Up to now this should be the only location allowed */

  assert(cs_flag_test(cx->dof_location, cs_flag_primal_cell));

  /* Call the DoF function to evaluate the function at xc */

  double  cell_eval;
  cx->func(1, &(cm->c_id), true,  /* dense output ? */
           cx->input,
           &cell_eval);

  values[cm->n_fc] += cell_eval * cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a density defined on primal cells by a DoF function.
 *         Case of vector-valued face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fb_pcvd_by_dof_func(const cs_xdef_t           *source,
                                   const cs_cell_mesh_t      *cm,
                                   cs_real_t                  time_eval,
                                   cs_cell_builder_t         *cb,
                                   void                      *input,
                                   double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(time_eval);
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  cs_xdef_dof_context_t  *cx = (cs_xdef_dof_context_t *)source->context;

  assert(values != nullptr && cm != nullptr);

  /* Up to now this should be the only location allowed */

  assert(cs_flag_test(cx->dof_location, cs_flag_primal_cell));

  /* Call the DoF function to evaluate the function at xc */

  cs_real_t  cell_eval[3];
  cx->func(1, &(cm->c_id), true,  /* compacted output ? */
           cx->input,
           cell_eval);

  for (int k = 0; k < 3; k++)
    values[3*cm->n_fc+k] += cell_eval[k] * cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a vector-valued density (vd) defined on primal cells
 *         by a value.
 *         Case of face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fb_pcvd_by_value(const cs_xdef_t           *source,
                                const cs_cell_mesh_t      *cm,
                                cs_real_t                  time_eval,
                                cs_cell_builder_t         *cb,
                                void                      *input,
                                double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);

  const cs_real_t *v_input = (const cs_real_t *)source->context;
  for (int k = 0; k < source->dim; k++)
    values[source->dim*cm->n_fc + k] = v_input[k] * cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar density defined at primal cells by an analytical
 *         function.
 *         Use the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *         Case of face-based and cell-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in]      input      null or pointer to a structure cast on-the-fly
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fcb_pcsd_bary_by_analytic(const cs_xdef_t           *source,
                                         const cs_cell_mesh_t      *cm,
                                         cs_real_t                  time_eval,
                                         cs_cell_builder_t         *cb,
                                         void                      *input,
                                         double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);

  cs_xdef_analytic_context_t  *ac =
    (cs_xdef_analytic_context_t *)source->context;

  /* Call the analytic function to evaluate the function at xc */

  double  eval_xc;
  ac->func(time_eval,
           1,
           nullptr,
           (const cs_real_t *)cm->xc,
           true, /* compacted output ? */
           ac->input,
           &eval_xc);

  values[cm->n_fc] += cm->vol_c * eval_xc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of a source term for a cell and add it to
 *         the given array of values.
 *         Case of a scalar density (sd) defined on primal cells by an analytic
 *         function.
 *         Case of face-based and cell-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fcb_pcsd_by_analytic(const cs_xdef_t           *source,
                                    const cs_cell_mesh_t      *cm,
                                    cs_real_t                  time_eval,
                                    cs_cell_builder_t         *cb,
                                    void                      *input,
                                    double                    *values)
{
  CS_NO_WARN_IF_UNUSED(input);
  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));
  assert(source->dim == 1);

  if (source->qtype == CS_QUADRATURE_BARY ||
      source->qtype == CS_QUADRATURE_NONE)
    cs_source_term_fcb_pcsd_bary_by_analytic(source, cm, time_eval, cb, input,
                                             values);

  else {

    constexpr cs_real_t c_1ov3 = 1./3.;
    const cs_real_t *xv = cm->xv;

    double  cell_values  = 0.0;

    cs_quadrature_tetra_integral_t  *qfunc =
      cs_quadrature_get_tetra_integral(1, source->qtype);

    cs_xdef_analytic_context_t *ac =
      (cs_xdef_analytic_context_t *)source->context;

    /* Switch according to the cell type: optimised version for tetra */

    switch (cm->type) {

    case FVM_CELL_TETRA:
    {
      assert(cm->n_fc == 4 && cm->n_vc == 4);
      qfunc(time_eval, xv, xv+3, xv+6, xv+9, cm->vol_c, ac->func, ac->input,
            &cell_values);
    } break;

    case FVM_CELL_PYRAM:
    case FVM_CELL_PRISM:
    case FVM_CELL_HEXA:
    case FVM_CELL_POLY:
    {
      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_quant_t  pfq = cm->face[f];
        const double  hf_coef = c_1ov3 * cm->hfc[f];
        const int  start = cm->f2e_idx[f];
        const int  end = cm->f2e_idx[f+1];
        const short int  n_vf = end - start;  /* #vertices (=#edges) */
        const short int  *f2e_ids = cm->f2e_ids + start;

        assert(n_vf > 2);
        switch(n_vf){

        case 3: /* triangle (optimized version, no subdivision) */
        {
          short int v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          qfunc(time_eval, cm->xc, xv+3*v0, xv+3*v1, xv+3*v2,
                hf_coef*pfq.meas,
                ac->func, ac->input,
                &cell_values);
        } break;

        default:
          {
          const double *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            /* Edge-related variables */

            const short int e0  = f2e_ids[e];
            const double   *xv0 = xv + 3 * cm->e2v_ids[2 * e0];
            const double   *xv1 = xv + 3 * cm->e2v_ids[2 * e0 + 1];

            qfunc(time_eval, cm->xc, pfq.center, xv0, xv1,
                  hf_coef*tef[e], ac->func, ac->input, &cell_values);
          }
        } break;

        } /* End of switch */

      } /* End of loop on faces */
    }
    break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Unknown cell-type.\n", __func__);
      break;

    } /* End of switch on the cell-type */

    values[cm->n_fc] += cell_values;

  }  /* If not a barycentric quadrature */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar density defined at primal cells by an array.
 *         Case of face-based and cell-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fcb_pcsd_by_array(const cs_xdef_t           *source,
                                 const cs_cell_mesh_t      *cm,
                                 cs_real_t                  time_eval,
                                 cs_cell_builder_t         *cb,
                                 void                      *input,
                                 double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);

  const cs_xdef_array_context_t *cx
    = static_cast<const cs_xdef_array_context_t *>(source->context);

  assert(cs_flag_test(cx->value_location, cs_flag_primal_cell));
  values[cm->n_fc] += cx->values[cm->c_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a vector-valued density defined at primal cells by an
 *         analytical function.
 *         Use the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *         Case of face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in]      input      null or pointer to a structure cast on-the-fly
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fb_pcvd_bary_by_analytic(const cs_xdef_t           *source,
                                        const cs_cell_mesh_t      *cm,
                                        cs_real_t                  time_eval,
                                        cs_cell_builder_t         *cb,
                                        void                      *input,
                                        double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(source->dim == 3);

  cs_xdef_analytic_context_t  *ac =
    (cs_xdef_analytic_context_t *)source->context;

  /* Call the analytic function to evaluate the function at xc */

  double  eval_xc[3];
  ac->func(time_eval,
           1,
           nullptr,
           (const cs_real_t *)cm->xc,
           true, /* compacted output ? */
           ac->input,
           eval_xc);

  cs_real_t  *c_val = values + 3*cm->n_fc;
  for (int k = 0; k < source->dim; k++)
    c_val[k] += cm->vol_c * eval_xc[k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a vector-valued density defined at primal cells by an
 *         analytical function.
 *         Use the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *         Case of MAC face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in]      input      null or pointer to a structure cast on-the-fly
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_macfb_pcvd_by_analytic(const cs_xdef_t      *source,
                                      const cs_cell_mesh_t *cm,
                                      cs_real_t             time_eval,
                                      cs_cell_builder_t    *cb,
                                      void                 *input,
                                      double               *values)
{
  CS_NO_WARN_IF_UNUSED(cb);

  const cs_macfb_builder_t *macb = (const cs_macfb_builder_t *)input;

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr && macb != nullptr);
  assert(source->dim == 3);

  cs_quadrature_hexa_integral_t *qfunc
    = cs_quadrature_get_hexa_integral(3, source->qtype);

  cs_xdef_analytic_context_t *ac
    = (cs_xdef_analytic_context_t *)source->context;

  switch (cm->type) {

  case FVM_CELL_HEXA: {
    assert(cm->n_fc == 6);
  } break;

  case FVM_CELL_TETRA:
  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_POLY:
  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Unknown cell-type.\n", __func__);
    break;

  } /* End of switch on the cell-type */

  /* Volume dual cell */
  cs_real_3_t cm_h;

  for (short int f = 0; f < cm->n_fc; f++) {

    cm_h[macb->f_axis[f]] = 2.0 * cm->hfc[f];

  } /* Loop on cell faces */

  for (short int f = 0; f < cm->n_fc; f++) {

    /* Face info */
    const cs_quant_t fq = cm->face[f];

    /* Barycenter of the dual cell */
    cs_real_3_t dc_xc;
    cs_math_3_average(cm->xc, fq.center, dc_xc);

    /* Size of dual cell */
    cs_real_3_t dc_h = { cm_h[0], cm_h[1], cm_h[2] };
    dc_h[macb->f_axis[f]] /= 2.0;

    /* Call the quadrature function to integrate the function */
    cs_real_t int_dc[3] = { 0.0, 0.0, 0.0 };
    qfunc(
      time_eval, dc_xc, dc_h[0], dc_h[1], dc_h[2], ac->func, ac->input, int_dc);

    values[f] = int_dc[macb->f_axis[f]] / macb->f_vol_cv[f];

  } /* Loop on cell faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of a source term for a cell and add it to
 *         the given array of values.
 *         Case of a vector density (vd) defined on primal cells by an analytic
 *         function.
 *         Case of face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fb_pcvd_by_analytic(const cs_xdef_t           *source,
                                   const cs_cell_mesh_t      *cm,
                                   cs_real_t                  time_eval,
                                   cs_cell_builder_t         *cb,
                                   void                      *input,
                                   double                    *values)
{
  CS_NO_WARN_IF_UNUSED(input);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));
  assert(source->dim == 3);

  if (source->qtype == CS_QUADRATURE_BARY)
    cs_source_term_fb_pcvd_bary_by_analytic(source, cm, time_eval, cb, input,
                                         values);

  else {

    constexpr cs_real_t c_1ov3 = 1./3.;
    const cs_real_t *xv = cm->xv;

    cs_real_3_t  cell_values = {0.0, 0.0, 0.0};

    cs_quadrature_tetra_integral_t  *qfunc =
      cs_quadrature_get_tetra_integral(3, source->qtype);

    cs_xdef_analytic_context_t  *ac =
      (cs_xdef_analytic_context_t *)source->context;

    /* Switch according to the cell type: optimized version for tetra */

    switch (cm->type) {

    case FVM_CELL_TETRA:
    {
      assert(cm->n_fc == 4 && cm->n_vc == 4);
            qfunc(time_eval, xv, xv+3, xv+6, xv+9, cm->vol_c,
            ac->func, ac->input,
            cell_values);

    } break;

    case FVM_CELL_PYRAM:
    case FVM_CELL_PRISM:
    case FVM_CELL_HEXA:
    case FVM_CELL_POLY:
    {
      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_quant_t  pfq = cm->face[f];
        const double  hf_coef = c_1ov3 * cm->hfc[f];
        const int  start = cm->f2e_idx[f];
        const int  end = cm->f2e_idx[f+1];
        const short int  n_vf = end - start; // #vertices (=#edges)
        const short int  *f2e_ids = cm->f2e_ids + start;

        assert(n_vf > 2);
        switch(n_vf){

        case 3: /* triangle (optimized version, no subdivision) */
        {
          short int v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          qfunc(time_eval, cm->xc, xv+3*v0, xv+3*v1, xv+3*v2,
                hf_coef*pfq.meas,
                ac->func, ac->input,
                cell_values);
        } break;

        default:
          {
          const double *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            /* Edge-related variables */

            const short int e0  = f2e_ids[e];
            const double   *xv0 = xv + 3 * cm->e2v_ids[2 * e0];
            const double   *xv1 = xv + 3 * cm->e2v_ids[2 * e0 + 1];

            qfunc(time_eval, cm->xc, pfq.center, xv0, xv1,
                  hf_coef*tef[e],
                  ac->func, ac->input,
                  cell_values);
          }
        } break;

        } /* End of switch */

      } /* End of loop on faces */

    }
    break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Unknown cell-type.\n", __func__);
      break;

    } /* End of switch on the cell-type */

    cs_real_t *c_val = values + 3*cm->n_fc;
    c_val[0] += cell_values[0];
    c_val[1] += cell_values[1];
    c_val[2] += cell_values[2];

  }  /* If not a barycentric quadrature */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of a source term for a cell and add it to
 *         the given array of values.
 *         Case of a vector density (vd) defined on primal cells by an array
 *         Case of CDO face-based schemes.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fb_pcvd_by_array(const cs_xdef_t           *source,
                                const cs_cell_mesh_t      *cm,
                                cs_real_t                  time_eval,
                                cs_cell_builder_t         *cb,
                                void                      *input,
                                double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(source->dim == 3);

  const cs_xdef_array_context_t *cx
    = static_cast<const cs_xdef_array_context_t *>(source->context);
  const double  *arr = cx->values + 3*cm->c_id;

  assert(cs_flag_test(cx->value_location, cs_flag_primal_cell));

  double  *val_c = values + 3*cm->n_fc;
  val_c[0] += arr[0];
  val_c[1] += arr[1];
  val_c[2] += arr[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of a source term for a cell and add it to
 *         the given array of values.
 *         Case of a scalar density (sd) defined on primal cells by a value.
 *         Case of HHO schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_hhosd_by_value(const cs_xdef_t           *source,
                              const cs_cell_mesh_t      *cm,
                              cs_real_t                  time_eval,
                              cs_cell_builder_t         *cb,
                              void                      *input,
                              double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr && input != nullptr);

  constexpr cs_real_t c_1ov3 = 1./3.;

  const cs_real_t  *const_val = (const cs_real_t *)source->context;

  cs_hho_builder_t  *hhob = (cs_hho_builder_t *)input;
  cs_real_t  *cell_values = values + cm->n_fc * hhob->face_basis[0]->size;

  const cs_basis_func_t  *cbf = hhob->cell_basis;

  switch (cbf->poly_order) {

  case 0:
  case 1:
    cbf->eval_all_at_point(cbf, cm->xc, cell_values);
    for (int i = 0; i < cbf->size; i++)
      cell_values[i] *= cm->vol_c * const_val[0];
    break;

  default:
    /* Reset cell values */

    memset(cell_values, 0, sizeof(cs_real_t)*cbf->size);

    /* Switch according to the cell type: optimised version for tetra */

    switch (cm->type) {

    case FVM_CELL_TETRA:
      assert(cm->n_fc == 4 && cm->n_vc == 4);
      _hho_add_tetra_by_val(const_val[0], cbf,
                            cm->xv, cm->xv + 3, cm->xv + 6, cm->xv + 9,
                            cm->vol_c, cb, cell_values);
      break;

    case FVM_CELL_PYRAM:
    case FVM_CELL_PRISM:
    case FVM_CELL_HEXA:
    case FVM_CELL_POLY:
      {
      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_quant_t pfq     = cm->face[f];
        const double     hf_coef = c_1ov3 * cm->hfc[f];
        const int        start   = cm->f2e_idx[f];
        const int        end     = cm->f2e_idx[f + 1];
        const short int  n_vf    = end - start; /* #vertices (=#edges) */
        const short int *f2e_ids = cm->f2e_ids + start;

        assert(n_vf > 2);
        switch (n_vf) {

        case CS_TRIANGLE_CASE: /* Optimized version, no subdivision */
        {
          short int v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          _hho_add_tetra_by_val(const_val[0], cbf,
                                cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2,
                                cm->xc,
                                hf_coef * pfq.meas, cb, cell_values);
        } break;

        default: {
          const double *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            /* Edge-related variables */

            const short int e0  = f2e_ids[e];
            const double   *xv0 = cm->xv + 3 * cm->e2v_ids[2 * e0];
            const double   *xv1 = cm->xv + 3 * cm->e2v_ids[2 * e0 + 1];

            _hho_add_tetra_by_val(const_val[0], cbf,
                                  xv0, xv1, pfq.center, cm->xc,
                                  hf_coef*tef[e], cb, cell_values);
          }
        } break;

        } /* End of switch */

      } /* End of loop on faces */
    } break;

    default:
      bft_error(__FILE__, __LINE__, 0,  _(" Unknown cell-type.\n"));
      break;

    } /* End of switch on the cell-type */
    break;

  } /* Switch on polynomial order */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of a source term for a cell and add it to
 *         the given array of values.
 *         Case of a scalar density (sd) defined on primal cells by an analytic
 *         function.
 *         Case of HHO schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_hhosd_by_analytic(const cs_xdef_t           *source,
                                 const cs_cell_mesh_t      *cm,
                                 cs_real_t                  time_eval,
                                 cs_cell_builder_t         *cb,
                                 void                      *input,
                                 double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr && input != nullptr);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  constexpr cs_real_t c_1ov3 = 1./3.;

  cs_hho_builder_t  *hhob = (cs_hho_builder_t *)input;
  cs_xdef_analytic_context_t  *ac =
    (cs_xdef_analytic_context_t *)source->context;
  cs_real_t  *cell_values = values + cm->n_fc * hhob->face_basis[0]->size;

  const cs_basis_func_t  *cbf = hhob->cell_basis;

  /* Reset cell values */

  memset(cell_values, 0, sizeof(cs_real_t)*cbf->size);

  /* Switch according to the cell type: optimised version for tetra */

  switch (cm->type) {

  case FVM_CELL_TETRA:
    assert(cm->n_fc == 4 && cm->n_vc == 4);
    _hho_add_tetra_by_ana(ac, cbf,
                          cm->xv, cm->xv+3, cm->xv+6, cm->xv+9,
                          cm->vol_c, time_eval,
                          cb, cell_values);
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:
    for (short int f = 0; f < cm->n_fc; f++) {

      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = c_1ov3 * cm->hfc[f];
      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_vf = end - start; /* #vertices (=#edges) */
      const short int *f2e_ids = cm->f2e_ids + start;

      assert(n_vf > 2);
      switch(n_vf){

      case CS_TRIANGLE_CASE: /* triangle (optimized version, no subdivision) */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          _hho_add_tetra_by_ana(ac, cbf,
                                cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2, cm->xc,
                                hf_coef * pfq.meas, time_eval,
                                cb, cell_values);
        }
        break;

      default:
        {
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            /* Edge-related variables */
            const short int e0  = f2e_ids[e];
            const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
            const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

            _hho_add_tetra_by_ana(ac, cbf,
                                  xv0, xv1, pfq.center, cm->xc,
                                  hf_coef*tef[e], time_eval,
                                  cb, cell_values);
          }
        }
        break;

      } /* End of switch */

    } /* End of loop on faces */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,  _(" Unknown cell-type.\n"));
    break;
  } /* End of switch on the cell-type */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of a source term for a cell and add it to
 *         the given array of values.
 *         Case of a vector field (vd) defined on primal cells by an analytic
 *         function.
 *         Case of HHO schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_hhovd_by_analytic(const cs_xdef_t           *source,
                                 const cs_cell_mesh_t      *cm,
                                 cs_real_t                  time_eval,
                                 cs_cell_builder_t         *cb,
                                 void                      *input,
                                 double                    *values)
{
  CS_NO_WARN_IF_UNUSED(cb);
  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr && input != nullptr);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  constexpr cs_real_t c_1ov3 = 1./3.;

  cs_hho_builder_t  *hhob = (cs_hho_builder_t *)input;
  cs_xdef_analytic_context_t  *ac =
    (cs_xdef_analytic_context_t *)source->context;
  cs_real_t  *cell_values = values + 3*cm->n_fc * hhob->face_basis[0]->size;

  const cs_basis_func_t  *cbf = hhob->cell_basis;

  /* Reset cell values */

  memset(cell_values, 0, 3*sizeof(cs_real_t)*cbf->size);

  /* Switch according to the cell type: optimised version for tetra */

  switch (cm->type) {

  case FVM_CELL_TETRA:
    assert(cm->n_fc == 4 && cm->n_vc == 4);
    _hho_add_tetra_by_ana_vd(ac, cbf,
                             cm->xv, cm->xv+3, cm->xv+6, cm->xv+9,
                             cm->vol_c, time_eval,
                             cb, cell_values);
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:
    for (short int f = 0; f < cm->n_fc; f++) {

      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = c_1ov3 * cm->hfc[f];
      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_vf = end - start; /* #vertices (=#edges) */
      const short int *f2e_ids = cm->f2e_ids + start;

      assert(n_vf > 2);
      switch(n_vf){

      case 3: /* triangle (optimized version, no subdivision) */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          _hho_add_tetra_by_ana_vd(ac, cbf,
                                   cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2,
                                   cm->xc, hf_coef * pfq.meas, time_eval,
                                   cb, cell_values);
        }
        break;

      default:
        {
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            /* Edge-related variables */
            const short int e0  = f2e_ids[e];
            const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
            const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

            _hho_add_tetra_by_ana_vd(ac, cbf,
                                     xv0, xv1, pfq.center, cm->xc,
                                     hf_coef*tef[e], time_eval,
                                     cb, cell_values);
          }
        }
        break;

      } /* End of switch */

    } /* End of loop on faces */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,  _(" Unknown cell-type.\n"));
    break;
  } /* End of switch on the cell-type */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it to the given array of values.
 *         Case of a scalar flux defined at dual faces by a constant value.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or nullptr)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dfsf_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_real_t                  time_eval,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values)
{
  CS_NO_WARN_IF_UNUSED(input);
  CS_NO_WARN_IF_UNUSED(cb);
  CS_NO_WARN_IF_UNUSED(time_eval);

  if (source == nullptr)
    return;

  assert(values != nullptr && cm != nullptr);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_DFQ));

  const cs_real_t *vector = (const cs_real_t *)source->context;

  /* Retrieve the values of the normal flux for each dual face */

  for (short int e = 0; e < cm->n_ec; e++)
    values[e] = cm->dface[e].meas *
      cs_math_3_dot_product(vector, cm->dface[e].unitv);
}

/*----------------------------------------------------------------------------*/
END_C_DECLS
