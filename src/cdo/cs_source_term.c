/*============================================================================
 * Functions and structures to deal with source term computations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_defs.h"

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

#include <bft_mem.h>

#include "cs_basis_func.h"
#include "cs_evaluate.h"
#include "cs_hho_builder.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_source_term.h"

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
static const cs_time_step_t  *cs_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a name (copy or generic name)
 *
 * \param[in] name       input name
 * \param[in] basename   generic name by default if input name is NULL
 * \param[in] id         id related to this name
 *
 * \return a pointer to a new allocated string
 */
/*----------------------------------------------------------------------------*/

static inline char *
_get_name(const char   *name,
          const char   *base_name,
          int           id)
{
  char *n = NULL;

  if (name == NULL) { /* Define a name by default */
    assert(id < 100);
    int len = strlen(base_name) + 4; // 1 + 3 = "_00\n"
    BFT_MALLOC(n, len, char);
    sprintf(n, "%s_%2d", base_name, id);
  }
  else {  /* Copy name */
    int  len = strlen(name) + 1;
    BFT_MALLOC(n, len, char);
    strncpy(n, name, len);
  }

  return n;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the mask associated to each cell from the mask related to
 *         the given source term structure
 *
 * \param[in]      st          pointer to a cs_xdef_t structure
 * \param[in]      st_mask     id related to this source term
 * \param[in, out] cell_mask   mask related to each cell to be updated
 */
/*----------------------------------------------------------------------------*/

static void
_set_mask(const cs_xdef_t     *st,
          int                  st_id,
          cs_mask_t           *cell_mask)
{
  if (st == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  const cs_mask_t  mask = (1 << st_id); // value of the mask for the source term

  if (st->meta & CS_FLAG_FULL_LOC) // All cells are selected
#   pragma omp parallel for if (cs_cdo_quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < cs_cdo_quant->n_cells; i++) cell_mask[i] |= mask;

  else {

    /* Retrieve information from the volume zone structure */
    const cs_volume_zone_t *z = cs_volume_zone_by_id(st->z_id);
    for (cs_lnum_t i = 0; i < z->n_cells; i++)
      cell_mask[z->cell_ids[i]] |= mask;

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
 * \param[in]       xv4        third vertex
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
 * \param[in]       anai     pointer to an analytical definition
 * \param[in]       cbf      pointer to a structure for face basis functions
 * \param[in]       xv1      first vertex
 * \param[in]       xv2      second vertex
 * \param[in]       xv3      third vertex
 * \param[in]       xv4      third vertex
 * \param[in]       vol      volume of the tetrahedron
 * \param[in, out]  cb       pointer to a cs_cell_builder_structure_t
 * \param[in, out]  array    array storing values to compute
 */
/*----------------------------------------------------------------------------*/

static void
_hho_add_tetra_by_ana(const cs_xdef_analytic_input_t  *anai,
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
  cs_real_t  *ana_eval = cb->values + 15;
  cs_real_t  *phi_eval = cb->values + 30;

  /* Compute Gauss points and related weights */
  cs_quadrature_tet_15pts(xv1, xv2, xv3, xv4, vol, gpts, gw);

  /* Evaluate the analytical function at the Gauss points */
  anai->func(cs_glob_time_step->t_cur, 15, NULL, (const cs_real_t *)gpts, true,
             anai->input, ana_eval);

  for (short int gp = 0; gp < 15; gp++) {

    cbf->eval_all_at_point(cbf, gpts[gp], phi_eval);

    const cs_real_t  w = gw[gp] * ana_eval[gp];
    for (short int i = 0; i < cbf->size; i++)
      array[i] += w * phi_eval[i];

  }  /* End of loop on Gauss points */
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      time_step  pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                   const cs_cdo_connect_t       *connect,
                                   const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */
  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
  cs_time_step = time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default flag related to a source term according to the
 *         numerical scheme chosen for discretizing an equation
 *
 * \param[in]       scheme    numerical scheme used for the discretization
 *
 * \return a default flag
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_source_term_set_default_flag(cs_param_space_scheme_t   scheme)
{
  cs_flag_t  meta_flag = 0;

  switch (scheme) {
  case CS_SPACE_SCHEME_CDOVB:
    meta_flag = CS_FLAG_DUAL | CS_FLAG_CELL; // Default
    break;

  case CS_SPACE_SCHEME_CDOFB:
    meta_flag = CS_FLAG_PRIMAL | CS_FLAG_CELL; // Default
    break;

  case CS_SPACE_SCHEME_CDOVCB:
  case CS_SPACE_SCHEME_HHO_P0:
  case CS_SPACE_SCHEME_HHO_P1:
  case CS_SPACE_SCHEME_HHO_P2:
    meta_flag = CS_FLAG_PRIMAL;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid numerical scheme to set a source term."));

  }

  return meta_flag;
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
  if (st == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  if (st->meta & flag)
    return; // Nothing to do

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
                " Stop modifying the source term flag.\n"
                " This case is not handled.");
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
  if (st == NULL)
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

cs_flag_t
cs_source_term_init(cs_param_space_scheme_t       space_scheme,
                    const int                     n_source_terms,
                    const cs_xdef_t             **source_terms,
                    cs_source_term_cellwise_t    *compute_source[],
                    cs_flag_t                    *sys_flag,
                    cs_mask_t                    *source_mask[])
{
  if (n_source_terms > CS_N_MAX_SOURCE_TERMS)
    bft_error(__FILE__, __LINE__, 0,
              " Limitation to %d source terms has been reached!",
              CS_N_MAX_SOURCE_TERMS);

  cs_flag_t  msh_flag = 0;
  *source_mask = NULL;
  for (short int i = 0; i < CS_N_MAX_SOURCE_TERMS; i++)
    compute_source[i] = NULL;

  if (n_source_terms == 0)
    return msh_flag;

  bool  need_mask = false;

  for (int st_id = 0; st_id < n_source_terms; st_id++) {

    const cs_xdef_t  *st_def = source_terms[st_id];

    if (st_def->meta & CS_FLAG_PRIMAL) {
      if (space_scheme == CS_SPACE_SCHEME_CDOVB ||
          space_scheme == CS_SPACE_SCHEME_CDOVCB) {
        msh_flag |= CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ |
          CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
        *sys_flag |= CS_FLAG_SYS_HLOC_CONF | CS_FLAG_SYS_SOURCES_HLOC;
      }
    }

    if ((st_def->meta & CS_FLAG_FULL_LOC) == 0) // Not defined on the whole mesh
      need_mask = true;

    switch (space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:

      msh_flag |= CS_CDO_LOCAL_PV;
      if (st_def->meta & CS_FLAG_DUAL) {

        switch (st_def->type) {

        case CS_XDEF_BY_VALUE:
          msh_flag |= CS_CDO_LOCAL_PVQ;
          compute_source[st_id] = cs_source_term_dcsd_by_value;
          break;

        case CS_XDEF_BY_ANALYTIC_FUNCTION:

          switch (st_def->qtype) {

          case CS_QUADRATURE_BARY:
            msh_flag |= CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_PFQ |
              CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FE  | CS_CDO_LOCAL_FEQ;
            compute_source[st_id] = cs_source_term_dcsd_bary_by_analytic;
            break;

          case CS_QUADRATURE_BARY_SUBDIV:
            msh_flag |= CS_CDO_LOCAL_EV | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_HFQ |
              CS_CDO_LOCAL_FE | CS_CDO_LOCAL_FEQ;
            compute_source[st_id] = cs_source_term_dcsd_q1o1_by_analytic;
            break;

          case CS_QUADRATURE_HIGHER:
            msh_flag |= CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FE |
              CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_PVQ |
              CS_CDO_LOCAL_PEQ;
            compute_source[st_id] = cs_source_term_dcsd_q10o2_by_analytic;
            break;

          case CS_QUADRATURE_HIGHEST:
            msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE |
              CS_CDO_LOCAL_EV;
            compute_source[st_id] = cs_source_term_dcsd_q5o3_by_analytic;
            break;

          default:
            bft_error(__FILE__, __LINE__, 0,
                      " Invalid type of quadrature for computing a source term"
                      " with CDOVB schemes");
          } // quad_type
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    " Invalid type of definition for a source term in CDOVB");
          break;
        } // switch def_type

      }
      else {
        assert(st_def->meta & CS_FLAG_PRIMAL);

        switch (st_def->type) {

        case CS_XDEF_BY_VALUE:
          msh_flag |= CS_CDO_LOCAL_PV;
          compute_source[st_id] = cs_source_term_pvsp_by_value;
          break;

        case CS_XDEF_BY_ANALYTIC_FUNCTION:
          msh_flag |= CS_CDO_LOCAL_PV;
          compute_source[st_id] = cs_source_term_pvsp_by_analytic;
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    " Invalid type of definition for a source term in CDOVB");
          break;

        } // switch def_type

      } // flag PRIMAL or DUAL
      break; // CDOVB

    case CS_SPACE_SCHEME_CDOVCB:
      if (st_def->meta & CS_FLAG_DUAL) {

        bft_error(__FILE__, __LINE__, 0,
                  " Invalid type of definition for a source term in CDOVB");

        /* TODO
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
          msh_flag |= CS_CDO_LOCAL_PV;
          compute_source[st_id] = cs_source_term_vcsp_by_value;
          break;

        case CS_XDEF_BY_ANALYTIC_FUNCTION:
          msh_flag |= CS_CDO_LOCAL_PV;
          compute_source[st_id] = cs_source_term_vcsp_by_analytic;
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    " Invalid type of definition for a source term in CDOVB");
          break;

        } // switch def_type

      }
      break; // CDOVCB

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
      switch (st_def->type) {

      case CS_XDEF_BY_VALUE:
        if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
          compute_source[st_id] = cs_source_term_fbvd_by_value;
        else
          compute_source[st_id] = cs_source_term_fbsd_by_value;
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        msh_flag |= CS_CDO_LOCAL_PV;
        if ((*sys_flag) & CS_FLAG_SYS_VECTOR)
          compute_source[st_id] = cs_source_term_fbvd_bary_by_analytic;
        else
          compute_source[st_id] = cs_source_term_fbsd_bary_by_analytic;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid type of definition for a source term in CDOVB");
        break;

      } // switch def_type
      break;

    case CS_SPACE_SCHEME_HHO_P1:
    case CS_SPACE_SCHEME_HHO_P2:
      switch (st_def->type) {

      case CS_XDEF_BY_VALUE:
        compute_source[st_id] = cs_source_term_hhosd_by_value;
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        compute_source[st_id] = cs_source_term_hhosd_by_analytic;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid type of definition for a source term in CDOVB");
        break;

      } // switch def_type
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "Invalid space scheme for setting the source term.");
      break;

    } // Switch on space scheme

  } // Loop on source terms

  if (need_mask) {

    const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;

    /* Initialize mask buffer */
    cs_mask_t  *mask = NULL;
    BFT_MALLOC(mask, n_cells, cs_mask_t);
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (int i = 0; i < n_cells; i++) mask[i] = 0;

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
 * \param[in, out] input           pointer to an element cast on-the-fly
 * \param[in, out] cb              pointer to a cs_cell_builder_t structure
 * \param[in, out] csys            cellwise algebraic system
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute_cellwise(const int                    n_source_terms,
                                const cs_xdef_t            **source_terms,
                                const cs_cell_mesh_t        *cm,
                                const cs_mask_t             *source_mask,
                                cs_source_term_cellwise_t   *compute_source[],
                                void                        *input,
                                cs_cell_builder_t           *cb,
                                cs_cell_sys_t               *csys)
{
  /* Reset local contributions */
  memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

  if (n_source_terms < 1)
    return;

  if (source_mask == NULL) { // All source terms are defined on the whole mesh

    for (short int st_id = 0; st_id < n_source_terms; st_id++) {

      cs_source_term_cellwise_t  *compute = compute_source[st_id];

      /* Contrib is updated inside */
      compute(source_terms[st_id], cm, cb, input, csys->source);

    } // Loop on source terms

  }
  else { /* Some source terms are only defined on a selection of cells */

    for (short int st_id = 0; st_id < n_source_terms; st_id++) {

      const cs_mask_t  st_mask = (1 << st_id);
      if (source_mask[cm->c_id] & st_mask) {

        cs_source_term_cellwise_t  *compute = compute_source[st_id];

        /* Contrib is updated inside */
        compute(source_terms[st_id], cm, cb, input, csys->source);

      } // Compute the source term on this cell

    } // Loop on source terms

  } // Source terms are defined on the whole domain or not ?

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a source term in the case of
 *         an input data which is a density
 *
 * \param[in]      loc        describe where is located the associated DoF
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in, out] p_values   pointer to the computed values (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute_from_density(cs_flag_t                loc,
                                    const cs_xdef_t         *source,
                                    double                  *p_values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  double  *values = *p_values;

  if (source == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  int  stride = 0;
  if (loc & CS_FLAG_SCALAR)
    stride = 1;
  else if (loc & CS_FLAG_VECTOR)
    stride = 3;
  else
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case\n", __func__);

  cs_lnum_t n_ent = 0;
  if (cs_flag_test(loc, cs_flag_dual_cell) ||
      cs_flag_test(loc, cs_flag_primal_vtx))
    n_ent = quant->n_vertices;
  else if (cs_flag_test(loc, cs_flag_primal_cell))
    n_ent = quant->n_cells;
  else
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case\n", __func__);

  /* Initialize values */
  if (values == NULL)
    BFT_MALLOC(values, n_ent*stride, double);
  memset(values, 0, n_ent*stride*sizeof(double));

  switch (source->type) {

  case CS_XDEF_BY_VALUE:
    cs_evaluate_density_by_value(loc, source, values);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_evaluate_density_by_analytic(loc, source, values);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _(" Invalid type of definition.\n"));

  } /* Switch according to def_type */

  /* Return values */
  *p_values = values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a source term in the case of
 *         an input data which is a potential
 *
 * \param[in]      loc        describe where is located the associated DoF
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in, out] p_values   pointer to the computed values (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute_from_potential(cs_flag_t                loc,
                                      const cs_xdef_t         *source,
                                      double                  *p_values[])
{
  const int  stride = 1; // Only this case is managed up to now
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  double  *values = *p_values;

  if (source == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_st));

  cs_lnum_t n_ent = 0;
  if (cs_flag_test(loc, cs_flag_dual_cell) ||
      cs_flag_test(loc, cs_flag_primal_vtx))
    n_ent = quant->n_vertices;
  else if (cs_flag_test(loc, cs_flag_primal_cell))
    n_ent = quant->n_cells;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid case. Not able to compute the source term.\n"));

  /* Initialize values */
  if (values == NULL)
    BFT_MALLOC(values, n_ent*stride, double);
  for (cs_lnum_t i = 0; i < n_ent*stride; i++)
    values[i] = 0.0;

  switch (source->type) {

  case CS_XDEF_BY_VALUE:
    cs_evaluate_potential_by_value(loc, source, values);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_evaluate_potential_by_analytic(loc, source, values);
    break;

    default:
      bft_error(__FILE__, __LINE__, 0, _(" Invalid type of definition.\n"));

  } /* Switch according to def_type */

  /* Return values */
  *p_values = values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar potential defined at primal vertices by a constant
 *         value.
 *         A discrete Hodge operator has to be computed before this call and
 *         stored inside a cs_cell_builder_t structure
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_pvsp_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values)
{
  if (source == NULL)
    return;

  CS_UNUSED(input);
  /* Sanity checks */
  assert(values != NULL && cm != NULL);
  assert(cb != NULL && cb->hdg != NULL);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PV));

  const cs_real_t *s_input = (const cs_real_t *)source->input;
  const cs_real_t  pot_value = s_input[0];

  /* Retrieve the values of the potential at each cell vertices */
  double  *eval = cb->values;
  for (short int v = 0; v < cm->n_vc; v++)
    eval[v] = pot_value;

  /* Multiply these values by a cellwise Hodge operator previously computed */
  double  *hdg_eval = cb->values + cm->n_vc;
  cs_sdm_square_matvec(cb->hdg, eval, hdg_eval);

  for (short int v = 0; v < cm->n_vc; v++)
    values[v] += hdg_eval[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar potential defined at primal vertices by an
 *         analytical function.
 *         A discrete Hodge operator has to be computed before this call and
 *         stored inside a cs_cell_builder_t structure
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_pvsp_by_analytic(const cs_xdef_t           *source,
                                const cs_cell_mesh_t      *cm,
                                cs_cell_builder_t         *cb,
                                void                      *input,
                                double                    *values)
{
  if (source == NULL)
    return;

  CS_UNUSED(input);
  /* Sanity checks */
  assert(values != NULL && cm != NULL);
  assert(cb != NULL && cb->hdg != NULL);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PV));

  const double  tcur = cs_time_step->t_cur;

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)source->input;

  /* Retrieve the values of the potential at each cell vertices */
  double  *eval = cb->values;
  anai->func(tcur, cm->n_vc, NULL, cm->xv,
             true, // compacted output ?
             anai->input,
             eval);

  /* Multiply these values by a cellwise Hodge operator previously computed */
  double  *hdg_eval = cb->values + cm->n_vc;
  cs_sdm_square_matvec(cb->hdg, eval, hdg_eval);

  for (short int v = 0; v < cm->n_vc; v++)
    values[v] += hdg_eval[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by a value.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values)
{
  CS_UNUSED(cb);
  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PVQ));

  const cs_real_t *s_input = (const cs_real_t *)source->input;
  const cs_real_t  density_value = s_input[0];

  for (int v = 0; v < cm->n_vc; v++)
    values[v] += density_value * cm->wvc[v] * cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_bary_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values)
{
  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_HFQ |
                      CS_CDO_LOCAL_FE  | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV));

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)source->input;

  /* Compute the barycenter of each portion of dual cells */
  cs_real_3_t  *xgv = cb->vectors;
  for (short int v = 0; v < cm->n_vc; v++)
    xgv[v][0] = xgv[v][1] = xgv[v][2] = 0.;

  for (short int f = 0; f < cm->n_fc; f++) {

    cs_real_3_t  xfc;

    const double  *xf = cm->face[f].center;
    const double  hf_coef = cs_math_onesix * cm->hfc[f];

    for (int k = 0; k < 3; k++) xfc[k] = 0.25*(xf[k] + cm->xc[k]);

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];
      const double  *xv1 = cm->xv + 3*v1, *xv2 = cm->xv + 3*v2;
      const double  tet_vol = cm->tef[i]*hf_coef;

      // xg = 0.25(xv1 + xe + xf + xc) where xe = 0.5*(xv1 + xv2)
      for (int k = 0; k < 3; k++)
        xgv[v1][k] += tet_vol*(xfc[k] + 0.375*xv1[k] + 0.125*xv2[k]);

      // xg = 0.25(xv2 + xe + xf + xc) where xe = 0.5*(xv1 + xv2)
      for (int k = 0; k < 3; k++)
        xgv[v2][k] += tet_vol*(xfc[k] + 0.375*xv2[k] + 0.125*xv1[k]);

    } // Loop on face edges

  } // Loop on cell faces

  /* Compute the source term contribution for each vertex */
  double  *vol_vc = cb->values;
  for (short int v = 0; v < cm->n_vc; v++) {
    vol_vc[v] = cm->vol_c * cm->wvc[v];
    const double  invvol = 1/vol_vc[v];
    for (int k = 0; k < 3; k++) xgv[v][k] *= invvol;
  }

  /* Call the analytic function to evaluate the function at xgv */
  const double  tcur = cs_time_step->t_cur;
  double  *eval_xgv = vol_vc + cm->n_vc;
  anai->func(tcur, cm->n_vc, NULL, (const cs_real_t *)xgv,
             true, // compacted output ?
             anai->input,
             eval_xgv);

  for (short int v = 0; v < cm->n_vc; v++)
    values[v] = cm->vol_c * cm->wvc[v] * eval_xgv[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use a the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q1o1_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values)
{
  CS_UNUSED(cb);
  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FE |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV));

  const double  tcur = cs_time_step->t_cur;

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)source->input;

  for (short int f = 0; f < cm->n_fc; f++) {

    cs_real_3_t  xg[2], xfc;
    cs_real_t  eval_xg[2];

    const double  *xf = cm->face[f].center;
    const double  hf_coef = cs_math_onesix * cm->hfc[f];

    for (int k = 0; k < 3; k++) xfc[k] = 0.25*(xf[k] + cm->xc[k]);

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];
      const double  *xv1 = cm->xv + 3*v1, *xv2 = cm->xv + 3*v2;
      const double  half_pef_vol = cm->tef[i]*hf_coef;

      // xg = 0.25(xv1 + xe + xf + xc) where xe = 0.5*(xv1 + xv2)
      for (int k = 0; k < 3; k++)
        xg[0][k] = xfc[k] + 0.375*xv1[k] + 0.125*xv2[k];

      // xg = 0.25(xv1 + xe + xf + xc) where xe = 0.5*(xv1 + xv2)
      for (int k = 0; k < 3; k++)
        xg[1][k] = xfc[k] + 0.375*xv2[k] + 0.125*xv1[k];

      anai->func(tcur, 2, NULL, (const cs_real_t *)xg,
                 true, // compacted output ?
                 anai->input,
                 eval_xg);

      values[v1] += half_pef_vol * eval_xg[0];
      values[v2] += half_pef_vol * eval_xg[1];

    } // Loop on face edges

  } // Loop on cell faces

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use a ten-point quadrature rule to evaluate the integral.
 *         Exact for quadratic function.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q10o2_by_analytic(const cs_xdef_t           *source,
                                      const cs_cell_mesh_t      *cm,
                                      cs_cell_builder_t         *cb,
                                      void                      *input,
                                      double                    *values)
{
  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);
  assert(cb != NULL);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FE  |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_PVQ |
                      CS_CDO_LOCAL_PEQ));

  const double  tcur = cs_time_step->t_cur;

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)source->input;

  /* Temporary buffers */
  double  *contrib = cb->values;              // size n_vc

  /* 1) Compute the contributions seen by the whole portion of dual cell */

  /* Cell evaluation */
  double  eval_c;
  anai->func(tcur, 1, NULL, cm->xc,
             true, // compacted output ?
             anai->input,
             &eval_c);

  /* Contributions related to vertices */
  double  *eval_v = cb->values + cm->n_vc; // size n_vc
  anai->func(tcur, cm->n_vc, NULL, cm->xv,
             true,  // compacted output ?
             anai->input,
             eval_v);

  cs_real_3_t  *xvc = cb->vectors;
  for (short int v = 0; v < cm->n_vc; v++) {
    const double  *xv = cm->xv + 3*v;
    for (int k = 0; k < 3; k++) xvc[v][k] = 0.5*(cm->xc[k] + xv[k]);
  }

  double  *eval_vc = cb->values + 2*cm->n_vc; // size n_vc
  anai->func(tcur, cm->n_vc, NULL, (const cs_real_t *)xvc,
             true,  // compacted output ?
             anai->input,
             eval_vc);

  for (short int v = 0; v < cm->n_vc; v++) {

    /* Set the initial values
       -1/20 on extremity points and 1/5 on midpoints */
      contrib[v] = cm->wvc[v]*cm->vol_c
      * (-0.05*(eval_c + eval_v[v]) + 0.2*eval_vc[v]);

  } // Loop on vertices

  /* 2) Compute the contribution related to edge
     The portion of dual cell seen by each vertex is 1/2 |pec| */
  cs_real_3_t  *x_e = cb->vectors;
  cs_real_3_t  *xec = cb->vectors + cm->n_ec; // size = n_ec (overwrite xvc)

  for (short int e = 0; e < cm->n_ec; e++) {
    for (int k = 0; k < 3; k++) {
      x_e[e][k] = cm->edge[e].center[k];
      xec[e][k] = 0.5*(cm->xc[k] + x_e[e][k]);
    }
  }

  // Evaluate the analytic function at xe and xec
  double  *eval_e = cb->values + cm->n_vc; // size=n_ec (overwrite eval_v)
  double  *eval_ec = eval_e + cm->n_ec;    // size=n_ec (overwrite eval_vc)
  anai->func(tcur, 2*cm->n_ec, NULL, (const cs_real_t *)cb->vectors,
             true,  // compacted output ?
             anai->input,
             eval_e);

  // xev (size = 2*n_ec)
  cs_real_3_t  *xve = cb->vectors;         // size=2*n_ec (overwrite xe and xec)
  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_real_t  *xe = cm->edge[e].center;
    const short int  v1 = cm->e2v_ids[2*e];
    const double  *xv1 = cm->xv + 3*v1;
    const short int  v2 = cm->e2v_ids[2*e+1];
    const double  *xv2 = cm->xv + 3*v2;

    for (int k = 0; k < 3; k++) {
      xve[2*e  ][k] = 0.5*(xv1[k] + xe[k]);
      xve[2*e+1][k] = 0.5*(xv2[k] + xe[k]);
    }

  } // Loop on edges

  double  *eval_ve = eval_ec + cm->n_ec; // size = 2*n_ec
  anai->func(tcur, 2*cm->n_ec, NULL, (const cs_real_t *)cb->vectors,
             true,  // compacted output ?
             anai->input,
             eval_ve);

  /* 3) Main loop on faces */
  double  *pvf_vol = eval_ve + 2*cm->n_ec;  // size n_vc

  for (short int f = 0; f < cm->n_fc; f++) {

    const double  *xf = cm->face[f].center;
    const double  hfc = cm->hfc[f];

    /* Reset volume of the face related to a vertex */
    for (short int v = 0; v < cm->n_vc; v++) pvf_vol[v] = 0;

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];
      const double  half_pef_vol = cs_math_onesix * cm->tef[i] * hfc;

      pvf_vol[v1] += half_pef_vol;
      pvf_vol[v2] += half_pef_vol;

      cs_real_3_t  xef;
      cs_real_t  eval_ef;
      for (int k = 0; k < 3; k++) xef[k] = 0.5*(cm->edge[e].center[k] + xf[k]);
      anai->func(tcur, 1, NULL, xef,
                 true,  // compacted output ?
                 anai->input,
                 &eval_ef);

      // 1/5 (EF + EC) -1/20 * (E)
      const double  common_ef_contrib =
        0.2*(eval_ef + eval_ec[e]) -0.05*eval_e[e];

      contrib[v1] += half_pef_vol*(common_ef_contrib + 0.2*eval_ve[2*e]);
      contrib[v2] += half_pef_vol*(common_ef_contrib + 0.2*eval_ve[2*e+1]);

    } // Loop on face edges

    /* Contributions related to this face */
    cs_real_3_t  *xvfc = cb->vectors;  // size=2+n_vc (overwrite xev)
    for (int k = 0; k < 3; k++) {
      xvfc[0][k] = xf[k];                    // xf
      xvfc[1][k] = 0.5*(xf[k] + cm->xc[k]);  // xfc
    }

    short int  n_vf = 0;
    for (short int v = 0; v < cm->n_vc; v++) {
      if (pvf_vol[v] > 0) {
        cb->ids[n_vf] = v;
        for (int k = 0; k < 3; k++)
          xvfc[2+n_vf][k] = 0.5*(xf[k] + cm->xv[3*v+k]);
        n_vf++;
      }
    }

    double  *eval_vfc = pvf_vol + cm->n_vc; // size=n_vf + 2
    anai->func(tcur, 2+n_vf, NULL, (const cs_real_t *)xvfc,
               true,  // compacted output ?
               anai->input,
               eval_vfc);

    for (short int i = 0; i < n_vf; i++) {
      short int  v = cb->ids[i];
      const double  val_vfc = -0.05*eval_vfc[0] + 0.2*eval_vfc[1];
      contrib[v] += pvf_vol[v] * (val_vfc + 0.2*eval_vfc[2+i]);
    }

  } // Loop on cell faces

  /* Add the computed contributions to the return values */
  for (short int v = 0; v < cm->n_vc; v++)
    values[v] += contrib[v];

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use a five-point quadrature rule to evaluate the integral.
 *         Exact for cubic function.
 *         This function may be expensive since many evaluations are needed.
 *         Please use it with care.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q5o3_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values)
{
  double  sum, weights[5], results[5];
  cs_real_3_t  gauss_pts[5];

  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);
  assert(cb != NULL);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE |
                      CS_CDO_LOCAL_EV));

  const double  tcur = cs_time_step->t_cur;
  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)source->input;

  /* Temporary buffers */
  double  *contrib = cb->values;
  for (short int v = 0; v < cm->n_vc; v++) contrib[v] = 0;

  /* Main loop on faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const double  *xf = cm->face[f].center;

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];
      const double  tet_vol = 0.5*cs_math_voltet(cm->xv + 3*v1,
                                                 cm->xv + 3*v2,
                                                 xf,
                                                 cm->xc);

      /* Compute Gauss points and its weights */
      cs_quadrature_tet_5pts(cm->xv + 3*v1, cm->edge[e].center, xf, cm->xc,
                             tet_vol,
                             gauss_pts, weights);

      anai->func(tcur, 5, NULL, (const cs_real_t *)gauss_pts,
                 true,  // compacted output ?
                 anai->input,
                 results);

      sum = 0.;
      for (int p = 0; p < 5; p++) sum += results[p] * weights[p];
      contrib[v1] += sum;

      /* Compute Gauss points and its weights */
      cs_quadrature_tet_5pts(cm->xv + 3*v2, cm->edge[e].center, xf, cm->xc,
                             tet_vol,
                             gauss_pts, weights);

      anai->func(tcur, 5, NULL, (const cs_real_t *)gauss_pts,
                 true,  // compacted output ?
                 anai->input,
                 results);

      sum = 0.;
      for (int p = 0; p < 5; p++) sum += results[p] * weights[p];
      contrib[v2] += sum;

    } // Loop on face edges

  } // Loop on cell faces

  /* Add the computed contributions to the return values */
  for (short int v = 0; v < cm->n_vc; v++)
    values[v] += contrib[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar potential defined at primal vertices and cells
 *         by a constant value.
 *         A discrete Hodge operator has to be computed before this call and
 *         stored inside a cs_cell_builder_t structure
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_vcsp_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values)
{
  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);
  assert(cb != NULL && cb->hdg != NULL);

  const cs_real_t *s_input = (const cs_real_t *)source->input;
  const cs_real_t  pot_value = s_input[0];

  /* Retrieve the values of the potential at each cell vertices */
  double  *eval = cb->values;
  for (short int v = 0; v < cm->n_vc; v++)
    eval[v] = pot_value;
  eval[cm->n_vc] = pot_value;

  /* Multiply these values by a cellwise Hodge operator previously computed */
  double  *hdg_eval = cb->values + cm->n_vc + 1;
  cs_sdm_square_matvec(cb->hdg, eval, hdg_eval);

  for (short int v = 0; v < cm->n_vc + 1; v++)
    values[v] += hdg_eval[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar potential defined at primal vertices and cells by
 *         an analytical function.
 *         A discrete Hodge operator has to be computed before this call and
 *         stored inside a cs_cell_builder_t structure
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_vcsp_by_analytic(const cs_xdef_t           *source,
                                const cs_cell_mesh_t      *cm,
                                cs_cell_builder_t         *cb,
                                void                      *input,
                                double                    *values)
{
  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);
  assert(cb != NULL && cb->hdg != NULL);

  const double  tcur = cs_time_step->t_cur;
  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)source->input;

  /* Retrieve the values of the potential at each cell vertices */
  double  *eval = cb->values;

  anai->func(tcur, cm->n_vc, NULL, cm->xv,
             true,  // compacted output ?
             anai->input,
             eval);

  anai->func(tcur, 1, NULL, cm->xc,
             true,  // compacted output ?
             anai->input,
             eval + cm->n_vc);

  /* Multiply these values by a cellwise Hodge operator previously computed */
  double  *hdg_eval = cb->values + cm->n_vc + 1;
  cs_sdm_square_matvec(cb->hdg, eval, hdg_eval);

  for (short int v = 0; v < cm->n_vc + 1; v++)
    values[v] += hdg_eval[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density (sd) defined on primal cells by a value.
 *         Case of face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fbsd_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values)
{
  CS_UNUSED(cb);
  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);

  const cs_real_t *s_input = (const cs_real_t *)source->input;

  values[cm->n_fc] = s_input[0] * cm->vol_c;
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
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fbvd_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values)
{
  CS_UNUSED(cb);
  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);

  const cs_real_t *v_input = (const cs_real_t *)source->input;
  for (int k = 0; k < source->dim; k++)
    values[source->dim*cm->n_fc + k] = v_input[k] * cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at primal cells by an analytical
 *         function.
 *         Use the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *         Case of face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fbsd_bary_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values)
{
  CS_UNUSED(cb);
  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)source->input;

  /* Call the analytic function to evaluate the function at xc */
  const double  tcur = cs_time_step->t_cur;
  double  eval_xc;

  anai->func(tcur, 1, NULL, (const cs_real_t *)cm->xc,
             true,  // compacted output ?
             anai->input,
             &eval_xc);

  values[cm->n_fc] = cm->vol_c * eval_xc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a vector-valued density defined at primal cells by an
 *         analytical function.
 *         Use the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *         Case of face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fbvd_bary_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values)
{
  CS_UNUSED(cb);
  CS_UNUSED(input);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL);
  assert(source->dim < 4);

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)source->input;

  /* Call the analytic function to evaluate the function at xc */
  const double  tcur = cs_time_step->t_cur;
  double  eval_xc[3];

  anai->func(tcur, 1, NULL, (const cs_real_t *)cm->xc,
             true,  // compacted output ?
             anai->input,
             eval_xc);

  for (int k = 0; k < source->dim; k++)
    values[source->dim*cm->n_fc + k] = cm->vol_c * eval_xc[k];
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
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_hhosd_by_value(const cs_xdef_t           *source,
                              const cs_cell_mesh_t      *cm,
                              cs_cell_builder_t         *cb,
                              void                      *input,
                              double                    *values)
{
  CS_UNUSED(cb);

  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL && input != NULL);

  const cs_real_t  *const_val = (const cs_real_t *)source->input;

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
        for (short int f = 0; f < cm->n_fc; ++f) {

          const cs_quant_t  pfq = cm->face[f];
          const double  hf_coef = cs_math_onethird * cm->hfc[f];
          const int  start = cm->f2e_idx[f];
          const int  end = cm->f2e_idx[f+1];
          const short int n_vf = end - start; // #vertices (=#edges)
          const short int *f2e_ids = cm->f2e_ids + start;

          assert(n_vf > 2);
          switch(n_vf){

          case 3: /* triangle (optimized version, no subdivision) */
            {
              short int  v0, v1, v2;
              cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids,
                                               &v0, &v1, &v2);

              _hho_add_tetra_by_val(const_val[0], cbf,
                                    cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2,
                                    cm->xc,
                                    hf_coef * pfq.meas, cb, cell_values);

            }
            break;

          default:
            {
              const double  *tef = cm->tef + start;

              for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

                // Edge-related variables
                const short int e0  = f2e_ids[e];
                const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
                const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

                _hho_add_tetra_by_val(const_val[0], cbf,
                                      xv0, xv1, pfq.center, cm->xc,
                                      hf_coef*tef[e], cb, cell_values);
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
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_hhosd_by_analytic(const cs_xdef_t           *source,
                                 const cs_cell_mesh_t      *cm,
                                 cs_cell_builder_t         *cb,
                                 void                      *input,
                                 double                    *values)
{
  CS_UNUSED(cb);
  if (source == NULL)
    return;

  /* Sanity checks */
  assert(values != NULL && cm != NULL && input != NULL);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV));

  cs_hho_builder_t  *hhob = (cs_hho_builder_t *)input;
  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)source->input;
  cs_real_t  *cell_values = values + cm->n_fc * hhob->face_basis[0]->size;

  const cs_basis_func_t  *cbf = hhob->cell_basis;

  /* Reset cell values */
  memset(cell_values, 0, sizeof(cs_real_t)*cbf->size);

  /* Switch according to the cell type: optimised version for tetra */
  switch (cm->type) {

  case FVM_CELL_TETRA:
    {
      assert(cm->n_fc == 4 && cm->n_vc == 4);
      _hho_add_tetra_by_ana(anai, cbf,
                            cm->xv, cm->xv+3, cm->xv+6, cm->xv+9,
                            cm->vol_c, cb, cell_values);
    }
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:
  {
    for (short int f = 0; f < cm->n_fc; ++f) {

      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = cs_math_onethird * cm->hfc[f];
      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_vf = end - start; // #vertices (=#edges)
      const short int *f2e_ids = cm->f2e_ids + start;

      assert(n_vf > 2);
      switch(n_vf){

      case 3: /* triangle (optimized version, no subdivision) */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          _hho_add_tetra_by_ana(anai, cbf,
                                cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2, cm->xc,
                                hf_coef * pfq.meas, cb, cell_values);

        }
        break;

      default:
        {
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            // Edge-related variables
            const short int e0  = f2e_ids[e];
            const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
            const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

            _hho_add_tetra_by_ana(anai, cbf,
                                  xv0, xv1, pfq.center, cm->xc,
                                  hf_coef*tef[e], cb, cell_values);
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

END_C_DECLS
