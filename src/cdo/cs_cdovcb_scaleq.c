/*============================================================================
 * Build an algebraic CDO vertex+cell-based system for unsteady convection
 * diffusion reaction scalar equations with source terms
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include <float.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_scheme_geometry.h"
#include "cs_equation_common.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_search.h"
#include "cs_source_term.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdovcb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVCB_SCALEQ_DBG  0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/* Algebraic system for CDO vertex+cell-based discretization */

struct _cs_cdovcb_scaleq_t {

  /* Pointer to a cs_equation_param_t structure shared with a cs_equation_t
     structure.  */

  const cs_equation_param_t  *eqp;

  /* System size */
  cs_lnum_t  n_vertices;
  cs_lnum_t  n_cells;
  cs_lnum_t  n_sys;        // n_vertices + n_cells
  int        n_max_vcbyc;  // n_max_vbyc + 1

  /* Shortcut to know what to build */
  bool       has[N_CDO_TERMS];
  cs_flag_t  flag;

  /* Store the matrix to invert after assembling and static condensation for
     upper left part
     Store in the right and left bottom part quantities attached to the system
     in order to recover the value at each cell centers */
  cs_sla_hmatrix_t      *hybrid_storage;

  /* Store the values of the field at cell centers */
  cs_real_t             *cell_values;

  /* Common members for all terms */
  double                *loc_vals; // local temporary values
  cs_hodge_builder_t    *hb;       // can be used by reaction, time or source

  /* Builder structure for diffusion term */
  bool                   diff_pty_uniform;
  cs_cdo_diff_t         *diff;

  /* Builder structure for advection term */
  cs_cdo_adv_t          *adv;

  /* Time term */
  bool                   time_pty_uniform;
  bool                   time_mat_is_diag;
  double                 time_pty_val;

  /* Reaction term */
  bool                  *reaction_pty_uniform;
  double                *reaction_pty_val;

  /* Source terms */
  cs_real_t             *source_terms;

  /* Boundary conditions:

     face_bc should not change during the simulation.
     The case of a definition of the BCs which changes of type during the
     simulation is possible but not implemented.
     You just have to call the initialization step each time the type of BCs
     is modified to define an updated cs_cdo_bc_t structure.

     We translate Dirichlet BCs to border vertices
     The values can be modified at each time step in case of transient
     simulation.
     For Neumann and Robin BCs, the treatment is more complex since the
     contributions of these BCs to a dual cell related to a border vertex is
     computed from the contribution of different portions of primal faces
     (these primal border faces form a closure of the dual cell).
     These contributions are computed on the fly.

   */

  cs_param_bc_enforce_t  enforce; // type of enforcement of BCs
  cs_cdo_bc_t           *face_bc; // list of faces sorted by type of BCs
  cs_cdo_bc_list_t      *vtx_dir; // list of vertices attached to a Dirichlet BC
  double                *dir_val; // size = vtx_dir->n_nhmg_elts

  /* In case of a weak enforcement of the Dirichlet BCs */
  cs_lnum_t             *c2bcbf_idx;  // size: n_cells + 1
  cs_lnum_t             *c2bcbf_ids;  // cell --> border faces ids

};

/*============================================================================
 * Private variables
 *============================================================================*/

static double  cs_cdovcb_threshold = 1e-12; // Set during initialization
static cs_locmat_t  *cs_cell_condmat = NULL;

/* Size = 1 if openMP is not used */
static cs_cdo_locsys_t  **cs_cdovcb_cell_systems = NULL;

/* Hodge^{VC,Conf} used for the computation of source terms
   Special matrix: hmatrix = hybrid matrix (cf. cs_sla.h)
   In all cases, this matrix along with its index are shared */
static cs_sla_hmatrix_t  *cs_cdovcb_hconf = NULL;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_time_step_t  *cs_shared_time_step;

/* Flag to indicate which members have to be built in a cs_cell_mesh_t
   structure */
static const cs_flag_t  cs_cdovcb_cmflag = CS_CDO_LOCAL_V | CS_CDO_LOCAL_E |
  CS_CDO_LOCAL_EV | CS_CDO_LOCAL_F | CS_CDO_LOCAL_FE;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contribution of source terms to the rhs for this
 *          time step
 *
 * \param[in, out] b          pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] full_rhs   right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_add_source_terms(cs_cdovcb_scaleq_t    *b,
                  cs_real_t              full_rhs[])
{
  const cs_equation_param_t  *eqp = b->eqp;

  if (!b->has[CDO_SOURCETERM]) // Test if there is at least one source term
    return;

  if (b->has[CDO_TIME]) {

    const cs_param_time_t  t_info = eqp->time_info;

    /* Previous values are stored inside b->source_terms i.e.
       values of the source terms related to t_prev */
    if (t_info.scheme == CS_TIME_SCHEME_EXPLICIT)
# pragma omp parallel for if (b->n_sys > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_sys; i++)
        full_rhs[i] += b->source_terms[i];

    else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
             t_info.scheme == CS_TIME_SCHEME_THETA) {

      const double  tcoef = 1 - t_info.theta;

# pragma omp parallel for if (b->n_sys > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_sys; i++)
        full_rhs[i] += tcoef * b->source_terms[i];

    }

    /* Update b->source_term with the value attached to t_cur */
    cs_cdovcb_scaleq_compute_source(b);

    if (t_info.scheme == CS_TIME_SCHEME_IMPLICIT)
# pragma omp parallel for if (b->n_sys > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_sys; i++)
        full_rhs[i] += b->source_terms[i];

    else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
             t_info.scheme == CS_TIME_SCHEME_THETA) {
# pragma omp parallel for if (b->n_sys > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_sys; i++)
        full_rhs[i] += t_info.theta * b->source_terms[i];

    }

  }
  else { /* Steady case: source terms have already been computed during
            the initialization step */

# pragma omp parallel for if (b->n_sys > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < b->n_sys; i++)
      full_rhs[i] += b->source_terms[i];

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define a discrete Hodge VC conforming op.
 *
 * \param[in, out] b      pointer to a cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_hvc_conf(cs_cdovcb_scaleq_t        *b)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Initialize a hybrid matrix structure */
  cs_cdovcb_hconf = cs_sla_hmatrix_create(b->n_vertices,
                                          b->n_cells,
                                          true,  // (0,1) (1,0) are transposed
                                          true,  // (0,0) is symmetric
                                          connect->v2v,
                                          connect->c2v);

  // To be modified for an fully integration of openMP
  cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(0);

  /* Cellwise construction ==> Loop on cells */
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Set the local mesh structure for the current cell */
    cs_cell_mesh_build(c_id, cs_cdovcb_cmflag, connect, quant, cm);

    /* Build the local dense matrix related to this operator */
    cs_locmat_t  *hloc = cs_hodge_build_cellwise(cm, b->hb);

    // TODO: Set this openMP section to critical
    /* Assemble the cellwise matrix into the "global" matrix */
    cs_sla_assemble_hmat_sym(hloc, cs_cdovcb_hconf);

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      field_val    pointer to the current value of the field
 * \param[in, out] builder      pointer to a cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_dir_values(const cs_mesh_t            *mesh,
                    const cs_real_t            *field_val,
                    const cs_cdovcb_scaleq_t   *b)
{
  const cs_cdo_bc_list_t  *vtx_dir = b->vtx_dir;
  const cs_equation_param_t  *eqp = b->eqp;

  if (vtx_dir->n_nhmg_elts == 0)
    return; // Nothing to do

  cs_flag_t  dof_flag = cs_cdo_primal_vtx | CS_FLAG_SCAL;

  /* Get the value of the Dirichlet for the current time */
  cs_cdo_bc_dirichlet_set(dof_flag,
                          cs_shared_time_step,
                          mesh,
                          eqp->bc,
                          vtx_dir,
                          b->dir_val);

  /* Update the values of dir_val in the following case */
  if (b->enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
      b->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
    if (b->has[CDO_TIME]) {

      const cs_param_time_t  t_info = eqp->time_info;

      /* Previous values of the unknown are stored inside field_val (iter n) */
      if (t_info.scheme == CS_TIME_SCHEME_EXPLICIT)
        for (cs_lnum_t i = 0; i < vtx_dir->n_nhmg_elts; i++)
          b->dir_val[i] = field_val[vtx_dir->elt_ids[i]];

      else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
               t_info.scheme == CS_TIME_SCHEME_THETA) {

        const double  tcoef = 1 - t_info.theta;

        for (cs_lnum_t i = 0; i < vtx_dir->n_nhmg_elts; i++) {
          b->dir_val[i] *= t_info.theta;
          b->dir_val[i] += tcoef * field_val[vtx_dir->elt_ids[i]];
        }

      }

    } /* Unsteady */

  } /* Enforcement is not strong or strongly penalized */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the time discretization to the local system
 *
 * \param[in, out] tpty_val   current value of the time property
 * \param[in]      vtx_val    pointer to the current value of the field
 * \param[in]      cm         pointer to a cs_locmesh_t structure
 * \param[in]      loc_hconf  pointer to a conforming discrete Hodge op.
 * \param[in, out] b          pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] loc_sys    pointer to a cs_locmat_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_apply_time_scheme(double                   tpty_val,
                   const cs_real_t         *vtx_val,
                   const cs_cell_mesh_t    *cm,
                   const cs_locmat_t       *loc_hconf,
                   cs_cdovcb_scaleq_t      *b,
                   cs_cdo_locsys_t         *loc_sys)
{
  const int  n_vc = cm->n_vc;
  const int  n_sys = n_vc + 1;
  const cs_param_time_t  t_info = b->eqp->time_info;
  const cs_real_t  *cell_val = b->cell_values;

  cs_locmat_t  *loc_mat = loc_sys->mat; // square matrix of size n_vc + 1
  double  *loc_rhs = loc_sys->rhs;      // vector of size n_vc + 1

  /* Temporary buffers of size equal to the number of cell vertices */
  double  *fval = b->loc_vals;
  double  *adr_pn = b->loc_vals + n_sys;

  /* Sanity check */
  assert(n_sys == loc_mat->n_ent);

  /* Set the values of the fields attached to this cell */
  for (short int v = 0; v < n_vc; v++)
    fval[v] = vtx_val[cm->v_ids[v]];
  fval[n_vc] = cell_val[cm->c_id];

  if (b->time_mat_is_diag) { // Lumped or Voronoi --> only diag.

    /* |c| * 0.75* wvc = 0.75 * |dual_cell(v) cap c| for cell vertices
       |c| * 0.25 for cell entry */
    tpty_val *= cm->vol_c;

    switch (t_info.scheme) {

    case CS_TIME_SCHEME_EXPLICIT:

      /* Compute (Adv + Dif + Rea)*p^(n) */
      cs_locmat_matvec(loc_mat, fval, adr_pn);

      /* Reset the local system matrix to assemble */
      for (short int vc = 0; vc < n_sys*n_sys; vc++)
        loc_mat->val[vc] = 0;

      /* Update the rhs with -(Adv + Dif + Rea)*p^(n) + TimeMat*p^(n)
         Set the local matrix to a (time) diagonal matrix */
      for (short int v = 0; v <  n_vc; v++) {

        const double dval_v = 0.75 * tpty_val * cm->wvc[v];

        loc_mat->val[v*n_sys + v] = dval_v;
        loc_rhs[v] += -adr_pn[v] + dval_v * fval[v];

      }

      /* Cell row and cell_rhs entries are updated */
      {
        const double  dval_c = 0.25*tpty_val;

        loc_mat->val[n_vc*n_sys + n_vc] = dval_c;
        loc_rhs[n_vc] += -adr_pn[n_vc] + dval_c * fval[n_vc];
      }
      break;

    case CS_TIME_SCHEME_IMPLICIT:

      /* Update the rhs with TimeMat*p^(n)
         Update the local matrix with the time diagonal matrix */
      for (short int v = 0; v < n_vc; v++) {

        const double dval_v = 0.75 * tpty_val * cm->wvc[v];

        loc_mat->val[v*n_sys + v] += dval_v;
        loc_rhs[v] += dval_v * fval[v];

      }

      /* Cell row and cell_rhs entries are updated */
      {
        const double  dval_c = 0.25*tpty_val;

        loc_mat->val[n_vc*n_sys + n_vc] += dval_c;
        loc_rhs[n_vc] += dval_c * fval[n_vc];
      }
      break;

    case CS_TIME_SCHEME_CRANKNICO:
    case CS_TIME_SCHEME_THETA:
      {
        const double  tcoef = t_info.theta - 1;

        /* Compute (Adv + Dif + Rea)*p^(n) */
        cs_locmat_matvec(loc_mat, fval, adr_pn);

        /* Update the rhs with (theta-1)(Adv + Dif + Rea)*p^(n) + TimeMat*p^(n)
           Update the matrix m -> theta*m + time contrib (only on diagonal) */
        for (short int vi = 0; vi < n_vc; vi++) {

          const double dval_v = 0.75 * tpty_val * cm->wvc[vi];

          loc_rhs[vi] += tcoef * adr_pn[vi] + dval_v * fval[vi];

          double  *mi = loc_mat->val + vi*n_sys;
          mi[n_vc] *= t_info.theta; // cell column
          for (short int vj = 0; vj < n_vc; vj++) {

            mi[vj] *= t_info.theta;
            if (vi == vj)
              mi[vj] += dval_v;

          } // Loop on cell vertices (vj)

        } // Loop on cell vertices (vi)


        /* Cell row and cell_rhs entries are updated */
        const double  dval_c = 0.25*tpty_val;

        double  *mc = loc_mat->val + n_vc*n_sys;

        loc_rhs[n_vc] += tcoef * adr_pn[n_vc] + dval_c * fval[n_vc];
        mc[n_vc] *= t_info.theta;
        mc[n_vc] += dval_c;
        for (short int vj = 0; vj < n_vc; vj++)
          mc[vj] *= t_info.theta;

      }
      break;

    default:
      break;

    } /* Switch on time scheme */

  } /* Time matrix is diagonal */
  else {

    /* Time matrix is not diagonal anymore.
       time matrix = tpty_val * Hconf */

    assert(loc_hconf != NULL);
    assert(loc_hconf->n_ent == loc_mat->n_ent);

    switch (t_info.scheme) {

    case CS_TIME_SCHEME_EXPLICIT:

      /* Compute (Adv + Dif + Rea)*p^(n)
         Update the rhs with -(Adv + Dif + Rea)*p^(n) */
      cs_locmat_matvec(loc_mat, fval, adr_pn);
      for (short int vc = 0; vc < n_sys; vc++)
        loc_rhs[vc] -= adr_pn[vc];

      /* Replace the local system matrix by that of time */
      for (short int vi = 0; vi < n_sys; vi++) {

        double  *mi = loc_mat->val + vi*n_sys;
        double  *hi = loc_hconf->val + vi*n_sys;

        for (short int vj = 0; vj < n_sys; vj++)
          mi[vj] = tpty_val * hi[vj];

      }

      /* Update the rhs with TimeMat*p^(n) */
      cs_locmat_matvec(loc_mat, fval, adr_pn);
      for (short int vc = 0; vc < n_sys; vc++)
        loc_rhs[vc] += adr_pn[vc];

      break;

    case CS_TIME_SCHEME_IMPLICIT:

      /* Update the rhs with TimeMat*p^(n) */
      cs_locmat_matvec(loc_hconf, fval, adr_pn);
      for (short int vc = 0; vc < n_sys; vc++)
        loc_rhs[vc] += tpty_val * adr_pn[vc];

      /* Update the local system with the time diagonal matrix */
      for (short int vi = 0; vi < n_sys; vi++) {

        double  *mi = loc_mat->val + vi*n_sys;
        double  *hi = loc_hconf->val + vi*n_sys;

        for (short int vj = 0; vj < n_sys; vj++)
          mi[vj] += tpty_val * hi[vj];

      }
      break;

    case CS_TIME_SCHEME_CRANKNICO:
    case CS_TIME_SCHEME_THETA:
      {
        const double  tcoef = t_info.theta - 1;

        /* Update full RHS with -(Adv + Dif + Rea)*p^(n) + TimeMat*p^(n) */
        /* Compute (Adv + Dif + Rea)*p^(n) */
        cs_locmat_matvec(loc_mat, fval, adr_pn);
        for (short int vc = 0; vc < n_sys; vc++)
          loc_rhs[vc] += tcoef * adr_pn[vc];

        /* Update the rhs with TimeMat*p^(n) */
        cs_locmat_matvec(loc_hconf, fval, adr_pn);
        for (short int vc = 0; vc < n_sys; vc++)
          loc_rhs[vc] += tpty_val * adr_pn[vc];

        /* Update the local system with the time matrix */
        for (short int vi = 0; vi < n_sys; vi++) {

          double  *mi = loc_mat->val + vi*n_sys;
          double  *hi = loc_hconf->val + vi*n_sys;

          for (short int vj = 0; vj < n_sys; vj++) {

            mi[vj] *= t_info.theta;
            mi[vj] += tpty_val * hi[vj];

          } // Loop on cell vertices (vj)

        } // Loop on cell vertices (vi)

      }
      break;

    default:
      break;

    } /* Switch on time scheme */

  } /* Time matrix is not diagonal */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Proceed to a static condensation of the local system and keep
 *          information inside hydrid_storage to be able to compute the values
 *          at cell centers
 *          rhs stored in cs_cell_sys is also updated.
 *
 * \param[in, out] builder     pointer to a cs_cdovcb_scaleq_t structure
 * \param[in, out] locsys      pointer to a cs_cdo_locsys_t structure
 * \param[in, out] condmat     pointer to the local condensed matrix
 */
/*----------------------------------------------------------------------------*/

static void
_condense_and_store(cs_cdovcb_scaleq_t        *b,
                    cs_cdo_locsys_t           *locsys,
                    cs_locmat_t               *condmat)
{
  const int  n_sys = locsys->mat->n_ent;
  const int  n_vc = n_sys - 1;
  const cs_lnum_t  c_id = locsys->mat->ids[n_vc];
  const double  *mc = locsys->mat->val + n_sys*n_vc; // Last row
  const double  inv_acc = 1/mc[n_vc];
  const double  sc_contrib = inv_acc * locsys->rhs[n_vc];

  cs_sla_hmatrix_t  *hmat = b->hybrid_storage;
  double  *cc_diag = hmat->cc_diag;
  double  *cx_vals = hmat->cx_vals + cs_shared_connect->c2v->idx[c_id];

  /* Define condmat, update rhs and store information inside hybrid_storage */
  condmat->n_ent = n_vc;
  cc_diag[c_id] = inv_acc;

  for (short int vi = 0; vi < n_vc; vi++) {

    condmat->ids[vi] = locsys->mat->ids[vi];

    double  *mi = locsys->mat->val + n_sys*vi;
    double  *condmi = condmat->val + n_vc*vi;

    /* Store c2x_vals */
    cx_vals[vi] = mc[vi];

    /* Update RHS: RHS_v = RHS_v - Avc*Acc^-1*s_c */
    locsys->rhs[vi] -= sc_contrib * mi[n_vc];

    /* Condensate the local matrix */
    for (short int vj = 0; vj < n_vc; vj++)
      condmi[vj] = mi[vj] - inv_acc*mc[vj]*mi[n_vc];

  } // Loop on vi cell vertices

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Update the value of stabilization coefficient in given situation
 *
 * \param[in] adv     pointer to a cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_set_cip_coef(const cs_adv_field_t     *adv)
{
  assert(adv != NULL); // Sanity check

  const double  gseed = 1e-2;  /* Default value to multiply according to the
                                  problem and the ratio of diameters */

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const double  hc_max = quant->cell_info.h_max;
  const double  hc_min = quant->cell_info.h_min;
  const double  hf_max = quant->face_info.h_max;
  const double  hf_min = quant->face_info.h_min;
  const double  hcMm = hc_max * hc_min;
  const double  hfMm = hf_min * hf_max;
  const double  rho_fc = hcMm / hfMm;

  double  gamma = gseed * hc_max * hc_max * rho_fc;

  cs_cdo_advection_set_cip_coef(gamma);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                     const cs_cdo_connect_t       *connect,
                                     const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to CDO
 *         vertex+cell-based schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_initialize(void)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdovcb_threshold = 0.01*cs_math_get_machine_epsilon();

  /* Structure used to build the final system by a cell-wise process */
  cs_cell_condmat = cs_locmat_create(connect->n_max_vbyc);

  /* Specific treatment for handling openMP */
  int  size = cs_glob_n_threads;
  BFT_MALLOC(cs_cdovcb_cell_systems, size, cs_cdo_locsys_t *);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cdovcb_cell_systems[t_id] = cs_cdo_locsys_create(connect->n_max_vbyc+1);
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cdovcb_cell_systems[0] = cs_cdo_locsys_create(connect->n_max_vbyc+1);
#endif /* openMP */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO vertex-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_finalize(void)
{
  /* Free the related discrete Hodge operator associated to a unity property */
  cs_cdovcb_hconf = cs_sla_hmatrix_free(cs_cdovcb_hconf);

  /* Free local structures */
  cs_cell_condmat = cs_locmat_free(cs_cell_condmat);

  for (int i = 0; i < cs_glob_n_threads; i++)
    cs_cdo_locsys_free(&(cs_cdovcb_cell_systems[i]));
  BFT_FREE(cs_cdovcb_cell_systems);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovcb_scaleq_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovcb_scaleq_init(const cs_equation_param_t   *eqp,
                      const cs_mesh_t             *mesh)
{
  /* Sanity checks */
  assert(eqp != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVCB &&
      eqp->var_type != CS_PARAM_VAR_SCAL)
    bft_error(__FILE__, __LINE__, 0, " Invalid type of equation.\n"
              " Expected: scalar-valued CDO vertex+cell-based equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_vertices = connect->v_info->n_elts;
  const cs_lnum_t  n_b_faces = connect->f_info->n_b_elts;
  const cs_lnum_t  n_cells = connect->c_info->n_elts;
  const cs_param_bc_t  *bc_param = eqp->bc;

  cs_cdovcb_scaleq_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_cdovcb_scaleq_t);

  /* Shared pointers */
  b->eqp = eqp;
  b->enforce = bc_param->enforcement;

  /* System dimension */
  b->n_vertices = n_vertices;
  b->n_cells = n_cells;
  b->n_sys = n_vertices + n_cells;
  b->n_max_vcbyc = connect->n_max_vbyc + 1;

  /* Store a direct access to which term one has to compute */
  b->has[CDO_DIFFUSION] = (eqp->flag & CS_EQUATION_DIFFUSION) ? true : false;
  b->has[CDO_ADVECTION] = (eqp->flag & CS_EQUATION_CONVECTION) ? true : false;
  b->has[CDO_REACTION] = (eqp->flag & CS_EQUATION_REACTION) ? true : false;
  b->has[CDO_TIME] = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;
  b->has[CDO_SOURCETERM] = (eqp->n_source_terms > 0) ? true : false;

  /* Store the last computed values of the field at cell centers */
  BFT_MALLOC(b->cell_values, n_cells, double);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++)
    b->cell_values[i] = 0;

  /* Way of storing information to compute the values at cell center and in
     the same time storing the matrix to inverse after assembling */
  b->hybrid_storage = cs_sla_hmatrix_create(n_vertices,
                                            n_cells,
                                            true,  // c2v and v2c are transposed
                                            false, // not symmetric by default
                                            connect->v2v,
                                            connect->c2v);

  /* Initialization of members common to several terms */
  b->flag = 0;
  b->hb = NULL;

  BFT_MALLOC(b->loc_vals, 2*b->n_max_vcbyc, double);
  for (int i = 0; i < 2*b->n_max_vcbyc; i++)
    b->loc_vals[i] = 0;

  /* Diffusion part */
  b->diff = NULL;
  b->diff_pty_uniform = false;
  if (b->has[CDO_DIFFUSION]) {

    bool is_uniform = cs_property_is_uniform(eqp->diffusion_property);

    b->diff_pty_uniform = is_uniform;
    b->diff = cs_cdo_diffusion_builder_init(connect,
                                            CS_SPACE_SCHEME_CDOVCB,
                                            is_uniform,
                                            eqp->diffusion_hodge,
                                            b->enforce);

  }

  /* Advection part */
  b->adv = NULL;
  if (b->has[CDO_ADVECTION])
    b->adv = cs_cdo_advection_builder_init(connect, eqp, b->has[CDO_DIFFUSION]);

  /* Reaction part */
  b->reaction_pty_val = NULL;
  b->reaction_pty_uniform = NULL;
  if (b->has[CDO_REACTION]) {

    if (eqp->reaction_hodge.algo == CS_PARAM_HODGE_ALGO_WBS)
      b->flag |= CS_CDO_BUILD_LOC_HCONF;
    else
      bft_error(__FILE__, __LINE__, 0,
                " Invalid choice of algorithm for the reaction term.");

    BFT_MALLOC(b->reaction_pty_uniform, eqp->n_reaction_terms, bool);
    BFT_MALLOC(b->reaction_pty_val, eqp->n_reaction_terms, double);
    for (int i = 0; i < eqp->n_reaction_terms; i++) {
      b->reaction_pty_val[i] = 0;
      b->reaction_pty_uniform[i] =
        cs_property_is_uniform(eqp->reaction_properties[i]);
    }

  }

  /* Time part */
  b->time_mat_is_diag = false;
  b->time_pty_uniform = false;
  b->time_pty_val = 0.;
  if (b->has[CDO_TIME]) {

    b->time_pty_uniform = cs_property_is_uniform(eqp->time_property);
    if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_VORONOI)
      b->time_mat_is_diag = true;
    else if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_WBS) {
      if (eqp->time_info.do_lumping)
        b->time_mat_is_diag = true;
      else
        b->flag |= CS_CDO_BUILD_LOC_HCONF;
    }

  }

  /* Source term part */
  b->source_terms = NULL;
  if (b->has[CDO_SOURCETERM]) {

    BFT_MALLOC(b->source_terms, b->n_sys, cs_real_t);

    for (int  st_id = 0; st_id < eqp->n_source_terms; st_id++) {
      const cs_source_term_t  *st = eqp->source_terms[st_id];
      if (cs_source_term_get_reduction(st) == CS_SOURCE_TERM_REDUC_PRIM)
        b->flag |= CS_CDO_PRIMAL_SOURCE | CS_CDO_BUILD_HCONF;
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid setting for CDO-V+C scheme.\n"
                  " Source terms have to be reduced on primal entities.");
    }

  } /* There is at least one source term */

  if (b->flag & CS_CDO_BUILD_HCONF ||
      b->flag & CS_CDO_BUILD_LOC_HCONF) {

    cs_param_hodge_t  hwbs_info = {.inv_pty = false,
                                   .type = CS_PARAM_HODGE_TYPE_VC,
                                   .algo = CS_PARAM_HODGE_ALGO_WBS,
                                   .coef = 1.0}; // not useful in this case

    b->hb = cs_hodge_builder_init(connect, hwbs_info);

    if ((b->flag & CS_CDO_BUILD_HCONF) && cs_cdovcb_hconf == NULL)
      _build_hvc_conf(b);

  }

  /* Set members and structures related to the management of the BCs

     Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.
     We also compute also the list of Dirichlet vertices along with their
     related definition.
  */
  b->face_bc = cs_cdo_bc_init(bc_param, n_b_faces);
  b->vtx_dir = cs_cdo_bc_vtx_dir_create(mesh, b->face_bc);

  /* Allocate and initialize dir_val */
  BFT_MALLOC(b->dir_val, b->vtx_dir->n_nhmg_elts, double);
# pragma omp parallel for if (b->vtx_dir->n_nhmg_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->vtx_dir->n_nhmg_elts; i++)
    b->dir_val[i] = 0.0;

  b->c2bcbf_idx = NULL;
  b->c2bcbf_ids = NULL;

  switch (b->enforce) {
  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    if (b->has[CDO_DIFFUSION])
      cs_cdo_diffusion_build_c2bcbf(connect,
                                    b->face_bc->dir,
                                    &(b->c2bcbf_idx),
                                    &(b->c2bcbf_ids));
    break;

  case CS_PARAM_BC_ENFORCE_STRONG:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid choice of enforcement of the boundary conditions.\n"
              " Strong enforcement is not for Vertex-Cell-based schemes.\n"
              " Please modify your settings.");
    break;

  default: // Nothing to do
    break;

  } /* Enforcement of BCs */

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovcb_scaleq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdovcb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovcb_scaleq_free(void   *builder)
{
  cs_cdovcb_scaleq_t  *b = (cs_cdovcb_scaleq_t *)builder;

  if (b == NULL)
    return b;

  /* eqp is only shared. This structure is freed later. */

  BFT_FREE(b->cell_values);
  b->hybrid_storage = cs_sla_hmatrix_free(b->hybrid_storage);

  /* Free builder members common to several terms */
  BFT_FREE(b->loc_vals);
  if (b->hb != NULL)
    b->hb = cs_hodge_builder_free(b->hb);

  /* Free builder sub-structures */
  if (b->has[CDO_DIFFUSION]) {
    b->diff = cs_cdo_diffusion_builder_free(b->diff);

    if (b->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM ||
        b->enforce ==  CS_PARAM_BC_ENFORCE_WEAK_NITSCHE) {
      BFT_FREE(b->c2bcbf_idx);
      BFT_FREE(b->c2bcbf_ids);
    }

  }

  if (b->has[CDO_ADVECTION])
    b->adv = cs_cdo_advection_builder_free(b->adv);

  if (b->has[CDO_REACTION]) {
    BFT_FREE(b->reaction_pty_uniform);
    BFT_FREE(b->reaction_pty_val);
  }

  if (b->has[CDO_SOURCETERM])
    BFT_FREE(b->source_terms);

  /* Free BC structure */
  if (b->vtx_dir->n_nhmg_elts > 0)
    BFT_FREE(b->dir_val);

  b->face_bc = cs_cdo_bc_free(b->face_bc);
  b->vtx_dir = cs_cdo_bc_list_free(b->vtx_dir);

  /* Last free */
  BFT_FREE(b);
  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_sla_matrix_t related to the system to solve
 *
 * \param[in, out]  builder   pointer to a builder structure
 * \param[in, out]  matrix    pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_free_sysmat(void              *builder,
                             cs_sla_matrix_t   *matrix)
{
  cs_cdovcb_scaleq_t  *b = (cs_cdovcb_scaleq_t *)builder;

  if (b == NULL)
    bft_error(__FILE__, __LINE__, 0, " builder structure is empty");

  assert(b->hybrid_storage->xx_block == matrix);

  /* Free matrix */
  matrix = cs_sla_matrix_free(matrix);
  b->hybrid_storage->xx_block = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in, out]  builder     pointer to a cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_compute_source(void   *builder)
{
  cs_desc_t  desc;

  const cs_lnum_t  n_vertices = cs_shared_quant->n_vertices;

  double  *work = cs_equation_get_tmpbuf();
  double  *eval_v = work;
  double  *eval_c = work + n_vertices;
  cs_cdovcb_scaleq_t  *b = (cs_cdovcb_scaleq_t *)builder;

  const cs_equation_param_t  *eqp = b->eqp;

  if (eqp->n_source_terms == 0)
    return;

# pragma omp parallel for if (b->n_sys > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->n_sys; i++)
    b->source_terms[i] = 0;

  /* Sanity check */
  assert(cs_cdovcb_hconf != NULL);

  for (int  st_id = 0; st_id < eqp->n_source_terms; st_id++) {

    const cs_source_term_t  *st = eqp->source_terms[st_id];

    assert(cs_source_term_get_reduction(st) == CS_SOURCE_TERM_REDUC_PRIM);

    double  *mv_v = work + b->n_sys;
    double  *mv_c = mv_v + n_vertices;

    /* Fill eval_v inside this function */
    desc.location = CS_FLAG_SCAL | cs_cdo_primal_vtx;
    desc.state = CS_FLAG_STATE_POTENTIAL;
    cs_source_term_compute(desc, st, &eval_v);

    /* Fill eval_c inside this function */
    desc.location = CS_FLAG_SCAL | cs_cdo_primal_cell;
    desc.state = CS_FLAG_STATE_POTENTIAL;
    cs_source_term_compute(desc, st, &eval_c);

    /* Compute Hodge*eval = st */
    cs_sla_hmatvec(cs_cdovcb_hconf, eval_v, eval_c, &mv_v, &mv_c, true);

    /* Update source term array */
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < b->n_vertices; i++)
      b->source_terms[i] += mv_v[i];

    double  *st_cell = b->source_terms + b->n_vertices;
# pragma omp parallel for if (b->n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < b->n_cells; i++)
      st_cell[i] += mv_c[i];

  } // Loop on source term definitions

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO vertex-based scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the vertex field
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to cs_cdovcb_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] sla_mat    pointer to cs_sla_matrix_t structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_build_system(const cs_mesh_t        *mesh,
                              const cs_real_t        *field_val,
                              double                  dt_cur,
                              void                   *builder,
                              cs_real_t             **rhs,
                              cs_sla_matrix_t       **sla_mat)
{
  cs_cdovcb_scaleq_t  *b = (cs_cdovcb_scaleq_t *)builder;

  const cs_equation_param_t  *eqp = b->eqp;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_flag_t  cmflag = cs_cdovcb_cmflag;

  /* Allocate and initialize a matrix with the larger stencil (that related
     to diffusion => all vertices of a cell are potentially in interaction)
     "adr" means Advection/Diffusion/Reaction
  */
  cs_sla_matrix_t  *sys_mat = b->hybrid_storage->xx_block;
  if (sys_mat == NULL)
    sys_mat = cs_sla_matrix_create_msr_from_index(connect->v2v,
                                                  false,  // symmetric
                                                  true,   // sorted
                                                  1);     // stride

  if (!b->has[CDO_ADVECTION] && b->enforce != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
    sys_mat->flag |= CS_SLA_MATRIX_SYM;

  /* Preparatory step for diffusion term */
  cs_real_33_t  diff_tensor = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  if (b->has[CDO_DIFFUSION]) {

    if (b->diff_pty_uniform)
      cs_property_get_cell_tensor(0, // cell_id
                                  eqp->diffusion_property,
                                  eqp->diffusion_hodge.inv_pty,
                                  diff_tensor);

  } // DIFFUSION

  if (b->has[CDO_ADVECTION]) {

    cmflag |= CS_CDO_LOCAL_EF;
    _set_cip_coef(eqp->advection_field);

  }

  /* Preparatory step for unsteady term */
  if (b->has[CDO_TIME])
    if (b->time_pty_uniform)
      b->time_pty_val = cs_property_get_cell_value(0, eqp->time_property);

  /* Preparatory step for reaction term */
  if (b->has[CDO_REACTION]) {

    for (int r = 0; r < eqp->n_reaction_terms; r++) {
      if (b->reaction_pty_uniform[r]) {
        cs_property_t  *r_pty = eqp->reaction_properties[r];
        b->reaction_pty_val[r] = cs_property_get_cell_value(0, r_pty);
      }
    }

  } // REACTION

  /* Initialize full rhs */
  cs_real_t  *full_rhs = *rhs;
  if (full_rhs == NULL)
    BFT_MALLOC(full_rhs, b->n_sys, cs_real_t);
  cs_real_t  *cell_rhs = full_rhs + b->n_vertices;

# pragma omp parallel for if (b->n_sys > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->n_sys; i++)
    full_rhs[i] = 0.0;

  /* Add the contribution of source terms to the full rhs for this time step */
  _add_source_terms(b, full_rhs);

  /* Compute the values of the Dirichlet BC.
     TODO: do the analogy for Neumann BC */
  _compute_dir_values(mesh, field_val, b);

  /* Temporary pre-allocated buffers (the following buffer is used in
     _add_source_terms (be careful if the order of calls is changed) */
  cs_real_t  *dir_bc_vals = cs_equation_get_tmpbuf();
  cs_flag_t  *cell_flag = connect->c_info->flag;

   /* Initialize arrays */
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->n_vertices; i++)
    dir_bc_vals[i] = 0;

  /* Store the Dirichlet values into an array of size n_vertices */
# pragma omp parallel for if (b->vtx_dir->n_nhmg_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->vtx_dir->n_nhmg_elts; i++)
    dir_bc_vals[b->vtx_dir->elt_ids[i]] = b->dir_val[i];

  // Preparatory step fur a future use of openMP
  int  t_id = 0;
  cs_cell_mesh_t *mesh_c = cs_cdo_local_get_cell_mesh(t_id);
  cs_cdo_locsys_t  *csys = cs_cdovcb_cell_systems[t_id];

  /* Main loop on cells to build the linear system */
  /* --------------------------------------------- */

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Set the local mesh structure for the current cell */
    cs_cell_mesh_build(c_id, cmflag, connect, quant, mesh_c);

    /* Cell-wise view of the linear system to build */
    const int  n_vc = mesh_c->n_vc;
    const int  n_sysc = n_vc + 1;

    cs_locmat_t  *hconf_cell = NULL;

    /* Store the local values attached to Dirichlet values if the current cell
       has at least one border face */
    if (cell_flag[c_id] & CS_CDO_CONNECT_BD)
      for (short int v = 0; v < n_vc; v++)
        csys->dir_bc[v] = dir_bc_vals[mesh_c->v_ids[v]];

    /* Initialize the local matrix storing advection/diffusion/reaction terms */
    csys->mat->n_ent = n_sysc;
    for (short int v = 0; v < n_vc; v++) {
      csys->mat->ids[v] = mesh_c->v_ids[v];
      csys->rhs[v] = 0.;
    }
    csys->mat->ids[n_vc] = c_id;
    csys->rhs[n_vc] = cell_rhs[c_id]; /* Store the rhs_c (used for static
                                       condensation) */

    for (short int v = 0; v < n_sysc*n_sysc; v++)
      csys->mat->val[v] = 0;

    /* DIFFUSION TERM */
    if (b->has[CDO_DIFFUSION]) { /* Define the local stiffness matrix */

      if (b->diff_pty_uniform == false)
        cs_property_get_cell_tensor(c_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    diff_tensor);

      cs_locmat_t  *diff_mat = // local matrix owned by the diffusion builder
        cs_cdo_diffusion_build_local(quant,
                                     mesh_c,
                                     (const cs_real_3_t (*))diff_tensor,
                                     b->diff);

      cs_locmat_add(csys->mat, diff_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
      bft_printf(">> Local diffusion matrix");
      cs_locmat_dump(c_id, diff_mat);
#endif

      /* Weakly enforced Dirichlet BCs for cells attached to the boundary */
      if (b->c2bcbf_idx != NULL && cell_flag[c_id] & CS_CDO_CONNECT_BD) {

        for (cs_lnum_t j = b->c2bcbf_idx[c_id];
             j < b->c2bcbf_idx[c_id+1]; j++) {

          /* csys is updated inside (matrix and rhs) */
          cs_cdo_diffusion_weak_bc(b->c2bcbf_ids[j], // border face id
                                   mesh_c,
                                   (const cs_real_3_t (*))diff_tensor,
                                   b->diff,
                                   csys); // csys is updated inside

        } // Loop on border faces

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 2
        bft_printf(">> Local diffusion matrix after weak enforcement");
        cs_locmat_dump(c_id, diff_mat);
#endif
      } /* Weak enforcement of Dirichlets BCs */

    } /* DIFFUSION */

    /* ADVECTION TERM */
    if (b->has[CDO_ADVECTION]) { /* Define the local advection matrix */

      cs_locmat_t  *adv_mat = cs_cdovcb_advection_build(mesh_c, eqp, b->adv);

      cs_locmat_add(csys->mat, adv_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
      bft_printf(">> Local advection matrix");
      cs_locmat_dump(c_id, adv_mat);
#endif

      /* Last treatment for the advection term: Apply boundary conditions
         csys is updated inside (matrix and rhs) */
      if (cell_flag[c_id] & CS_CDO_CONNECT_BD)
        cs_cdovcb_advection_add_bc(mesh_c, eqp, b->adv, csys);

    } /* ADVECTION */

    if (b->flag & CS_CDO_BUILD_LOC_HCONF)
      hconf_cell = cs_hodge_build_cellwise(mesh_c, b->hb);

    /* REACTION TERM */
    if (b->has[CDO_REACTION]) { /* Define the local reaction matrix */

      double  rpty_val = 0;
      for (int r = 0; r < eqp->n_reaction_terms; r++) // Loop on reaction terms
        if (b->reaction_pty_uniform[r])
          rpty_val += b->reaction_pty_val[r];
        else
          rpty_val += cs_property_get_cell_value(c_id,
                                                 eqp->reaction_properties[r]);

      /* Update local system matrix with the reaction term */
      cs_locmat_mult_add(csys->mat, rpty_val, hconf_cell);

    } /* REACTION */

    /* TIME CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
    if (b->has[CDO_TIME]) {

      /* Get the value of the time property */
      double  tpty_val = 0;
      if (b->time_pty_uniform)
        tpty_val = b->time_pty_val/dt_cur;
      else
        tpty_val = cs_property_get_cell_value(c_id, eqp->time_property)/dt_cur;

      /* Apply the time discretization to the local system.
         Update csys (matrix and rhs) */
      _apply_time_scheme(tpty_val, field_val, mesh_c, hconf_cell, b, csys);

    } /* Time contribution */

    /* Sore and condense the local system matrix of size n_vc + 1 into
       a matrix of size n_vc. Store information in hybrid_storage in order to
       be able to compute the values at cell centers.
       rhs is updated before assembling. */
    _condense_and_store(b, csys, cs_cell_condmat);

    /* Assemble the matrix related to the advcetion/diffusion/reaction terms
       If advection is activated, the resulting system is not symmetric
       Otherwise, the system is symmetric with extra-diagonal terms. */
    if (sys_mat->flag & CS_SLA_MATRIX_SYM)
      cs_sla_assemble_msr_sym(cs_cell_condmat, sys_mat, false);
    else
      cs_sla_assemble_msr(cs_cell_condmat, sys_mat);

    /* Assemble the right-hand side (rhs) */
    for (short int v = 0; v < n_vc; v++)
      full_rhs[mesh_c->v_ids[v]] += csys->rhs[v];
    cell_rhs[c_id] = csys->rhs[n_vc];

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 0
    bft_printf(">> (FINAL) Local system matrix");
    cs_locmat_dump(c_id, csys->mat);
#endif

  } // Main loop on cells

  /* Final step in BC management.
     Must be call after the application of the time scheme. */
  if (b->enforce == CS_PARAM_BC_ENFORCE_WEAK_PENA) {

    const cs_cdo_bc_list_t  *vtx_dir = b->vtx_dir;

    // Advanced parameters
    const cs_real_t  penalization_coef = 1e-2/cs_math_get_machine_epsilon();

# pragma omp parallel for if (vtx_dir->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < vtx_dir->n_elts; i++)
      sys_mat->diag[vtx_dir->elt_ids[i]] += penalization_coef;

# pragma omp parallel for if (vtx_dir->n_nhmg_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < vtx_dir->n_nhmg_elts; i++)
      full_rhs[vtx_dir->elt_ids[i]] += penalization_coef * b->dir_val[i];

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
  cs_sla_system_dump("system-BeforeClean.log", NULL, sys_mat, full_rhs);
#endif

  /* Clean matrix (set entries to zero if there are below a given threshold
     in order to improve the matrix conditionning */
  cs_sla_matrix_clean(eqp->verbosity, cs_cdovcb_threshold, sys_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
  cs_sla_system_dump("system-AfterClean.log", NULL, sys_mat, full_rhs);
#endif

  /* Return pointers */
  *rhs = full_rhs;
  *sla_mat = sys_mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in, out] builder    pointer to builder structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_update_field(const cs_real_t     *solu,
                              const cs_real_t     *rhs,
                              void                *builder,
                              cs_real_t           *field_val)
{
  cs_cdovcb_scaleq_t  *b = (cs_cdovcb_scaleq_t  *)builder;

  /* Sanity check */
  assert(b->enforce != CS_PARAM_BC_ENFORCE_STRONG);

  /* Set the values at vertices */
  for (cs_lnum_t i = 0; i < b->n_vertices; i++)
    field_val[i] = solu[i];

  /* Compute values at cells pc = acc^-1*(sc - acv*pv) */
  const cs_connect_index_t  *c2v = cs_shared_connect->c2v;
  const double  *cv_vals = b->hybrid_storage->cx_vals;
  const double  *cc_vals = b->hybrid_storage->cc_diag;
  const double  *cell_rhs = rhs + b->n_vertices;

  for (cs_lnum_t c_id = 0; c_id < b->n_cells; c_id++) {

    double  v_contrib = 0.;
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
      v_contrib += cv_vals[j]*field_val[c2v->ids[j]];

    b->cell_values[c_id] = cc_vals[c_id] * (cell_rhs[c_id] - v_contrib);

  } // Loop on cells
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at cell centers (DoF used in the linear
 *         system are located at primal vertices and field related to the
 *         structure equation is also attached to primal vertices
 *
 * \param[in]  builder    pointer to a builder structure
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

double *
cs_cdovcb_scaleq_get_cell_values(const void          *builder)
{
  const cs_cdovcb_scaleq_t  *b = (const cs_cdovcb_scaleq_t  *)builder;

  if (b == NULL)
    return NULL;
  else
    return b->cell_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a list of faces
 *
 * \param[in]       builder    pointer to a builder structure
 * \param[in]       pdi        pointer to an array of field values
 * \param[in]       ml_id      id related to a cs_mesh_location_t struct.
 * \param[in]       direction  indicate in which direction flux is > 0
 * \param[in, out]  diff_flux  pointer to the value of the diffusive flux
 * \param[in, out]  conv_flux  pointer to the value of the convective flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_compute_flux_across_plane(const void          *builder,
                                           const cs_real_t     *pdi,
                                           int                  ml_id,
                                           const cs_real_t      direction[],
                                           double              *diff_flux,
                                           double              *conv_flux)
{
  double  flx, p_f;
  cs_real_33_t  pty_tens;
  cs_nvec3_t  adv_c;

  *diff_flux = 0.;
  *conv_flux = 0.;

  if (pdi == NULL)
    return;

  cs_mesh_location_type_t  ml_t = cs_mesh_location_get_type(ml_id);

  if (ml_t != CS_MESH_LOCATION_INTERIOR_FACES &&
      ml_t != CS_MESH_LOCATION_BOUNDARY_FACES) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_(" Mesh location type is incompatible with the computation\n"
                 " of the flux across faces.\n"));
    return;
  }

  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  if (n_elts[0] > 0 && elt_ids == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Computing the flux across all interior or border faces is not"
                " managed yet."));

  const cs_cdovcb_scaleq_t  *b = (const cs_cdovcb_scaleq_t  *)builder;
  const cs_equation_param_t  *eqp = b->eqp;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_sla_matrix_t  *f2c = connect->f2c;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  // To be modified for an fully integration of openMP
  cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(0);
  double  *p_v = cs_equation_get_tmpbuf(); // used as a temporary buffer

  if (ml_t == CS_MESH_LOCATION_BOUNDARY_FACES) {

    const cs_lnum_t  n_i_faces = connect->f_info->n_i_elts;
    const cs_lnum_t  shift_if = 2*n_i_faces;

    for (cs_lnum_t id = 0; id < n_elts[0]; id++) {

      const cs_lnum_t  bf_id = elt_ids[id];
      const cs_lnum_t  f_id = n_i_faces + bf_id;
      const cs_lnum_t  c_id = f2c->col_id[shift_if + bf_id];

      /* Build a face-wise view of the mesh */
      cs_face_mesh_build(c_id, f_id, connect, quant, fm);

      const short int  sgn = (_dp3(fm->face.unitv, direction) < 0) ? -1 : 1;

      /* Store values related to this face */
      for (short int v = 0; v < fm->n_vf; v++)
        p_v[v] = pdi[fm->v_ids[v]];

      cs_reco_pv_at_face_center(fm, p_v, &p_f);

      if (b->has[CDO_DIFFUSION]) { /* Compute the local diffusive flux */

        cs_property_get_cell_tensor(c_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);

        flx = cs_cdo_diffusion_face_flux(fm,
                                         (const cs_real_3_t (*))pty_tens,
                                         p_v, p_f, b->cell_values[c_id],
                                         b->diff);
        *diff_flux += sgn * flx;

      } // Diffusive flux

      if (b->has[CDO_ADVECTION]) { /* Compute the local advective flux */

        const double  coef = sgn * fm->face.meas * p_f;

        cs_advection_field_get_cell_vector(c_id, eqp->advection_field, &adv_c);
        *conv_flux += coef * adv_c.meas * _dp3(adv_c.unitv, fm->face.unitv);

      }

    } // Loop on selected border faces

  }
  else if (ml_t == CS_MESH_LOCATION_INTERIOR_FACES) {

    for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

      const cs_lnum_t  f_id = elt_ids[i];
      const cs_lnum_t  shift_f = 2*f_id;
      const cs_lnum_t  c1_id = f2c->col_id[shift_f];
      const cs_lnum_t  c2_id = f2c->col_id[shift_f+1];

      /* Build a face-wise view of the mesh */
      cs_face_mesh_build(c1_id, f_id, connect, quant, fm);

      const short int  sgn = (_dp3(fm->face.unitv, direction) < 0) ? -1 : 1;

      /* Store values related to this face */
      for (short int v = 0; v < fm->n_vf; v++)
        p_v[v] = pdi[fm->v_ids[v]];

      cs_reco_pv_at_face_center(fm, p_v, &p_f);

      if (b->has[CDO_DIFFUSION]) { /* Compute the local diffusive flux */

        /* Compute the diffusive flux seen from cell c1 */
        cs_property_get_cell_tensor(c1_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);

        flx = cs_cdo_diffusion_face_flux(fm,
                                         (const cs_real_3_t (*))pty_tens,
                                         p_v, p_f, b->cell_values[c1_id],
                                         b->diff);

        /* Compute the diffusive flux seen from cell c2 */
        cs_face_mesh_build(c2_id, f_id, connect, quant, fm);

        cs_property_get_cell_tensor(c2_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);

        flx += cs_cdo_diffusion_face_flux(fm,
                                          (const cs_real_3_t (*))pty_tens,
                                          p_v, p_f, b->cell_values[c2_id],
                                          b->diff);

        *diff_flux += sgn * 0.5 * flx;

      } // Diffusive flux

      if (b->has[CDO_ADVECTION]) { /* Compute the local advective flux */

        /* Advective flux seen from cell 1 */
        cs_advection_field_get_cell_vector(c1_id, eqp->advection_field, &adv_c);
        flx = adv_c.meas * _dp3(adv_c.unitv, fm->face.unitv);

        /* Advective flux seen from cell 2 */
        cs_advection_field_get_cell_vector(c2_id, eqp->advection_field, &adv_c);
        flx += adv_c.meas * _dp3(adv_c.unitv, fm->face.unitv);

        *conv_flux += sgn * 0.5 * flx  * p_f * fm->face.meas;

      } // Advective flux

    } // Loop on selected interior faces

  } // Set of interior or border faces

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across all dual faces.
 *
 * \param[in]       values      discrete values for the potential
 * \param[in, out]  builder     pointer to builder structure
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_cellwise_diff_flux(const cs_real_t   *values,
                                    void              *builder,
                                    cs_real_t         *diff_flux)
{
  cs_cdovcb_scaleq_t  *b = (cs_cdovcb_scaleq_t  *)builder;

  const cs_equation_param_t  *eqp = b->eqp;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_connect_index_t  *c2e = connect->c2e;

  if (!b->has[CDO_DIFFUSION]) { // No diffusion
    for (cs_lnum_t i = 0; i < c2e->idx[quant->n_cells]; i++)
      diff_flux = 0;
    return;
  }

  /* Sanity checks */
  assert(eqp->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_WBS);

  /* Diffusion tensor */
  cs_real_33_t  diff_tensor = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  if (b->diff_pty_uniform)
    cs_property_get_cell_tensor(0,  // cell_id
                                eqp->diffusion_property,
                                eqp->diffusion_hodge.inv_pty,
                                diff_tensor);

  // TODO: Prepare a full openMP implementation
  cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(0);

  /* Retrieve temporary buffers */
  double  *p_v = cs_equation_get_tmpbuf(); // used as a temporary buffer

  /* Define the flux by cellwise contributions */
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    if (b->diff_pty_uniform == false)
      cs_property_get_cell_tensor(c_id,
                                  eqp->diffusion_property,
                                  eqp->diffusion_hodge.inv_pty,
                                  diff_tensor);

    /* Set the local mesh structure for the current cell */
    cs_cell_mesh_build(c_id, cs_cdovcb_cmflag, connect, quant, cm);

    /* Define a local buffer keeping the value of the discrete potential
       for the current cell */
    for (short int v = 0; v < cm->n_vc; v++)
      p_v[v] = values[cm->v_ids[v]];

    cs_cdo_diffusion_cellwise_flux(cm,
                                   quant->dface + c2e->idx[c_id],
                                   (const cs_real_3_t (*))diff_tensor,
                                   p_v,
                                   b->cell_values[c_id], // p_c
                                   b->diff,
                                   diff_flux + c2e->idx[c_id]);

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field strufcture
 * \param[in, out]  builder    pointer to builder structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_extra_op(const char            *eqname,
                          const cs_field_t      *field,
                          void                  *builder)
{
  cs_cdovcb_scaleq_t  *b = (cs_cdovcb_scaleq_t  *)builder;

  const cs_equation_param_t  *eqp = b->eqp;

  // TODO
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
