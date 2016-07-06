/*============================================================================
 * Build an algebraic CDO vertex-based system for unsteady convection diffusion
 * reaction scalar equations with source terms
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

#include "cs_cdovb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVB_SCALEQ_DBG  0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/* Algebraic system for CDO vertex-based discretization */

struct _cs_cdovb_scaleq_t {

  /* Pointer to a cs_equation_param_t structure shared with a cs_equation_t
     structure.  */

  const cs_equation_param_t  *eqp;

  /* System size (known boundary entities may be removed if BCs are strongly
     enforced) */
  cs_lnum_t  n_vertices;
  cs_lnum_t  n_dof_vertices; /* n_rows = n_cols = n_vertices - dir. vertices */

  /* Shortcut to know what to build */
  bool                   has[N_CDO_TERMS];
  cs_flag_t              flag;

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

  /* Indirection between zipped numbering (without BC) and initial numbering
     Allocated only if the boundary conditions are strongly enforced. */
  cs_lnum_t     *v_z2i_ids;  // Mapping n_dof_vertices -> n_vertices
  cs_lnum_t     *v_i2z_ids;  // Mapping n_vertices     -> n_dof_vertices

};

/*============================================================================
 * Private variables
 *============================================================================*/

static double  cs_cdovb_threshold = 1e-12; // Set during initialization
static cs_sla_matrix_t  *cs_cdovb_hconf = NULL;
static cs_cdo_locsys_t  *cs_cell_sys = NULL;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_time_step_t  *cs_shared_time_step;

/* Flag to indicate which members have to be built in a cs_cell_mesh_t
   structure */
static const cs_flag_t  cs_cdovb_cmflag =
  CS_CDO_LOCAL_V | CS_CDO_LOCAL_E | CS_CDO_LOCAL_EV;

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
_add_source_terms(cs_cdovb_scaleq_t     *b,
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
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_vertices; i++)
        full_rhs[i] += b->source_terms[i];

    else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
             t_info.scheme == CS_TIME_SCHEME_THETA) {

      const double  tcoef = 1 - t_info.theta;

# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_vertices; i++)
        full_rhs[i] += tcoef * b->source_terms[i];

    }

    /* Update b->source_term with the value attached to t_cur */
    cs_cdovb_scaleq_compute_source(b);

    if (t_info.scheme == CS_TIME_SCHEME_IMPLICIT)
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_vertices; i++)
        full_rhs[i] += b->source_terms[i];

    else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
             t_info.scheme == CS_TIME_SCHEME_THETA) {
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_vertices; i++)
        full_rhs[i] += t_info.theta * b->source_terms[i];

    }

  }
  else { /* Steady case: source terms have already been computed during
            the initialization step */

# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < b->n_vertices; i++)
      full_rhs[i] += b->source_terms[i];

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a discrete Hodge op. Vp-->Cd using conforming reco. op.
 *
 * \param[in, out] b          pointer to a cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_hvpcd_conf(cs_cdovb_scaleq_t        *b)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Default flag value for vertex-based scalar equations */
  cs_flag_t  cm_flag = cs_cdovb_cmflag | CS_CDO_LOCAL_F | CS_CDO_LOCAL_FE;

  // To be modified for an fully integration of openMP
  cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(0);

  /* Initialize a matrix structure */
  cs_cdovb_hconf = cs_sla_matrix_create_msr_from_index(connect->v2v,
                                                       true, true, 1);

  /* Cellwise construction ==> Loop on cells */
  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Set the local mesh structure for the current cell */
    cs_cell_mesh_build(c_id, cm_flag, connect, quant, cm);

    /* Build the local dense matrix related to this operator */
    cs_locmat_t  *hloc = cs_hodge_build_cellwise(cm, b->hb);

    /* Assemble the cellwise matrix into the "global" matrix */
    cs_sla_assemble_msr_sym(hloc, cs_cdovb_hconf, false);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the time discretization to the local system
 *
 * \param[in, out] tpty_val   current value of the time property
 * \param[in]      loc_fval   pointer to the current value of the field
 * \param[in]      cm         pointer to a cs_locmesh_t structure
 * \param[in]      loc_hconf  pointer to a conforming discrete Hodge op.
 * \param[in, out] b          pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] loc_sys    pointer to a cs_locmat_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_apply_time_scheme(double                   tpty_val,
                   const cs_real_t         *field_val,
                   const cs_cell_mesh_t    *cm,
                   const cs_locmat_t       *loc_hconf,
                   cs_cdovb_scaleq_t       *b,
                   cs_cdo_locsys_t         *loc_sys)

{
  const cs_param_time_t  t_info = b->eqp->time_info;

  cs_locmat_t  *loc_mat = loc_sys->mat;
  double  *loc_rhs = loc_sys->rhs;

  /* Temporary buffers of size equal to the number of cell vertices */
  double  *fval = b->loc_vals;
  double  *adr_pn = b->loc_vals + cm->n_vc;

  /* Set the values of the fields attached to this cell */
  for (short int v = 0; v < cm->n_vc; v++)
    fval[v] = field_val[cm->v_ids[v]];

  if (b->time_mat_is_diag) { // Lumped or Voronoi --> only diag.

    tpty_val *= cm->vol_c; // |c| * wvc = |dual_cell(v) cap c|

    switch (t_info.scheme) {

    case CS_TIME_SCHEME_EXPLICIT:

      /* Compute (Adv + Dif + Rea)*p^(n) */
      cs_locmat_matvec(loc_mat, fval, adr_pn);

      /* Reset the local system matrix to assemble */
      for (short int v = 0; v < cm->n_vc*cm->n_vc; v++)
        loc_mat->val[v] = 0;

      /* Update the rhs with -(Adv + Dif + Rea)*p^(n) + TimeMat*p^(n)
         Set the local matrix to a (time) diagonal matrix */
      for (short int v = 0; v < cm->n_vc; v++) {

        const double dval_v = tpty_val * cm->wvc[v]; // |c|*wvc=|dcell(v) cap c|

        loc_mat->val[v*loc_mat->n_ent + v] = dval_v;
        loc_rhs[v] += -adr_pn[v] + dval_v * fval[v];

      }
      break;

    case CS_TIME_SCHEME_IMPLICIT:

      /* Update the rhs with TimeMat*p^(n)
         Update the local matrix with the time diagonal matrix */
      for (short int v = 0; v < cm->n_vc; v++) {

        const double dval_v = tpty_val * cm->wvc[v];

        loc_mat->val[v*loc_mat->n_ent + v] += dval_v;
        loc_rhs[v] += dval_v * fval[v];

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
        for (short int vi = 0; vi < cm->n_vc; vi++) {

          const double dval_v = tpty_val * cm->wvc[vi];

          loc_rhs[vi] += tcoef * adr_pn[vi] + dval_v * fval[vi];

          double  *mi = loc_mat->val + vi*loc_mat->n_ent;
          for (short int vj = 0; vj < cm->n_vc; vj++) {

            mi[vj] *= t_info.theta;
            if (vi == vj)
              mi[vj] += dval_v;

          } // Loop on cell vertices (vj)

        } // Loop on cell vertices (vi)

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
      for (short int v = 0; v < cm->n_vc; v++)
        loc_rhs[v] -= adr_pn[v];

      /* Replace the local system matrix by that of time */
      for (short int vi = 0; vi < cm->n_vc; vi++) {

        double  *mi = loc_mat->val + vi*loc_mat->n_ent;
        double  *hi = loc_hconf->val + vi*loc_hconf->n_ent;

        for (short int vj = 0; vj < cm->n_vc; vj++)
          mi[vj] = tpty_val * hi[vj];

      }

      /* Update the rhs with TimeMat*p^(n) */
      cs_locmat_matvec(loc_mat, fval, adr_pn);
      for (short int v = 0; v < cm->n_vc; v++)
        loc_rhs[v] += adr_pn[v];

      break;

    case CS_TIME_SCHEME_IMPLICIT:

      /* Update the rhs with TimeMat*p^(n) */
      cs_locmat_matvec(loc_hconf, fval, adr_pn);
      for (short int v = 0; v < cm->n_vc; v++)
        loc_rhs[v] += tpty_val * adr_pn[v];

      /* Update the local system with the time matrix */
      for (short int vi = 0; vi < cm->n_vc; vi++) {

        double  *mi = loc_mat->val + vi*loc_mat->n_ent;
        double  *hi = loc_hconf->val + vi*loc_hconf->n_ent;

        for (short int vj = 0; vj < cm->n_vc; vj++)
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
        for (short int v = 0; v < cm->n_vc; v++)
          loc_rhs[v] += tcoef * adr_pn[v];

        /* Update the rhs with TimeMat*p^(n) */
        cs_locmat_matvec(loc_hconf, fval, adr_pn);
        for (short int v = 0; v < cm->n_vc; v++)
          loc_rhs[v] += tpty_val * adr_pn[v];

        /* Update the local system with the time diagonal matrix */
        for (short int vi = 0; vi < cm->n_vc; vi++) {

          double  *mi = loc_mat->val + vi*loc_mat->n_ent;
          double  *hi = loc_hconf->val + vi*loc_hconf->n_ent;

          for (short int vj = 0; vj < cm->n_vc; vj++) {

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
 * \brief   Compute the values of the Dirichlet BCs
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      field_val    pointer to the current value of the field
 * \param[in, out] builder      pointer to a cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_dir_values(const cs_mesh_t            *mesh,
                    const cs_real_t            *field_val,
                    const cs_cdovb_scaleq_t    *builder)
{
  cs_lnum_t  i;

  const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;
  const cs_equation_param_t  *eqp = builder->eqp;

  if (vtx_dir->n_nhmg_elts == 0)
    return; // Nothing to do

  cs_flag_t  dof_flag = cs_cdo_primal_vtx | CS_FLAG_SCAL;

  /* Get the value of the Dirichlet for the current time */
  cs_cdo_bc_dirichlet_set(dof_flag,
                          cs_shared_time_step,
                          mesh,
                          eqp->bc,
                          vtx_dir,
                          builder->dir_val);

  if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
      builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
    if (builder->has[CDO_TIME]) {

      const cs_param_time_t  t_info = eqp->time_info;

      /* Previous values of the unknown are stored inside field_val (iter n) */
      if (t_info.scheme == CS_TIME_SCHEME_EXPLICIT)
        for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
          builder->dir_val[i] = field_val[vtx_dir->elt_ids[i]];

      else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
               t_info.scheme == CS_TIME_SCHEME_THETA) {

        const double  tcoef = 1 - t_info.theta;

        for (i = 0; i < vtx_dir->n_nhmg_elts; i++) {
          builder->dir_val[i] *= t_info.theta;
          builder->dir_val[i] += tcoef * field_val[vtx_dir->elt_ids[i]];
        }

      }

    } /* Unsteady */

  } /* Enforcement is not strong or penalized */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply boundary conditions. Update right-hand side and the system
 *          matrix
 *
 * \param[in, out] bld          pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] full_rhs     right-hand side
 * \param[in, out] full_matrix  matrix of the linear system
 */
/*----------------------------------------------------------------------------*/

static void
_strong_bc_enforcement(cs_cdovb_scaleq_t       *bld,
                       cs_real_t              **rhs,
                       cs_sla_matrix_t        **matrix)
{
  const cs_cdo_bc_list_t  *vtx_dir = bld->vtx_dir;

  if (vtx_dir->n_nhmg_elts == 0)
    return;

  /* Sanity check */
  assert(bld->n_vertices > bld->n_dof_vertices);

  cs_sla_matrix_t  *full_matrix = *matrix, *final_matrix = NULL;
  double  *full_rhs = *rhs, *final_rhs = NULL;
  double  *work = cs_equation_get_tmpbuf();
  double  *tmp_rhs = work;
  double  *x_bc = work + bld->n_vertices;
  double  *contrib = work + 2*bld->n_vertices;

# pragma omp parallel for if (bld->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < bld->n_vertices; i++)
    x_bc[i] = 0.0;
# pragma omp parallel for if (vtx_dir->n_nhmg_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < vtx_dir->n_nhmg_elts; i++)
    x_bc[vtx_dir->elt_ids[i]] = bld->dir_val[i];

  /* Compute full_matrix*Tbc: rhs = rhs - full_matrix*Tbc */
  cs_sla_matvec(full_matrix, x_bc, &contrib, true);
# pragma omp parallel for if (bld->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < bld->n_vertices; i++)
    full_rhs[i] -= contrib[i];

  /* Reduce the rhs size. Remove vertices with Dirichlet BC */
  memcpy(tmp_rhs, full_rhs, bld->n_vertices*sizeof(double));
  final_rhs = *rhs;
# pragma omp parallel for if (bld->n_dof_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < bld->n_dof_vertices; i++)
    final_rhs[i] = tmp_rhs[bld->v_z2i_ids[i]];

  /* Reduce the system size from n_vertices to n_dof_vertices.
     Vertices attached to a Dirichlet BC are removed.
     Extract block with degrees of freedom */
  final_matrix = cs_sla_matrix_pack(bld->n_dof_vertices,
                                    bld->n_dof_vertices,
                                    full_matrix,
                                    bld->v_z2i_ids,
                                    bld->v_i2z_ids,
                                    true); // keep sym.

  /* Free buffers */
  full_matrix = cs_sla_matrix_free(full_matrix);

  /* Return pointers */
  *matrix = final_matrix;
  *rhs = final_rhs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Modify the matrix of the linear system and its right hand side to
 *          take into account a strong enforcement or a large penalization of
 *          the boundary conditions.
 *          Nothing to do in case of weak enforcement.
 *
 * \param[in, out]  bld       pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out]  rhs       right-hand side
 * \param[in, out]  matrix    matrix of the linear system
 */
/*----------------------------------------------------------------------------*/

static void
_enforce_bc(cs_cdovb_scaleq_t          *bld,
            cs_real_t                 **rhs,
            cs_sla_matrix_t           **matrix)
{
  /* Sanity check */
  if (bld->enforce != CS_PARAM_BC_ENFORCE_STRONG &&
      bld->n_vertices != bld->n_dof_vertices)
    bft_error(__FILE__, __LINE__, 0,
              " Error detected: Boundary conditions are not strongly enforced"
              " but there are some removed vertices.");

  cs_sla_matrix_t  *full_matrix = *matrix;
  double  *full_rhs = *rhs;

  const cs_cdo_bc_list_t  *vtx_dir = bld->vtx_dir;

  /* Treatment differs according to the way of enforcing BCs.
     In vertex-based scheme, Dirichlet BC are essential and Neuman BC natural */
  switch (bld->enforce) {

  case CS_PARAM_BC_ENFORCE_STRONG:
    _strong_bc_enforcement(bld, rhs, matrix);
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_PENA:
    {
      // Advanced parameters
      const cs_real_t  penalization_coef = 1e-2/cs_math_get_machine_epsilon();

# pragma omp parallel for if (vtx_dir->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < vtx_dir->n_elts; i++)
        full_matrix->diag[vtx_dir->elt_ids[i]] += penalization_coef;

# pragma omp parallel for if (vtx_dir->n_nhmg_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < vtx_dir->n_nhmg_elts; i++)
        full_rhs[vtx_dir->elt_ids[i]] += penalization_coef * bld->dir_val[i];

    }
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    // Nothing to do
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " This kind of BC enforcement is not implemented yet.\n"
              " Please modify your settings.");

  } // End of switch (enforcement)

  /* TODO: Add contribution for Neumann BC (if homogeneous nothing to do)
     and Robin BC */

  if (bld->n_vertices == bld->n_dof_vertices) { // Keep the full system
    *matrix = full_matrix;
    *rhs = full_rhs;
  }

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
cs_cdovb_scaleq_set_shared_pointers(const cs_cdo_quantities_t    *quant,
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
 *         vertex-based schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_initialize(void)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdovb_threshold = 0.01*cs_math_get_machine_epsilon();

  /* Structure used to build the final system by a cell-wise process */
  cs_cell_sys = cs_cdo_locsys_create(connect->n_max_vbyc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO vertex-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_finalize(void)
{
  /* Free Hodge operator defined from conforming reconstruction op. */
  cs_cdovb_hconf = cs_sla_matrix_free(cs_cdovb_hconf);

  /* Free local structures */
  cs_cdo_locsys_free(&cs_cell_sys);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_scaleq_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_scaleq_init(const cs_equation_param_t   *eqp,
                     const cs_mesh_t             *mesh)
{
  /* Sanity checks */
  assert(eqp != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVB &&
      eqp->var_type != CS_PARAM_VAR_SCAL)
    bft_error(__FILE__, __LINE__, 0, " Invalid type of equation.\n"
              " Expected: scalar-valued CDO vertex-based equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_vertices = connect->v_info->n_elts;
  const cs_lnum_t  n_b_faces = connect->f_info->n_b_elts;
  const cs_param_bc_t  *bc_param = eqp->bc;

  cs_cdovb_scaleq_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_cdovb_scaleq_t);

  /* Shared pointers */
  b->eqp = eqp;
  b->enforce = bc_param->enforcement;

  /* Dimensions: By default, we set number of DoFs as if there is a weak
     enforcement of the boundary conditions */
  b->n_vertices = n_vertices;
  b->n_dof_vertices = n_vertices;

  /* Store a direct access to which term one has to compute */
  b->has[CDO_DIFFUSION] = (eqp->flag & CS_EQUATION_DIFFUSION) ? true : false;
  b->has[CDO_ADVECTION] = (eqp->flag & CS_EQUATION_CONVECTION) ? true : false;
  b->has[CDO_REACTION] = (eqp->flag & CS_EQUATION_REACTION) ? true : false;
  b->has[CDO_TIME] = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;
  b->has[CDO_SOURCETERM] = (eqp->n_source_terms > 0) ? true : false;

  /* Initialization of members common to several terms */
  b->flag = 0;
  b->hb = NULL;

  BFT_MALLOC(b->loc_vals, 2*connect->n_max_vbyc, double);
  for (int i = 0; i < 2*connect->n_max_vbyc; i++)
    b->loc_vals[i] = 0;

  /* Diffusion part */
  b->diff = NULL;
  b->diff_pty_uniform = false;
  if (b->has[CDO_DIFFUSION]) {

    bool is_uniform = cs_property_is_uniform(eqp->diffusion_property);

    b->diff_pty_uniform = is_uniform;
    b->diff = cs_cdo_diffusion_builder_init(connect,
                                            CS_SPACE_SCHEME_CDOVB,
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

    BFT_MALLOC(b->source_terms, b->n_vertices, cs_real_t);

    for (int  st_id = 0; st_id < eqp->n_source_terms; st_id++) {
      const cs_source_term_t  *st = eqp->source_terms[st_id];
      if (cs_source_term_get_reduction(st) == CS_SOURCE_TERM_REDUC_PRIM)
        b->flag |= CS_CDO_PRIMAL_SOURCE | CS_CDO_BUILD_HCONF;
      else
        b->flag |= CS_CDO_DUAL_SOURCE;
    }

  } /* There is at least one source term */

  if (b->flag & CS_CDO_BUILD_HCONF ||
      b->flag & CS_CDO_BUILD_LOC_HCONF) {

    cs_param_hodge_t  hwbs_info = {.inv_pty = false,
                                   .type = CS_PARAM_HODGE_TYPE_VPCD,
                                   .algo = CS_PARAM_HODGE_ALGO_WBS,
                                   .coef = 1.0}; // not useful in this case

    b->hb = cs_hodge_builder_init(connect, hwbs_info);

    if ((b->flag & CS_CDO_BUILD_HCONF) && cs_cdovb_hconf == NULL)
      _build_hvpcd_conf(b);

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

  /* Strong enforcement means that we need an indirection list between the
     compress (or zip) and initial numbering of vertices */

  b->v_z2i_ids = NULL; // zipped --> initial ids
  b->v_i2z_ids = NULL; // initial --> zipped ids

  switch (b->enforce) {
  case CS_PARAM_BC_ENFORCE_STRONG:
    if (b->vtx_dir->n_elts > 0) {

      if (b->has[CDO_ADVECTION] || b->has[CDO_TIME])
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid choice of enforcement of the boundary conditions.\n"
                  " Strong enforcement is not implemented when convection or"
                  " unsteady terms are activated.\n"
                  " Please modify your settings.");

      bool  *is_kept = NULL;

      b->n_dof_vertices = n_vertices - b->vtx_dir->n_elts;

      /* Build b->v_z2i_ids and b->i2i_ids */
      BFT_MALLOC(b->v_z2i_ids, b->n_dof_vertices, cs_lnum_t);
      BFT_MALLOC(b->v_i2z_ids, b->n_vertices, cs_lnum_t);
      BFT_MALLOC(is_kept, n_vertices, bool);

# pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_vertices; i++)
        is_kept[i] = true, b->v_i2z_ids[i] = -1; // by default, set to remove

# pragma omp parallel for if (b->vtx_dir->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->vtx_dir->n_elts; i++)
        is_kept[b->vtx_dir->elt_ids[i]] = false;

      cs_lnum_t  cur_id = 0;
      for (cs_lnum_t i = 0; i < b->n_vertices; i++) {
        if (is_kept[i]) {
          b->v_i2z_ids[i] = cur_id;
          b->v_z2i_ids[cur_id++] = i;
        }
      }
      assert(cur_id == b->n_dof_vertices);

      BFT_FREE(is_kept);

    }
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    if (b->has[CDO_DIFFUSION])
      cs_cdo_diffusion_build_c2bcbf(connect,
                                    b->face_bc->dir,
                                    &(b->c2bcbf_idx),
                                    &(b->c2bcbf_ids));
    break;

  default: // Nothing to do
    break;

  } /* Strong enforcement of BCs */

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovb_scaleq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdovb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_scaleq_free(void   *builder)
{
  cs_cdovb_scaleq_t  *b = (cs_cdovb_scaleq_t *)builder;

  if (b == NULL)
    return b;

  /* eqp is only shared. Thies structure is freed later. */

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

  /* Renumbering (if strong enforcement of BCs for instance) */
  if (b->n_vertices > b->n_dof_vertices) {
    BFT_FREE(b->v_z2i_ids);
    BFT_FREE(b->v_i2z_ids);
  }

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
cs_cdovb_scaleq_free_sysmat(void              *builder,
                            cs_sla_matrix_t   *matrix)
{
  CS_UNUSED(builder);

  /* Free matrix */
  matrix = cs_sla_matrix_free(matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in, out]  builder     pointer to a cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_compute_source(void   *builder)
{
  cs_desc_t  desc;

  double  *work = cs_equation_get_tmpbuf();
  double  *st_eval = NULL, *primal_cumul = NULL;
  cs_cdovb_scaleq_t  *b = (cs_cdovb_scaleq_t *)builder;
  const cs_equation_param_t  *eqp = b->eqp;

  if (eqp->n_source_terms == 0)
    return;

  st_eval = work;

  /* Reset source term array */
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->n_vertices; i++)
    b->source_terms[i] = 0;

  if (b->flag & CS_CDO_PRIMAL_SOURCE) {

    primal_cumul = work + b->n_vertices;

# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < b->n_vertices; i++)
      primal_cumul[i] = 0;
  }

  /* Loop on source term definitions */
  for (int st_id = 0; st_id < eqp->n_source_terms; st_id++) {

    const cs_source_term_t  *st = eqp->source_terms[st_id];

    if (cs_source_term_get_reduction(st) == CS_SOURCE_TERM_REDUC_DUAL) {

      desc.location = CS_FLAG_SCAL | cs_cdo_dual_cell;
      desc.state = CS_FLAG_STATE_DENSITY;

      /* st_eval is reset and updated inside this function */
      cs_source_term_compute(desc, st, &st_eval);

      /* Update source term array */
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_vertices; i++)
        b->source_terms[i] += st_eval[i];

    }
    else { // Source term is defined thanks to a reduction on primal entities

      assert(cs_source_term_get_reduction(st) == CS_SOURCE_TERM_REDUC_PRIM);

      desc.location = CS_FLAG_SCAL | cs_cdo_primal_vtx;
      desc.state = CS_FLAG_STATE_POTENTIAL;

      /* st_eval is reset and updated inside this function */
      cs_source_term_compute(desc, st, &st_eval);

# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_vertices; i++)
        primal_cumul[i] += st_eval[i];

    }

  } // Loop on source terms

  /* Last step for primal source terms */
  if (b->flag & CS_CDO_PRIMAL_SOURCE) {

    assert(cs_cdovb_hconf != NULL);

    /* Update source term array */
    cs_sla_matvec(cs_cdovb_hconf, primal_cumul, &(b->source_terms), false);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO vertex-based scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to cs_cdovb_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] sla_mat    pointer to cs_sla_matrix_t structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_build_system(const cs_mesh_t             *mesh,
                             const cs_real_t             *field_val,
                             double                       dt_cur,
                             void                        *builder,
                             cs_real_t                  **rhs,
                             cs_sla_matrix_t            **sla_mat)
{
  cs_cdovb_scaleq_t  *b = (cs_cdovb_scaleq_t *)builder;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_equation_param_t  *eqp = b->eqp;

  /* Default flag value for vertex-based scalar equations */
  cs_flag_t  cm_flag = cs_cdovb_cmflag;

  /* Allocate and initialize a matrix with the larger stencil (that related
     to diffusion => all vertices of a cell are potentially in interaction)
     Pattern is shared between all CDO equations using vertices as DoFs thanks
     to connect->v2v.
  */
  cs_sla_matrix_t  *sys_mat =
    cs_sla_matrix_create_msr_from_index(connect->v2v,
                                        false,  // symmetric
                                        true,   // sorted
                                        1);     // stride

  if (!b->has[CDO_ADVECTION] && b->enforce != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
    sys_mat->flag |= CS_SLA_MATRIX_SYM;

  /* Preparatory step for diffusion term */
  cs_real_33_t  diff_tensor = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  if (b->has[CDO_DIFFUSION]) {

    cs_hodge_builder_t  *hb = cs_cdo_diffusion_get_hodge_builder(b->diff);

    cs_hodge_builder_unset(hb);
    if (b->diff_pty_uniform)
      cs_property_get_cell_tensor(0, // cell_id
                                  eqp->diffusion_property,
                                  eqp->diffusion_hodge.inv_pty,
                                  diff_tensor);

    if (eqp->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_WBS ||
        b->enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
        b->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM)
      cm_flag |= CS_CDO_LOCAL_F | CS_CDO_LOCAL_FE;

  }

  /* Preparatory step for advection term */
  if (b->has[CDO_ADVECTION])
    cm_flag |=  CS_CDO_LOCAL_F | CS_CDO_LOCAL_FE;

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
    BFT_MALLOC(full_rhs, b->n_vertices, cs_real_t);

# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->n_vertices; i++)
    full_rhs[i] = 0.0;

  /* Add the contribution of source terms to the full rhs for this time step */
  _add_source_terms(b, full_rhs);

  /* Compute the values of the Dirichlet BC.
     TODO: do the analogy for Neumann BC */
  _compute_dir_values(mesh, field_val, b);

  /* Temporary pre-allocated buffers (the following buffer is used in
     _add_source_terms thus be careful if the order of calls is changed) */
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

  /* Main loop on cells to build the linear system */
  /* --------------------------------------------- */

  /* Prepare evolution for a full openMP implementation
     Get the cell-wise view of the mesh and the algebraic system */
  cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(0);

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Set the local mesh structure for the current cell */
    cs_cell_mesh_build(c_id, cm_flag, connect, quant, cm);

    /* Cell-wise view of the linear system to build */
    const int  n_vc = cm->n_vc;
    cs_locmat_t  *hconf_c = NULL;

    /* Store the local values attached to Dirichlet values if the current cell
       has at least one border face */
    if (cell_flag[c_id] & CS_CDO_CONNECT_BD)
      for (short int v = 0; v < n_vc; v++)
        cs_cell_sys->dir_bc[v] = dir_bc_vals[cm->v_ids[v]];

    /* Initialize the local system */
    cs_cell_sys->mat->n_ent = n_vc;
    for (short int v = 0; v < n_vc; v++) {
      cs_cell_sys->mat->ids[v] = cm->v_ids[v];
      cs_cell_sys->rhs[v] = 0.;
    }
    for (short int v = 0; v < n_vc*n_vc; v++)
      cs_cell_sys->mat->val[v] = 0;

    /* DIFFUSION TERM */
    if (b->has[CDO_DIFFUSION]) { /* Define the local stiffness matrix */

      if (b->diff_pty_uniform == false)
        cs_property_get_cell_tensor(c_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    diff_tensor);

      cs_locmat_t  *diff_mat = // local matrix owned by the diffusion builder
        cs_cdo_diffusion_build_local(quant,
                                     cm,
              (const cs_real_3_t (*))diff_tensor,
                                     b->diff);

      cs_locmat_add(cs_cell_sys->mat, diff_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      bft_printf(">> Local diffusion matrix");
      cs_locmat_dump(c_id, diff_mat);
#endif

      /* Weakly enforced Dirichlet BCs for cells attached to the boundary */
      if (b->c2bcbf_idx != NULL && cell_flag[c_id] & CS_CDO_CONNECT_BD) {

        for (cs_lnum_t j = b->c2bcbf_idx[c_id];
             j < b->c2bcbf_idx[c_id+1]; j++) {

          /* cs_cell_sys is updated inside (matrix and rhs) */
          cs_cdo_diffusion_weak_bc(b->c2bcbf_ids[j], // border face id
                                   cm,
                                   (const cs_real_3_t (*))diff_tensor,
                                   b->diff,
                                   cs_cell_sys);

        } // Loop on border faces

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
      bft_printf(">> Local diffusion matrix after weak enforcement");
      cs_locmat_dump(c_id, cs_cell_sys->mat);
#endif
      } /* Weak enforcement of Dirichlets BCs */

    } /* DIFFUSION */

    /* ADVECTION TERM */
    if (b->has[CDO_ADVECTION]) { /* Define the local advection matrix */

      cs_locmat_t  *adv_mat =
        cs_cdovb_advection_build(cm, eqp,
                                 (const cs_real_3_t (*))diff_tensor,
                                 b->adv);

      cs_locmat_add(cs_cell_sys->mat, adv_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      bft_printf(">> Local advection matrix");
      cs_locmat_dump(c_id, adv_mat);
#endif

      /* Last treatment for the advection term: Apply boundary conditions */
      if (cell_flag[c_id] & CS_CDO_CONNECT_BD) {

        /* cs_cell_sys is updated inside (matrix and rhs) */
        cs_cdovb_advection_add_bc(cm, eqp, b->adv, cs_cell_sys);

      } // Apply BC

    } /* ADVECTION */

    if (b->flag & CS_CDO_BUILD_LOC_HCONF)
      hconf_c = cs_hodge_build_cellwise(cm, b->hb);

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
      cs_locmat_mult_add(cs_cell_sys->mat, rpty_val, hconf_c);

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
         Update cs_cell_sys (matrix and rhs) */
      _apply_time_scheme(tpty_val, field_val, cm, hconf_c, b, cs_cell_sys);

    } /* Time contribution */

    /* Assemble the matrix related to the advcetion/diffusion/reaction terms
       If advection is activated, the resulting system is not symmetric
       Otherwise, the system is symmetric with extra-diagonal terms. */
    if (sys_mat->flag & CS_SLA_MATRIX_SYM)
      cs_sla_assemble_msr_sym(cs_cell_sys->mat, sys_mat, false);
    else
      cs_sla_assemble_msr(cs_cell_sys->mat, sys_mat);

    /* Assemble the right-hand side (rhs) */
    for (short int v = 0; v < n_vc; v++)
      full_rhs[cm->v_ids[v]] += cs_cell_sys->rhs[v];

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
    bft_printf(">> (FINAL) Local system matrix");
    cs_locmat_dump(c_id, cs_cell_sys->mat);
#endif

  } // Main loop on cells

  /* Final step in BC management.
     Apply the strong or penalized enforcement. In case of Nitsche enforcement,
     there is nothing to do (already done).
     Must be call after the application of the time scheme */
  _enforce_bc(b, &full_rhs, &sys_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
  cs_sla_system_dump("system-BeforeClean.log", NULL, sys_mat, full_rhs);
#endif

  /* Clean matrix (set entries to zero if there are below a given threshold
     in order to improve the matrix conditionning */
  cs_sla_matrix_clean(eqp->verbosity, cs_cdovb_threshold, sys_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
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
 * \param[in, out] builder    pointer to cs_cdovb_scaleq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_update_field(const cs_real_t     *solu,
                             const cs_real_t     *rhs,
                             void                *builder,
                             cs_real_t           *field_val)
{
  CS_UNUSED(rhs);

  cs_cdovb_scaleq_t  *b = (cs_cdovb_scaleq_t  *)builder;

  /* Set computed solution in field array */
  if (b->n_dof_vertices < b->n_vertices) {

    const cs_cdo_bc_list_t  *v_dir = b->vtx_dir;

    /* Sanity check */
    assert(b->enforce == CS_PARAM_BC_ENFORCE_STRONG);

# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < b->n_vertices; i++)
      field_val[i] = 0; // To tag unset values

# pragma omp parallel for if (b->n_dof_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < b->n_dof_vertices; i++)
      field_val[b->v_z2i_ids[i]] = solu[i];

    for (cs_lnum_t i = 0; i < v_dir->n_nhmg_elts; i++)
      field_val[v_dir->elt_ids[i]] = b->dir_val[i];

  }
  else
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < b->n_vertices; i++)
      field_val[i] = solu[i];

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
cs_cdovb_scaleq_compute_flux_across_plane(const void          *builder,
                                          const cs_real_t     *pdi,
                                          int                  ml_id,
                                          const cs_real_t      direction[],
                                          double              *diff_flux,
                                          double              *conv_flux)
{
  const cs_cdovb_scaleq_t  *b = (const cs_cdovb_scaleq_t  *)builder;
  const cs_equation_param_t  *eqp = b->eqp;

  cs_mesh_location_type_t  ml_t = cs_mesh_location_get_type(ml_id);

  *diff_flux = 0.;
  *conv_flux = 0.;

  if (pdi == NULL)
    return;

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

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_sla_matrix_t  *f2c = connect->f2c;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  double  pf;
  cs_real_3_t  gc, pty_gc;
  cs_real_33_t  pty_tens;
  cs_nvec3_t  adv_c;

  if (ml_t == CS_MESH_LOCATION_BOUNDARY_FACES) {

    const cs_lnum_t  n_i_faces = connect->f_info->n_i_elts;
    const cs_lnum_t  shift_if = 2*n_i_faces;

    for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

      const cs_lnum_t  bf_id = elt_ids[i];
      const cs_lnum_t  f_id = n_i_faces + bf_id;
      const cs_lnum_t  c_id = f2c->col_id[shift_if + bf_id];
      const cs_quant_t  f = quant->face[f_id];
      const short int  sgn = (_dp3(f.unitv, direction) < 0) ? -1 : 1;
      const double  coef = sgn * f.meas;

      if (b->has[CDO_DIFFUSION]) { /* Compute the local diffusive flux */

        cs_reco_grd_cell_from_pv(c_id, connect, quant, pdi, gc);
        cs_property_get_cell_tensor(c_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);
        cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, gc, pty_gc);
        *diff_flux += -coef * _dp3(f.unitv, pty_gc);

      }

      if (b->has[CDO_ADVECTION]) { /* Compute the local advective flux */

        cs_advection_field_get_cell_vector(c_id, eqp->advection_field, &adv_c);
        cs_reco_pf_from_pv(f_id, connect, quant, pdi, &pf);
        *conv_flux += coef * adv_c.meas * _dp3(adv_c.unitv, f.unitv) * pf;

      }

    } // Loop on selected border faces

  }
  else if (ml_t == CS_MESH_LOCATION_INTERIOR_FACES) {

    for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

      const cs_lnum_t  f_id = elt_ids[i];
      const cs_lnum_t  shift_f = 2*f_id;
      const cs_lnum_t  c1_id = f2c->col_id[shift_f];
      const cs_lnum_t  c2_id = f2c->col_id[shift_f+1];
      const cs_quant_t  f = quant->face[f_id];
      const short int  sgn = (_dp3(f.unitv, direction) < 0) ? -1 : 1;
      const double  coef = 0.5 * sgn * f.meas;

      if (b->has[CDO_DIFFUSION]) { /* Compute the local diffusive flux */

        cs_reco_grd_cell_from_pv(c1_id, connect, quant, pdi, gc);
        cs_property_get_cell_tensor(c1_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);
        cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, gc, pty_gc);
        *diff_flux += -coef * _dp3(f.unitv, pty_gc);

        cs_reco_grd_cell_from_pv(c2_id, connect, quant, pdi, gc);
        cs_property_get_cell_tensor(c2_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);
        cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, gc, pty_gc);
        *diff_flux += -coef * _dp3(f.unitv, pty_gc);

      }

      if (b->has[CDO_ADVECTION]) { /* Compute the local advective flux */

        cs_reco_pf_from_pv(f_id, connect, quant, pdi, &pf);

        cs_advection_field_get_cell_vector(c1_id, eqp->advection_field, &adv_c);
        *conv_flux += coef * adv_c.meas * _dp3(adv_c.unitv, f.unitv) * pf;

        cs_advection_field_get_cell_vector(c2_id, eqp->advection_field, &adv_c);
        *conv_flux += coef * adv_c.meas * _dp3(adv_c.unitv, f.unitv) * pf;

      }

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
cs_cdovb_scaleq_cellwise_diff_flux(const cs_real_t   *values,
                                   void              *builder,
                                   cs_real_t         *diff_flux)
{
  cs_cdovb_scaleq_t  *b = (cs_cdovb_scaleq_t  *)builder;

  // To be modified for an fully integration of openMP
  cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(0);

  const cs_equation_param_t  *eqp = b->eqp;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_connect_index_t  *c2e = connect->c2e;

  if (!b->has[CDO_DIFFUSION]) { // No diffusion

    for (int i = 0; i < c2e->idx[quant->n_cells]; i++)
      diff_flux = 0;

    return;
  }

  /* Sanity check */
  assert(cm->n_max_ebyc == connect->n_max_ebyc);

  /* Retrieve temporary buffers */
  double  *work = cs_equation_get_tmpbuf();
  double  *p_v = work; // used as a temporary buffer

  /* Default flag value for vertex-based scalar equations */
  cs_flag_t  cm_flag = cs_cdovb_cmflag;

  /* Diffusion tensor */
  cs_hodge_builder_t  *hbd = cs_cdo_diffusion_get_hodge_builder(b->diff);
  cs_real_33_t  diff_tensor = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  cs_hodge_builder_unset(hbd);
  if (b->diff_pty_uniform)
    cs_property_get_cell_tensor(0,  // cell_id
                                eqp->diffusion_property,
                                eqp->diffusion_hodge.inv_pty,
                                diff_tensor);

  if (eqp->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_WBS)
    cm_flag |= CS_CDO_LOCAL_F | CS_CDO_LOCAL_FE;

  /* Define the flux by cellwise contributions
     loc_flux = - loc_hodge * loc_gradient(h) */

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    cs_real_t  *flux = diff_flux + c2e->idx[c_id];

    /* Set the local mesh structure for the current cell */
    cs_cell_mesh_build(c_id, cm_flag, connect, quant, cm);

    /* Define a local buffer keeping the value of the discrete potential
       for the current cell */
    for (short int v = 0; v < cm->n_vc; v++)
      p_v[v] = values[cm->v_ids[v]];

    if (b->diff_pty_uniform == false)
      cs_property_get_cell_tensor(c_id,
                                  eqp->diffusion_property,
                                  eqp->diffusion_hodge.inv_pty,
                                  diff_tensor);

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
    case CS_PARAM_HODGE_ALGO_VORONOI:
      {
        double  *vec = work + cm->n_vc;

        assert(cs_equation_get_tmpbuf_size() >=
               (size_t)cm->n_vc + (size_t)cm->n_ec); // Sanity check

        if (b->diff_pty_uniform == false ||
            cs_hodge_builder_get_setting_flag(hbd) == false)
          cs_hodge_builder_set_tensor(hbd, (const cs_real_t (*)[3])diff_tensor);

        /* Build the local dense matrix related to this operator */
        const cs_locmat_t  *hloc = cs_hodge_build_cellwise(cm, hbd);

        for (short int e = 0; e < cm->n_ec; e++) {

          const short int  sgn_v1 = cm->e2v_sgn[2*e]; // sgn_v2 = -sgn_v1

          /* Used this buffer to store (temporary) the values of the local
             discrete gradient. Then, flux = - Hloc * grd_c(pdi_c) */
          vec[e] = sgn_v1 * (p_v[cm->e2v_ids[2*e+1]] - p_v[cm->e2v_ids[2*e]]);

        } // Loop on cell edges

        /* Store the local fluxes into diff_flux */
        cs_locmat_matvec(hloc, vec, flux);

      }
      break;

    case CS_PARAM_HODGE_ALGO_WBS:
      {
        /* Compute reconstructed values at cell center: p_c */
        double  p_c = 0;
        for (short int v = 0; v < cm->n_vc; v++)
          p_c += cm->wvc[v]*p_v[v];

        /* Compute the flux across dual faces for this cell */
        cs_cdo_diffusion_cellwise_flux(cm,
                                       quant->dface + c2e->idx[c_id],
                                       (const cs_real_3_t (*))diff_tensor,
                                       p_v,
                                       p_c,
                                       b->diff,
                                       diff_flux + c2e->idx[c_id]);

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid algorithm for computing the diffusive flux.");
      break;

    } // switch on algo.

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
cs_cdovb_scaleq_extra_op(const char            *eqname,
                         const cs_field_t      *field,
                         void                  *builder)
{
  CS_UNUSED(field);

  cs_cdovb_scaleq_t  *b = (cs_cdovb_scaleq_t  *)builder;

  const cs_equation_param_t  *eqp = b->eqp;

  if (b->has[CDO_ADVECTION] &&
      (eqp->process_flag & CS_EQUATION_POST_UPWIND_COEF)) {

    cs_real_t  *work_c = cs_equation_get_tmpbuf();
    char *postlabel = NULL;
    int  len = strlen(eqname) + 8 + 1;

    BFT_MALLOC(postlabel, len, char);
    sprintf(postlabel, "%s.UpwCoef", eqname);

    /* Compute in each cell an evaluation of upwind weight value */
    cs_cdo_advection_get_upwind_coef_cell(cs_shared_quant,
                                          eqp->advection_info,
                                          work_c);

    cs_post_write_var(-1,                   // id du maillage de post
                      postlabel,
                      1,
                      true,                 // interlace
                      true,                 // true = original mesh
                      CS_POST_TYPE_cs_real_t,
                      work_c,               // values on cells
                      NULL,                 // values at internal faces
                      NULL,                 // values at border faces
                      cs_shared_time_step); // time step management struct.

    BFT_FREE(postlabel);

  } // Post a Peclet attached to cells

}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
