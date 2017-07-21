/*============================================================================
 * Build an algebraic CDO face-based system for unsteady convection/diffusion
 * reaction of scalar equations with source terms
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
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
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_SCALEQ_DBG      4
#define CS_CDOFB_SCALEQ_MODULO  10

/* Algebraic system for CDO face-based discretization */

struct  _cs_cdofb_scaleq_t {

  /* Pointer to a cs_equation_param_t structure shared with a cs_equation_t
     structure.  */

  const cs_equation_param_t  *eqp;

  /* System size */
  cs_lnum_t    n_dofs;        // n_faces + n_cells

  /* Shortcut to know what to build */
  cs_flag_t    msh_flag;     // Information related to cell mesh
  cs_flag_t    bd_msh_flag;  // Information related to cell mesh (boundary)
  cs_flag_t    st_msh_flag;  // Information related to cell mesh (source term)
  cs_flag_t    sys_flag;     // Information related to the sytem

  /* Solution of the algebraic system at the last iteration */
  cs_real_t   *face_values;  /* DoF unknowns (x) + BCs */

  cs_real_t   *cell_rhs;  // right-hand side related to cell dofs
  cs_real_t   *acc_inv;   // inverse of a diagonal matrix (block cell-cell)
  cs_real_t   *acf;       /* Lower-Left block of the full matrix
                             (block cell-vertices). Access to the values thanks
                             to the c2f connectivity */

  /* Metadata related to associated properties */
  bool         diff_pty_uniform;
  bool         time_pty_uniform;
  bool         reac_pty_uniform[CS_CDO_N_MAX_REACTIONS];

  /* Source terms */
  cs_real_t    *source_terms; /* Array storing the value arising from the
                                 contribution of all source terms */
  cs_mask_t    *source_mask;  /* NULL if at least one source term is not
                                 defined for all cells (size = n_cells) */

  /* Pointer to functions which compute the value of the source term */
  cs_source_term_cellwise_t  *compute_source[CS_N_MAX_SOURCE_TERMS];

  /* Boundary conditions:

     face_bc should not change during the simulation.
     The case of a definition of the BCs which changes of type during the
     simulation is possible but not implemented.
     You just have to call the initialization step each time the type of BCs
     is modified to define an updated cs_cdo_bc_t structure.

     We translate Dirichlet BCs to the center of border faces.
     The values can be modified at each time step in case of transient
     simulation.
     For Neumann and Robin BCs,
     BCs contributions are computed on the fly.

   */

  cs_cdo_bc_t           *face_bc; // list of faces sorted by type of BCs

  /* Pointer of function to build the diffusion term */
  cs_hodge_t                      *get_stiffness_matrix;
  cs_hodge_t                      *get_diffusion_hodge;
  cs_cdo_diffusion_enforce_dir_t  *enforce_dirichlet;
  cs_cdo_diffusion_flux_trace_t   *boundary_flux_op;

  /* Pointer of function to build the advection term */
  cs_cdo_advection_t              *get_advection_matrix;
  cs_cdo_advection_bc_t           *add_advection_bc;

  /* Pointer of function to apply the time scheme */
  cs_cdo_time_scheme_t            *apply_time_scheme;

  /* Monitoring the efficiency */
  cs_equation_monitor_t           *monitor;

};

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */
static cs_cell_sys_t      **cs_cdofb_cell_sys = NULL;
static cs_cell_bc_t       **cs_cdofb_cell_bc = NULL;
static cs_cell_builder_t  **cs_cdofb_cell_bld = NULL;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_time_step_t  *cs_shared_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Retrieve the flag to give for building a cs_cell_mesh_t structure
 *
 * \param[in]      cell_flag   flag related to the current cell
 * \param[in]      v_msh_flag  default mesh flag for the volumic terms
 * \param[in]      b_msh_flag  default mesh flag for the boundary terms
 * \param[in]      s_msh_flag  default mesh flag for the source terms
 *
 * \return the flag to set for the current cell
 */
/*----------------------------------------------------------------------------*/

static inline cs_flag_t
_get_cell_mesh_flag(cs_flag_t       cell_flag,
                    cs_flag_t       v_msh_flag,
                    cs_flag_t       b_msh_flag,
                    cs_flag_t       s_msh_flag)
{
  cs_flag_t  msh_flag = v_msh_flag | s_msh_flag;

  if (cell_flag & CS_FLAG_BOUNDARY)
    msh_flag |= b_msh_flag;

  return msh_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *
 * \param[in]      cell_flag   flag related to the current cell
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      b           pointer to a cs_cdofb_scaleq_t structure
 * \param[in]      dir_values  Dirichlet values associated to each face
 * \param[in]      neu_tags    Definition id related to each Neumann face
 * \param[in]      field_tn    values of the field at the last computed time
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cbc         pointer to a cellwise view of the BCs
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_init_cell_structures(const cs_flag_t             cell_flag,
                      const cs_cell_mesh_t       *cm,
                      const cs_cdofb_scaleq_t    *b,
                      const cs_real_t             dir_values[],
                      const short int             neu_tags[],
                      const cs_real_t             field_tn[],
                      cs_cell_sys_t              *csys,
                      cs_cell_bc_t               *cbc,
                      cs_cell_builder_t          *cb)
{
  CS_UNUSED(cb);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_bc_t  *face_bc = b->face_bc;

  /* Cell-wise view of the linear system to build */
  const int  n_fc = cm->n_fc;
  const int  n_dofs = n_fc + 1;

  /* Initialize the local system */
  csys->c_id = cm->c_id;
  csys->n_dofs = n_dofs;
  csys->mat->n_ent = n_dofs;
  for (short int f = 0; f < n_fc; f++) {
    csys->mat->ids[f] = cm->f_ids[f];
    csys->val_n[f] = b->face_values[cm->f_ids[f]];
    csys->rhs[f] = csys->source[f] = 0.;
  }
  csys->mat->ids[n_fc] = cm->c_id;
  csys->val_n[n_fc] = field_tn[cm->c_id];
  csys->rhs[n_fc] = b->cell_rhs[cm->c_id]; // Used for static condensation
  csys->source[n_fc] = 0;

  for (short int i = 0; i < n_dofs*n_dofs; i++) csys->mat->val[i] = 0;

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */
  if (cell_flag & CS_FLAG_BOUNDARY) {

    /* Reset values */
    cbc->n_bc_faces = 0;
    cbc->n_dirichlet = cbc->n_nhmg_neuman = cbc->n_robin = 0;
    cbc->n_dofs = n_fc;
    for (short int i = 0; i < n_fc; i++) {
      cbc->dof_flag[i] = 0;
      cbc->dir_values[i] = cbc->neu_values[i] = 0;
      cbc->rob_values[2*i] = cbc->rob_values[2*i+1] = 0.;
    }

    /* Identify which face is a boundary face */
    for (short int f = 0; f < cm->n_fc; f++) {

      const cs_lnum_t  f_id = cm->f_ids[f] - connect->n_faces[2]; // n_i_faces
      if (f_id > -1) { // Border face

        cs_flag_t  face_flag = face_bc->flag[f_id];

        cbc->face_flag[cbc->n_bc_faces] = face_flag;
        cbc->bf_ids[cbc->n_bc_faces++] = f;

        if (face_flag & CS_CDO_BC_HMG_DIRICHLET)
            cbc->dof_flag[f] |= CS_CDO_BC_HMG_DIRICHLET;
        else if (face_flag & CS_CDO_BC_DIRICHLET) {
          cbc->dir_values[f] = dir_values[f_id];
          cbc->dof_flag[f] |= CS_CDO_BC_DIRICHLET;
        }
        else if (face_flag & CS_CDO_BC_NEUMANN) {

          cbc->dof_flag[f] |= CS_CDO_BC_NEUMANN;
          cs_equation_compute_neumann_sf(neu_tags[f_id], f, b->eqp, cm, cbc);

        }

      } // Border faces

    } // Loop on cell faces

    /* Update counters. Only DoFs related to faces are possibly attached
       to a boundary condition */
    for (short int f = 0; f < n_fc; f++) {

      if (cbc->dof_flag[f] & CS_CDO_BC_HMG_DIRICHLET ||
          cbc->dof_flag[f] & CS_CDO_BC_DIRICHLET)
        cbc->n_dirichlet += 1;

      if (cbc->dof_flag[f] & CS_CDO_BC_NEUMANN)
        cbc->n_nhmg_neuman += 1;

      if (cbc->dof_flag[f] & CS_CDO_BC_ROBIN)
        cbc->n_robin += 1;

    } // Loop on cell faces

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    for (short int f = 0; f < n_fc; f++) {
      if (cbc->dof_flag[f] & CS_CDO_BC_HMG_DIRICHLET)
        if (fabs(cbc->dir_values[f]) > 10*DBL_MIN)
          bft_error(__FILE__, __LINE__, 0,
                    "Invalid enforcement of Dirichlet BCs on faces");
    }
#endif

  } /* Border cell */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Proceed to a static condensation of the local system and keep
 *          information inside the builder to be able to compute the values
 *          at cell centers
 *
 * \param[in]      c2f         pointer to a cs_adjacency_t structure
 * \param[in, out] builder     pointer to a cs_cdovcb_scaleq_t structure
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 * \param[in, out] csys        pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_condense_and_store(const cs_adjacency_t    *c2f,
                    cs_cdofb_scaleq_t       *b,
                    cs_cell_builder_t       *cb,
                    cs_cell_sys_t           *csys)
{
  const int  n_dofs = csys->n_dofs;
  const int  n_fc = n_dofs - 1;

  /* Store information to compute the cell values after the resolution */
  const double  *row_c = csys->mat->val + n_dofs*n_fc; // Last row
  const double  inv_acc = 1/row_c[n_fc];
  const double  cell_rhs = csys->rhs[n_fc];

  b->acc_inv[csys->c_id] = inv_acc;
  b->cell_rhs[csys->c_id] = cell_rhs;

  /* Store the cell row (the last one) */
  double  *acf = b->acf + c2f->idx[csys->c_id];
  for (int f = 0; f < n_fc; f++) acf[f] = row_c[f];

  /* Store the cell column temporary (the last one) */
  double  *afc = cb->values;
  for (int f = 0; f < n_fc; f++) afc[f] = csys->mat->val[n_dofs*f + n_fc];

  /* Update csys */
  csys->n_dofs = n_fc;
  csys->mat->n_ent = n_fc;
  for (short int i = 0; i < n_fc; i++) {

    double  *old_i = csys->mat->val + n_dofs*i; // old row_i
    double  *new_i = csys->mat->val + n_fc*i;   // new row_i

    /* Condensate the local matrix */
    const double  coef_i = inv_acc * afc[i];
    for (short int j = 0; j < n_fc; j++)
      new_i[j] = old_i[j] - acf[j]*coef_i;

    /* Update RHS: RHS_f = RHS_f - Afc*Acc^-1*s_c */
    csys->rhs[i] -= cell_rhs*coef_i;

  } // Loop on fi cell faces

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
cs_cdofb_scaleq_set_shared_pointers(const cs_cdo_quantities_t    *quant,
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
 *         face-based schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_initialize(void)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Specific treatment for handling openMP */
  BFT_MALLOC(cs_cdofb_cell_sys, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(cs_cdofb_cell_bc, cs_glob_n_threads, cs_cell_bc_t *);
  BFT_MALLOC(cs_cdofb_cell_bld, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_cdofb_cell_sys[i] = NULL;
    cs_cdofb_cell_bc[i] = NULL;
    cs_cdofb_cell_bld[i] = NULL;
  }

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cdofb_cell_sys[t_id] = cs_cell_sys_create(connect->n_max_fbyc + 1);
    cs_cdofb_cell_bc[t_id] = cs_cell_bc_create(connect->n_max_fbyc + 1,
                                               connect->n_max_fbyc);
    cs_cdofb_cell_bld[t_id] = cs_cell_builder_create(CS_SPACE_SCHEME_CDOFB,
                                                     connect);
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cdofb_cell_sys[0] = cs_cell_sys_create(connect->n_max_fbyc + 1);
  cs_cdofb_cell_bc[0] = cs_cell_bc_create(connect->n_max_fbyc + 1,
                                           connect->n_max_fbyc);
  cs_cdofb_cell_bld[0] = cs_cell_builder_create(CS_SPACE_SCHEME_CDOFB,
                                                connect);
#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO face-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_finalize(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_cell_sys_free(&(cs_cdofb_cell_sys[t_id]));
    cs_cell_bc_free(&(cs_cdofb_cell_bc[t_id]));
    cs_cell_builder_free(&(cs_cdofb_cell_bld[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(cs_cdofb_cell_sys[0]));
  cs_cell_bc_free(&(cs_cdofb_cell_bc[0]));
  cs_cell_builder_free(&(cs_cdofb_cell_bld[0]));
#endif /* openMP */

  BFT_FREE(cs_cdofb_cell_sys);
  BFT_FREE(cs_cdofb_cell_bc);
  BFT_FREE(cs_cdofb_cell_bld);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_scaleq_t structure
 *
 * \param[in]  eqp        pointer to a cs_equation_param_t structure
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_scaleq_init(const cs_equation_param_t   *eqp,
                     const cs_mesh_t             *mesh)
{
  CS_UNUSED(mesh);
  /* Sanity checks */
  assert(eqp != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOFB && eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Invalid type of equation.\n"
              " Expected: scalar-valued CDO face-based equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_faces = connect->n_faces[0];
  const cs_lnum_t  n_b_faces = connect->n_faces[1];

  cs_cdofb_scaleq_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_cdofb_scaleq_t);

  /* Shared pointers */
  b->eqp = eqp;

  /* Dimensions of the algebraic system */
  b->n_dofs = n_faces + n_cells;

  /* Store a direct access to which term one has to compute */
  b->sys_flag = 0;
  if (eqp->flag & CS_EQUATION_DIFFUSION)
    b->sys_flag |= CS_FLAG_SYS_DIFFUSION;
  if (eqp->flag & CS_EQUATION_CONVECTION)
    b->sys_flag |= CS_FLAG_SYS_ADVECTION;
  if (eqp->flag & CS_EQUATION_REACTION)
    b->sys_flag |= CS_FLAG_SYS_REACTION;
  if (eqp->flag & CS_EQUATION_UNSTEADY)
    b->sys_flag |= CS_FLAG_SYS_TIME;
  if (eqp->n_source_terms > 0)
    b->sys_flag |= CS_FLAG_SYS_SOURCETERM;

  /* Flag to indicate what to build in a cell mesh */
  b->msh_flag = CS_CDO_LOCAL_PF | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  b->bd_msh_flag = 0;
  for (int i = 0; i < eqp->n_bc_desc; i++) {
    const cs_xdef_t  *def = eqp->bc_desc[i];
    if (def->meta & CS_CDO_BC_NEUMANN)
      if (def->qtype == CS_QUADRATURE_BARY_SUBDIV ||
          def->qtype == CS_QUADRATURE_HIGHER ||
          def->qtype == CS_QUADRATURE_HIGHEST)
        b->bd_msh_flag |= CS_CDO_LOCAL_EV | CS_CDO_LOCAL_EF | CS_CDO_LOCAL_EFQ;
  }

  /* Set members and structures related to the management of the BCs
     Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.  */
  b->face_bc = cs_cdo_bc_define(eqp->default_bc,
                                eqp->n_bc_desc,
                                eqp->bc_desc,
                                n_b_faces);

  /* Values at each face (interior and border) i.e. take into account BCs */
  BFT_MALLOC(b->face_values, n_faces, cs_real_t);
# pragma omp parallel for if (n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_faces; i++) b->face_values[i] = 0;

  /* Store the last computed values of the field at cell centers and the data
     needed to compute the cell values from the vertex values.
     No need to synchronize all these quantities since they are only cellwise
     quantities. */
  BFT_MALLOC(b->cell_rhs, n_cells, cs_real_t);
  BFT_MALLOC(b->acc_inv, n_cells, cs_real_t);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++)
    b->cell_rhs[i] = b->acc_inv[i] = 0;

  BFT_MALLOC(b->acf, connect->c2f->idx[n_cells], cs_real_t);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < connect->c2f->idx[n_cells]; i++) b->acf[i] = 0.;

  /* Diffusion part */
  /* -------------- */

  b->diff_pty_uniform = true;
  b->get_stiffness_matrix = NULL;
  b->boundary_flux_op = NULL;
  b->enforce_dirichlet = NULL;

  if (eqp->flag & CS_EQUATION_DIFFUSION) {

    b->sys_flag |= CS_FLAG_SYS_DIFFUSION;
    b->diff_pty_uniform = cs_property_is_uniform(eqp->diffusion_property);

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      b->get_stiffness_matrix = cs_hodge_fb_cost_get_stiffness;
      b->boundary_flux_op = NULL; //cs_cdovb_diffusion_cost_flux_op;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      b->get_stiffness_matrix = cs_hodge_fb_voro_get_stiffness;
      b->boundary_flux_op = NULL; //cs_cdovb_diffusion_cost_flux_op;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                (" Invalid type of algorithm to build the diffusion term."));

    } // Switch on Hodge algo.

    switch (eqp->enforcement) {

    case CS_PARAM_BC_ENFORCE_WEAK_PENA:
      b->enforce_dirichlet = cs_cdo_diffusion_pena_dirichlet;
      break;

    case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    default:
      bft_error(__FILE__, __LINE__, 0,
                (" Invalid type of algorithm to enforce Dirichlet BC."));

    }

  } // Diffusion part

  /* Advection part */

  b->get_advection_matrix = NULL;
  b->add_advection_bc = NULL;
  // TODO

  /* Reaction part */
  if (eqp->n_reaction_terms > CS_CDO_N_MAX_REACTIONS)
    bft_error(__FILE__, __LINE__, 0,
              " Number of reaction terms for an equation is too high.\n"
              " Modify your settings aor contact the developpement team.");

  for (int i = 0; i < eqp->n_reaction_terms; i++)
    b->reac_pty_uniform[i]
      = cs_property_is_uniform(eqp->reaction_properties[i]);

  /* Time part */
  if (b->sys_flag & CS_FLAG_SYS_TIME)
    b->sys_flag |= CS_FLAG_SYS_TIME_DIAG;

  b->time_pty_uniform = cs_property_is_uniform(eqp->time_property);
  b->apply_time_scheme = cs_cdo_time_get_scheme_function(b->sys_flag, eqp);

  /* Source term part */

  /* Default intialization */
  b->st_msh_flag = cs_source_term_init(CS_SPACE_SCHEME_CDOFB,
                                       eqp->n_source_terms,
                   (const cs_xdef_t **)eqp->source_terms,
                                       b->compute_source,
                                       &(b->sys_flag),
                                       &(b->source_mask));

  b->source_terms = NULL;
  if (b->sys_flag & CS_FLAG_SYS_SOURCETERM) {

    BFT_MALLOC(b->source_terms, n_cells, cs_real_t);
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells; i++) b->source_terms[i] = 0;

  } /* There is at least one source term */

  /* Monitoring */
  b->monitor = cs_equation_init_monitoring();

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_scaleq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdofb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_scaleq_free(void   *builder)
{
  cs_cdofb_scaleq_t   *b  = (cs_cdofb_scaleq_t *)builder;

  if (b == NULL)
    return b;

  /* eqp, connect, quant, and time_step are only shared.
     These quantities are freed later. */

  /* Free BC structure */
  b->face_bc = cs_cdo_bc_free(b->face_bc);

  /* Monitoring structure */
  BFT_FREE(b->monitor);

  /* Free temporary buffers */
  BFT_FREE(b->source_terms);
  BFT_FREE(b->face_values);
  BFT_FREE(b->cell_rhs);
  BFT_FREE(b->acc_inv);
  BFT_FREE(b->acf);

  BFT_FREE(b);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Display information related to the monitoring of the current system
 *
 * \param[in]  eqname    name of the related equation
 * \param[in]  builder   pointer to a cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_monitor(const char   *eqname,
                        const void   *builder)
{
  const cs_cdofb_scaleq_t  *b = (const cs_cdofb_scaleq_t *)builder;

  if (b == NULL)
    return;

  cs_equation_write_monitoring(eqname, b->monitor);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in, out] builder     pointer to a cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_compute_source(void            *builder)
{
  if (builder == NULL)
    return;

  cs_cdofb_scaleq_t  *b = (cs_cdofb_scaleq_t *)builder;

  const cs_equation_param_t  *eqp = b->eqp;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  if (eqp->n_source_terms == 0)
    return;
  cs_timer_t  t0 = cs_timer_time();

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < quant->n_cells; i++) b->source_terms[i] = 0.0;

  /* Compute the source term in each cell */
  double  *contrib = cs_equation_get_tmpbuf();
  cs_flag_t  dof_loc = CS_FLAG_SCALAR | cs_cdo_primal_cell;

  for (int  st_id = 0; st_id < eqp->n_source_terms; st_id++) {

    const cs_xdef_t  *st = eqp->source_terms[st_id];

    /* contrib is updated inside this function */
    cs_source_term_compute_from_density(dof_loc, st, &contrib);

    /* Update source term array */
    for (cs_lnum_t i = 0; i < quant->n_cells; i++)
      b->source_terms[i] += contrib[i];

  } // Loop on source terms

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 2
  cs_dump_array_to_listing("INIT_SOURCE_TERM",
                           quant->n_cells, b->source_terms, 8);
#endif

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(b->monitor->tcs), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         builder structure
 *
 * \param[in, out] builder        pointer to generic builder structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_initialize_system(void           *builder,
                                  cs_matrix_t   **system_matrix,
                                  cs_real_t     **system_rhs)
{
  if (builder == NULL)
    return;
  assert(*system_matrix == NULL && *system_rhs == NULL);

  cs_cdofb_scaleq_t  *b = (cs_cdofb_scaleq_t *)builder;
  cs_timer_t  t0 = cs_timer_time();

  /* Create the matrix related to the current algebraic system */
  const cs_matrix_structure_t  *ms =
    cs_equation_get_matrix_structure(CS_SPACE_SCHEME_CDOFB);

  *system_matrix = cs_matrix_create(ms);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Allocate and initialize the related right-hand side */
  BFT_MALLOC(*system_rhs, quant->n_faces, cs_real_t);
# pragma omp parallel for if  (quant->n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < quant->n_faces; i++) (*system_rhs)[i] = 0.0;

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(b->monitor->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO face-based scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the vertex field
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to cs_cdofb_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_build_system(const cs_mesh_t       *mesh,
                             const cs_real_t       *field_val,
                             double                 dt_cur,
                             void                  *builder,
                             cs_real_t             *rhs,
                             cs_matrix_t           *matrix)
{
  /* Sanity checks */
  assert(rhs != NULL && matrix != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_timer_t  t0 = cs_timer_time();

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  cs_cdofb_scaleq_t  *b = (cs_cdofb_scaleq_t *)builder;

  /* Dirichlet values at boundary faces are first computed */
  cs_real_t  *dir_values =
    cs_equation_compute_dirichlet_sf(mesh,
                                     b->eqp,
                                     b->face_bc->dir,
                                     cs_cdofb_cell_bld[0]);

  /* Tag faces with a non-homogeneous Neumann BC */
  short int  *neu_tags = cs_equation_tag_neumann_face(b->eqp);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)         \
  shared(dt_cur, quant, connect, b, rhs, matrix, mav, dir_values, neu_tags,  \
         field_val, cs_cdofb_cell_sys, cs_cdofb_cell_bld, cs_cdofb_cell_bc)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif
    const cs_equation_param_t  *eqp = b->eqp;

    /* Test to remove */
    if (eqp->flag & CS_EQUATION_CONVECTION)
      bft_error(__FILE__, __LINE__, 0,
                _(" Convection term is not handled yet.\n"));
    if (eqp->flag & CS_EQUATION_UNSTEADY)
      bft_error(__FILE__, __LINE__, 0,
                _(" Unsteady terms are not handled yet.\n"));

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdofb_cell_sys[t_id];
    cs_cell_bc_t  *cbc = cs_cdofb_cell_bc[t_id];
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];

    /* Set inside the OMP section so that each thread has its own value */

    /* Initialization of the values of properties */
    double  time_pty_val = 1.0;
    double  reac_pty_vals[CS_CDO_N_MAX_REACTIONS];
    for (int i = 0; i < CS_CDO_N_MAX_REACTIONS; i++) reac_pty_vals[i] = 1.0;

    /* Initialize members of the builder related to the current system
       Preparatory step for diffusion term */
    if (b->sys_flag & CS_FLAG_SYS_DIFFUSION) {
      if (b->diff_pty_uniform) {

        cs_property_get_cell_tensor(0, // cell_id
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    cb->pty_mat);

        if (eqp->diffusion_hodge.is_iso)
          cb->pty_val = cb->pty_mat[0][0];

      } /* Diffusion property is uniform */
    } /* Diffusion */

    /* Preparatory step for unsteady term */
    if (b->sys_flag & CS_FLAG_SYS_TIME)
      if (b->time_pty_uniform)
        time_pty_val = cs_property_get_cell_value(0, eqp->time_property);

    /* Preparatory step for reaction term */
    if (b->sys_flag & CS_FLAG_SYS_REACTION) {

      for (int r = 0; r < eqp->n_reaction_terms; r++) {
        if (b->reac_pty_uniform[r]) {
          cs_property_t  *r_pty = eqp->reaction_properties[r];
          reac_pty_vals[r] = cs_property_get_cell_value(0, r_pty);
        }
      } // Loop on reaction properties

    } // Reaction properties

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];
      const cs_flag_t  msh_flag = _get_cell_mesh_flag(cell_flag,
                                                      b->msh_flag,
                                                      b->bd_msh_flag,
                                                      b->st_msh_flag);

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _init_cell_structures(cell_flag, cm, b,
                            dir_values, neu_tags, field_val,  // in
                            csys, cbc, cb);                   // out

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 2
      if (c_id % CS_CDOFB_SCALEQ_MODULO == 0) cs_cell_mesh_dump(cm);
#endif

      /* DIFFUSION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
      /* ============================================== */

      if (b->sys_flag & CS_FLAG_SYS_DIFFUSION) {

        /* Define the local stiffness matrix */
        if (!(b->diff_pty_uniform)) {

          cs_property_tensor_in_cell(cm,
                                     eqp->diffusion_property,
                                     eqp->diffusion_hodge.inv_pty,
                                     cb->pty_mat);

          if (eqp->diffusion_hodge.is_iso)
            cb->pty_val = cb->pty_mat[0][0];

        }

        // local matrix owned by the cellwise builder (store in cb->loc)
        b->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

        // Add the local diffusion operator to the local system
        cs_locmat_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
        if (c_id % CS_CDOFB_SCALEQ_MODULO == 0)
          cs_cell_sys_dump("\n>> Local system after diffusion", c_id, csys);
#endif
      } /* END OF DIFFUSION */

      /* SOURCE TERM COMPUTATION */
      /* ======================= */

      if (b->sys_flag & CS_FLAG_SYS_SOURCETERM) {

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */
        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (const cs_xdef_t **)eqp->source_terms,
                                        cm,
                                        b->sys_flag,
                                        b->source_mask,
                                        b->compute_source,
                                        cb,    // mass matrix is cb->hdg
                                        csys); // Fill csys->source

        if ((b->sys_flag & CS_FLAG_SYS_TIME) == 0) {
          /* Same strategy as if one applies a implicit scheme */
          csys->rhs[cm->n_fc] += csys->source[cm->n_fc];
        }

        /* Reset the value of the source term for the cell DoF
           Source term is only hold by the cell DoF in face-based schemes */
        b->source_terms[c_id] = csys->source[cm->n_fc];

      } /* End of term source contribution */


#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
      if (c_id % CS_CDOFB_SCALEQ_MODULO == 0)
        cs_cell_sys_dump(">> Local system matrix before condensation",
                         c_id, csys);
#endif

      /* Neumann boundary conditions */
      if ((cell_flag & CS_FLAG_BOUNDARY) && cbc->n_nhmg_neuman > 0) {
        for (short int f  = 0; f < cm->n_fc; f++)
          csys->rhs[cm->n_fc] += cbc->neu_values[f];
      }

      /* Static condensation of the local system matrix of size n_vc + 1 into
         a matrix of size n_vc. Store information in the builder structure in
         order to be able to compute the values at cell centers. */
      _condense_and_store(connect->c2f, b, cb, csys);

      /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
      /* ======================================================= */

      if (eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_PENA) {

        /* Weakly enforced Dirichlet BCs for cells attached to the boundary
           csys is updated inside (matrix and rhs) */
        if (cell_flag & CS_FLAG_BOUNDARY)
          b->enforce_dirichlet(eqp->diffusion_hodge, cbc, cm, // in
                               b->boundary_flux_op,           // function (in)
                               fm, cb, csys);                 // in/out

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 0
      if (c_id % CS_CDOFB_SCALEQ_MODULO == 0)
        cs_cell_sys_dump(">> (FINAL) Local system matrix", c_id, csys);
#endif

      /* Assemble the local system (related to vertices only since one applies
         a static condensation) to the global system */
      cs_equation_assemble_f(csys, connect->f_rs, b->sys_flag, // in
                             rhs, mav);       // out

    } // Main loop on cells

  } // OPENMP Block

  cs_matrix_assembler_values_done(mav); // optional

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 2
  cs_dump_array_to_listing("FINAL RHS_FACE", quant->n_faces, rhs, 8);
  if (b->source_terms != NULL)
    cs_dump_array_to_listing("FINAL RHS_CELL",
                             quant->n_cells, b->source_terms, 8);
#endif

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(neu_tags);
  cs_matrix_assembler_values_finalize(&mav);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(b->monitor->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in, out] builder    pointer to cs_cdofb_scaleq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_update_field(const cs_real_t            *solu,
                             const cs_real_t            *rhs,
                             void                       *builder,
                             cs_real_t                  *field_val)
{
  CS_UNUSED(rhs);
  cs_cdofb_scaleq_t  *b = (cs_cdofb_scaleq_t *)builder;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_adjacency_t  *c2f = cs_shared_connect->c2f;

  /* Set computed solution in builder->face_values */
  memcpy(b->face_values, solu, quant->n_faces*sizeof(cs_real_t));

  /* Build the remaining discrete operators */
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    cs_real_t  f_contrib = 0.;
    for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++)
      f_contrib += b->face_values[c2f->ids[i]] * b->acf[i];

    /* acf stores minus the row sum of the discrete Hodge operator */
    if (b->sys_flag & CS_FLAG_SYS_SOURCETERM)
      field_val[c_id] = -b->acc_inv[c_id]*(b->source_terms[c_id] + f_contrib);
    else
      field_val[c_id] = -b->acc_inv[c_id]*f_contrib;

  } // loop on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in, out]  builder    pointer to builder structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_extra_op(const char            *eqname,
                         const cs_field_t      *field,
                         void                  *builder)
{
  CS_UNUSED(eqname); // avoid a compilation warning

  char *postlabel = NULL;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_i_faces = connect->n_faces[2];
  const cs_real_t  *face_pdi = cs_cdofb_scaleq_get_face_values(builder);

  /* Field post-processing */
  int  len = strlen(field->name) + 8 + 1;
  BFT_MALLOC(postlabel, len, char);
  sprintf(postlabel, "%s.Border", field->name);

  cs_post_write_var(CS_POST_MESH_BOUNDARY,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    postlabel,
                    field->dim,
                    true,
                    true,                  // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    NULL,                  // values on cells
                    NULL,                  // values at internal faces
                    face_pdi + n_i_faces,  // values at border faces
                    cs_shared_time_step);  // time step management structure


  BFT_FREE(postlabel);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at each face
 *
 * \param[in]  builder    pointer to a cs_cdofb_scaleq_t structure
 *
 * \return  a pointer to an array of double (size n_faces)
 */
/*----------------------------------------------------------------------------*/

double *
cs_cdofb_scaleq_get_face_values(const void    *builder)
{
  const cs_cdofb_scaleq_t  *b = (const cs_cdofb_scaleq_t *)builder;

  if (b == NULL)
    return NULL;
  else
    return b->face_values;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
