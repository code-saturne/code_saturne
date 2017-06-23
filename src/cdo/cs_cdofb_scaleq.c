/*============================================================================
 * Build an algebraic CDO face-based system for convection/diffusion equation
 * with source terms
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

#define CS_CDOFB_SCALEQ_DBG 0

/* Algebraic system for CDO face-based discretization */

struct  _cs_cdofb_scaleq_t {

  /* Pointer to a cs_equation_param_t structure shared with a cs_equation_t
     structure.  */

  const cs_equation_param_t  *eqp;

  /* Dimensions of the algebraic system */
  cs_lnum_t            n_cells;
  cs_lnum_t            n_faces;
  short int            max_sys_size;

  /* Shortcut to know what to build */
  cs_flag_t            cm_flag;   // Information related to cell mesh
  cs_flag_t            sys_flag;  // Information related to the sytem

  /* Common members for all terms */

  /* Builder structure for advection term */

  /* Builder structure for diffusion term */
  bool                 diff_pty_uniform;

  /* Time term */
  bool           time_pty_uniform;
  double         time_pty_val;

  /* Reaction terms */
  bool          *reaction_pty_uniform;
  double        *reaction_pty_val;

  /* Source terms */
  cs_real_t     *source_terms; /* Array storing the value arising from the
                                  contribution of all source terms */
  cs_mask_t     *source_mask;  /* NULL if at least one source term is not
                                  defined for all cells (size = n_cells) */
  /* Metadata related to where and how is defined a source term */
  cs_flag_t      st_flags[CS_N_MAX_SOURCE_TERMS];
  cs_xdef_t      st_desc[CS_N_MAX_SOURCE_TERMS];

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

  /* Solution of the algebraic system at the last iteration */
  cs_real_t             *face_values;  /* DoF unknowns (x) + BCs */

  /* Temporary buffers */
  double                 *loc_vals; // local values

  /* Temporary ==> TO REMOVE */
  cs_real_t        *dir_val; // TO BE REMOVED
  cs_lnum_t         n_dof_faces;  // TO BE REMOVED
  cs_lnum_t        *f_z2i_ids;
  cs_lnum_t        *f_i2z_ids;

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

  /* If one needs to build a local hodge op. for time and reaction */
  cs_param_hodge_t                 hdg_wbs;
  cs_hodge_t                      *get_mass_matrix;

  /* Monitoring the efficiency */
  cs_equation_monitor_t           *monitor;

};

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_time_step_t  *cs_shared_time_step;

static cs_cell_sys_t  **cs_cdofb_cell_sys = NULL;
static cs_cell_bc_t  **cs_cdofb_cell_bc = NULL;

/* Flag to indicate which members have to be built in a cs_cell_mesh_t
   structure */
static const cs_flag_t  cs_cdofb_cmflag =
  CS_CDO_LOCAL_PF | CS_CDO_LOCAL_PE | CS_CDO_LOCAL_FE;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Manage potential threading and return a cs_cell_mesh_t structure
 *
 * \return a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_cell_mesh_t *
_get_cell_mesh(void)
{
  cs_cell_mesh_t  *cm = NULL;
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    cm = cs_cdo_local_get_cell_mesh(omp_get_thread_num());
  }
#else
  return  cs_cdo_local_get_cell_mesh(0);
#endif /* openMP ? */
  return  cm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Manage potential threading and return a cs_face_mesh_t structure
 *
 * \return a cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_face_mesh_t *
_get_face_mesh(void)
{
  cs_face_mesh_t  *fm = NULL;
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int  t_id = omp_get_thread_num();
    fm = cs_cdo_local_get_face_mesh(t_id);
  }
#else
  return  cs_cdo_local_get_face_mesh(0);
#endif /* openMP ? */
  return  fm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Manage potential threading and return a cs_cell_sys_t structure
 *
 * \return a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_cell_sys_t *
_get_cell_system(void)
{
  cs_cell_sys_t *csys = NULL;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int  t_id = omp_get_thread_num();
    csys = cs_cdofb_cell_sys[t_id];
  }
#else
  csys = cs_cdofb_cell_sys[0];
#endif /* openMP ? */

  return csys;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Manage potential threading and return a cs_cell_bc_t structure
 *
 * \return a cs_cell_bc_t structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_cell_bc_t *
_get_cell_bc(void)
{
  cs_cell_bc_t  *cbc = NULL;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int  t_id = omp_get_thread_num();
    cbc = cs_cdofb_cell_bc[t_id];
  }
#else
  cbc = cs_cdofb_cell_bc[0];
#endif /* openMP ? */

  return cbc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize the matrix related to the diffusion op.
 *          Note: values are filled in a second step
 *
 * \param[in]    connect   pointer to a cs_cdo_connect_t structure
 * \param[in]    quant     pointer to a cs_cdo_quantities_t structure
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_diffusion_matrix(const cs_cdo_connect_t     *connect,
                       const cs_cdo_quantities_t  *quant)
{
  int  i, j, shift;

  cs_connect_index_t  *f2f = NULL, *c2f = NULL, *f2c = NULL;

  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_sla_matrix_t *mc2f = connect->c2f;
  const cs_sla_matrix_t *mf2c = connect->f2c;

  /* Allocate and initialize the matrix */
  cs_sla_matrix_t  *mat = cs_sla_matrix_create(n_faces, n_faces, 1,
                                               CS_SLA_MAT_MSR,
                                               false);

  /* Build a face -> face connectivity */
  f2c = cs_index_map(mf2c->n_rows, mf2c->idx, mf2c->col_id);
  c2f = cs_index_map(mc2f->n_rows, mc2f->idx, mc2f->col_id);
  f2f = cs_index_compose(n_faces, f2c, c2f);
  cs_index_sort(f2f);
  mat->flag |= CS_SLA_MATRIX_SORTED;

  /* Update index: f2f has the diagonal entry. Remove it for the Hodge index */
  mat->idx[0] = 0;
  for (i = 0; i < n_faces; i++)
    mat->idx[i+1] = mat->idx[i] + f2f->idx[i+1]-f2f->idx[i]-1;

  /* Fill column ids */
  BFT_MALLOC(mat->col_id, mat->idx[n_faces], cs_lnum_t);
  shift = 0;
  for (i = 0; i < n_faces; i++)
    for (j = f2f->idx[i]; j < f2f->idx[i+1]; j++)
      if (f2f->ids[j] != i)
        mat->col_id[shift++] = f2f->ids[j];

  /* Sanity check */
  assert(shift == mat->idx[n_faces]);

  /* Free temporary memory */
  cs_index_free(&f2f);
  cs_index_free(&f2c);
  cs_index_free(&c2f);

  /* Allocate and initialize value array */
  for (i = 0; i < n_faces; i++)
    mat->diag[i] = 0.0;

  BFT_MALLOC(mat->val, mat->idx[n_faces], double);
  for (i = 0; i < mat->idx[n_faces]; i++)
    mat->val[i] = 0.0;

  return mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the final (reduced) matrix for diffusion and its right hand
 *          side (RHS). RHS is the sum of three contributions
 *           - Source terms
 *           - Neumann boundary conditions
 *           - Dirichlet boundary conditions
 *
 * \param[in]      m           pointer to a cs_mesh_t structure
 * \param[in, out] rhs         right-hand side
 * \param[in, out] builder     pointer to a cs_cdofb_scaleq_t struct.
 *
 * \return a pointer to the full stiffness matrix
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_build_diffusion_system(const cs_mesh_t             *m,
                        cs_real_t                   *rhs,
                        cs_cdofb_scaleq_t           *builder)
{
  int  i, j, ij, c_id;
  double  dsum, rowsum, invdsum;

  double  *BHCtc = NULL; // local size arrays
  cs_sla_matrix_t  *final_matrix = NULL;
  cs_locmat_t  *_h = NULL;

  const cs_cdo_bc_list_t  *dir_faces = builder->face_bc->dir;
  const cs_equation_param_t  *eqp = builder->eqp;
  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const cs_time_step_t  *time_step = cs_shared_time_step;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_cells = quant->n_cells;

  cs_sla_matrix_t  *full_matrix = _init_diffusion_matrix(connect, quant);
  cs_locmat_t  *_a = cs_locmat_create(connect->n_max_fbyc);

  /* Define a builder for the related discrete Hodge operator */
  //  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect, false, h_info);

  /* Sanity check */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EDFP);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);

  /* Buffers stored in _fbscal_work */
  double  *work = cs_equation_get_tmpbuf();
  double  *contrib = work;                      // size: n_faces
  double  *face_rhs = work + builder->n_faces;  // size: n_faces

  for (i = 0; i < 2*builder->n_faces; i++)
    contrib[i] = 0;

  /*  Build full-size operators:

           n_cells      n_i_faces
         <--------->   <--------->
        | \            * * * * * * |
        |   \           *  BtHC *  |
        |   DcHcUc     * * * * * * |
        |       \       * * * * *  |
    A = |         \                |
        |* * * * * *               |
        | * CHBt * *         H     |
        |* * * * * *               |
        | * * * * *                |

  */

  /* Allocate local operators */
  BFT_MALLOC(BHCtc, connect->n_max_fbyc, double);

  /* Build the remaining discrete operators */
  for (c_id = 0; c_id < n_cells; c_id++) {

    /* Build a local discrete Hodge operator and return a local dense matrix */
    //    _h = cs_hodge_build_local(c_id, connect, quant, hb);

    /* Compute dsum = Dc*_H*Uc where Uc = transpose(Dc) */
    dsum = 0;
    _a->n_ent = _h->n_ent;

    for (i = 0; i < _h->n_ent; i++) {
      rowsum = 0;
      _a->ids[i] = _h->ids[i];

      for (j = 0; j < _h->n_ent; j++)
        rowsum += _h->val[i*_h->n_ent+j];

      dsum += rowsum;
      BHCtc[i] = -rowsum;

    }
    invdsum = 1/dsum;

    /* Define local diffusion matrix */
    for (i = 0; i < _a->n_ent; i++) {
      for (j = 0; j < _a->n_ent; j++) {
        ij = i*_a->n_ent+j;
        _a->val[ij] = -BHCtc[j]*invdsum*BHCtc[i];
        _a->val[ij] += _h->val[ij];
      }
    }

    /* Assemble local stiffness matrix */
    cs_sla_assemble_msr_sym(_a, full_matrix, false); // Not only diag. terms

    /* Assemble RHS (source term contribution) */
    for (i = 0; i < _a->n_ent; i++)
      face_rhs[_a->ids[i]] -= BHCtc[i]*invdsum*builder->source_terms[c_id];

  } /* End of loop on cells */

  /* Clean entries of the operators */
  // cs_sla_matrix_clean(full_matrix, cs_math_get_machine_epsilon());

  /* Free memory */
  BFT_FREE(BHCtc);
  _a = cs_locmat_free(_a);
  //  hb = cs_hodge_builder_free(hb);

  /* Take into account Dirichlet BCs to update RHS */
  if (dir_faces->n_nhmg_elts > 0) {

    cs_flag_t  dof_flag = CS_FLAG_FACE | CS_FLAG_PRIMAL | CS_FLAG_SCALAR;

    /* cs_cdo_bc_dirichlet_set(dof_flag, */
    /*                         time_step, */
    /*                         quant, */
    /*                         eqp->bc, */
    /*                         dir_faces, */
    /*                         builder->dir_val); */

  } // Dirichlet BCs with non-homogeneous values

  /* Modify the system according to the type of boundary enforcement */
  switch (eqp->enforcement) {

  case CS_PARAM_BC_ENFORCE_STRONG:
    {
      if (dir_faces->n_nhmg_elts > 0) {

        double  *x_bc = work + 2*builder->n_faces;  // size: n_faces

        for (i = 0; i < builder->n_faces; i++)
          x_bc[i] = 0;
        for (i = 0; i < dir_faces->n_nhmg_elts; i++) // interior then border
          x_bc[m->n_i_faces + dir_faces->elt_ids[i]] = builder->dir_val[i];

        cs_sla_matvec(full_matrix, x_bc, &contrib, true);
        for (i = 0; i < builder->n_faces; i++)
          face_rhs[i] -= contrib[i];

      } // Dirichlet BCs

      if (builder->n_dof_faces < builder->n_faces) { /* Reduce system size */

        for (i = 0; i < builder->n_dof_faces; i++)
          rhs[i] = face_rhs[builder->f_z2i_ids[i]];

        final_matrix = cs_sla_matrix_pack(builder->n_dof_faces,  // n_block_rows
                                          builder->n_dof_faces,  // n_block_cols
                                          full_matrix,           // full matrix
                                          builder->f_z2i_ids,    // row_z2i_ids
                                          builder->f_i2z_ids,    // col_i2z_ids
                                          true);                 // keep sym.

        /* Free buffers */
        full_matrix = cs_sla_matrix_free(full_matrix);

      }

    }
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_PENA:
    {
      assert(builder->n_faces == builder->n_dof_faces); /* Sanity check */

      cs_real_t  pena_coef = 0.01/cs_math_get_machine_epsilon();

      for (i = 0; i < dir_faces->n_nhmg_elts; i++)
        face_rhs[dir_faces->elt_ids[i]] += pena_coef * builder->dir_val[i];

      for (i = 0; i < dir_faces->n_nhmg_elts; i++)
        full_matrix->diag[dir_faces->elt_ids[i]] += pena_coef;

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " This kind of BC enforcement is not implemented yet.\n"
              " Please modify your settings.");

  } // End of switch on bc enforcement

  if (builder->n_faces == builder->n_dof_faces) {
    final_matrix = full_matrix;
    memcpy(rhs, face_rhs, sizeof(cs_real_t)*builder->n_faces);
  }

  return final_matrix;
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
  return; // Nothing to do up to now
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
  return;
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
  /* Sanity checks */
  assert(eqp != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOFB && eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Invalid type of equation.\n"
              " Expected: scalar-valued CDO face-based equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_faces = connect->n_faces[0];
  const cs_lnum_t  n_b_faces = connect->n_faces[1];
  const cs_lnum_t  n_i_faces = connect->n_faces[2];

  cs_cdofb_scaleq_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_cdofb_scaleq_t);

  /* Shared pointers */
  b->eqp = eqp;

  /* Dimensions of the algebraic system */
  b->n_cells = n_cells;
  b->n_faces = n_faces;
  b->max_sys_size = connect->n_max_fbyc + 1;

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


  /* Set members and structures related to the management of the BCs
     Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.  */
  b->face_bc = cs_cdo_bc_define(eqp->default_bc,
                                eqp->n_bc_desc,
                                eqp->bc_desc,
                                n_b_faces);

  if (eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_PENA)
    bft_error(__FILE__, __LINE__, 0,
              " CDO face-based schemes and weak enforcement by a strong"
              " penalization are not compatible yet.\n"
              " Please modify your settings.");

  /* Initialization of members common to several terms */
  b->sys_flag = 0;
  b->cm_flag = cs_cdofb_cmflag;

  BFT_MALLOC(b->loc_vals, 3*b->max_sys_size, double);
  for (int i = 0; i < 3*b->max_sys_size; i++)
    b->loc_vals[i] = 0;

  /* Diffusion part */
  /* -------------- */

  b->diff_pty_uniform = false;
  if (b->sys_flag & CS_FLAG_SYS_DIFFUSION) {

    bool is_uniform = cs_property_is_uniform(eqp->diffusion_property);
    b->diff_pty_uniform = is_uniform;

    bool is_isotropic = false;
    if (cs_property_get_type(eqp->diffusion_property) == CS_PROPERTY_ISO)
      is_isotropic = true;


  }

  /* Advection part */
  /* -------------- */

  /* b->adv = NULL; */
  /* if (b->sys_flag & CS_FLAG_SYS_ADVECTION) */
  /*   b->adv = cs_cdo_advection_builder_init(connect, eqp, b->has[CS_FLAG_SYS_DIFFUSION]); */
  /* else { */
  /*   if (eqp->enforcement != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE) */
  /*     b->flag |= CS_FLAG_SYS_SYM; // Algebraic system is symmetric */
  /* } */

  /* Time part */
  /* --------- */

  b->time_pty_uniform = false;
  b->time_pty_val = 0.;
  if (b->sys_flag & CS_FLAG_SYS_TIME) {

    b->time_pty_uniform = cs_property_is_uniform(eqp->time_property);
    if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_VORONOI)
      b->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
  }

  /* Source term part */
  /* ---------------- */

/*   /\* Default intialization *\/ */
/*   b->st_msh_flag = cs_source_term_init(CS_SPACE_SCHEME_CDOVCB, */
/*                       eqp->n_source_terms, */
/*                       eqp->source_terms, */
/*                       b->compute_source, */
/*                       &(b->sys_flag), */
/*                       &(b->source_mask)); */

/*   b->source_terms = NULL; */
/*   if (b->sys_flag & CS_FLAG_SYS_SOURCETERM) { */

/*     BFT_MALLOC(b->source_terms, b->n_dofs, cs_real_t); */
/* # pragma omp parallel for if (b->n_dofs > CS_THR_MIN) */
/*     for (cs_lnum_t i = 0; i < b->n_dofs; i++) */
/*       b->source_terms[i] = 0; */

/*   } /\* There is at least one source term *\/ */

  /* Values at each face (interior and border) i.e. take into account BCs */
  BFT_MALLOC(b->face_values, b->n_faces, cs_real_t);
# pragma omp parallel for if (b->n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->n_faces; i++)
    b->face_values[i] = 0;

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
  cs_cdofb_scaleq_t   *_builder  = (cs_cdofb_scaleq_t *)builder;

  if (_builder == NULL)
    return _builder;

  /* eqp, connect, quant, and time_step are only shared.
     These quantities are freed later. */

  /* Free BC structure */
  if (_builder->face_bc->dir->n_nhmg_elts > 0)
    BFT_FREE(_builder->dir_val);
  _builder->face_bc = cs_cdo_bc_free(_builder->face_bc);

  /* Renumbering */
  if (_builder->n_faces > _builder->n_dof_faces) {
    BFT_FREE(_builder->f_z2i_ids);
    BFT_FREE(_builder->f_i2z_ids);
  }

  /* Monitoring structure */
  BFT_FREE(_builder->monitor);

  /* Free temporary buffers */
  BFT_FREE(_builder->source_terms);
  BFT_FREE(_builder->face_values);

  BFT_FREE(_builder);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Display information related to the monitoring of the current system
 *
 * \param[in]  eqname    name of the related equation
 * \param[in]  builder   pointer to a cs_cdovcb_scaleq_t structure
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
  cs_flag_t  cm_flag = cs_cdofb_cmflag;
  cs_cell_mesh_t  *cm = _get_cell_mesh();
  cs_cell_sys_t  *csys = _get_cell_system();

  for (cs_lnum_t i = 0; i < b->n_cells; i++)
    b->source_terms[i] = 0;

  const cs_equation_param_t  *eqp = b->eqp;

  if (eqp->n_source_terms == 0)
    return;

  double  *contrib = cs_equation_get_tmpbuf();
  cs_flag_t  dof_loc = CS_FLAG_SCALAR | cs_cdo_primal_cell;

  for (int  st_id = 0; st_id < eqp->n_source_terms; st_id++) {

    const cs_xdef_t  *st = eqp->source_terms[st_id];

    /* contrib is updated inside this function */
    cs_source_term_compute_from_density(dof_loc, st, &contrib);

    /* Update source term array */
    for (cs_lnum_t i = 0; i < b->n_cells; i++)
      b->source_terms[i] += contrib[i];

  } // Loop on source terms

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

  /* Allocate and initialize the related right-hand side */
  BFT_MALLOC(*system_rhs, b->n_faces, cs_real_t);
#pragma omp parallel for if  (b->n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->n_faces; i++) (*system_rhs)[i] = 0.0;

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
 * \param[in, out] builder    pointer to cs_cdovcb_scaleq_t structure
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
  // To be removed (avoid compilation warnings)
  CS_UNUSED(field_val);
  CS_UNUSED(dt_cur);

  cs_timer_t  t0 = cs_timer_time();

  cs_sla_matrix_t  *diffusion_mat = NULL;
  cs_cdofb_scaleq_t  *b = (cs_cdofb_scaleq_t *)builder;

  const cs_equation_param_t  *eqp = b->eqp;

  /* Test to remove */
  if (eqp->flag & CS_EQUATION_CONVECTION)
    bft_error(__FILE__, __LINE__, 0,
              _(" Convection term is not handled yet.\n"));
  if (eqp->flag & CS_EQUATION_UNSTEADY)
    bft_error(__FILE__, __LINE__, 0,
              _(" Unsteady terms are not handled yet.\n"));

  /* Build diffusion system: stiffness matrix */
  if (eqp->flag & CS_EQUATION_DIFFUSION)
    diffusion_mat = _build_diffusion_system(mesh, rhs, b);

  /* Build convection system */
  // TODO

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
 * \param[in, out] builder    pointer to cs_cdovb_scaleq_t structure
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

  int  i, j, l, c_id, f_id;

  cs_cdofb_scaleq_t  *b = (cs_cdofb_scaleq_t *)builder;

  const cs_cdo_bc_list_t  *dir_faces = b->face_bc->dir;
  const cs_equation_param_t  *eqp = b->eqp;
  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Set computed solution in builder->face_values */
  if (b->n_dof_faces < b->n_faces) {
    for (i = 0; i < b->n_faces; i++)
      b->face_values[i] = 0;
    for (i = 0; i < b->n_dof_faces; i++)
      b->face_values[b->f_z2i_ids[i]] = solu[i];
  }
  else
    memcpy(b->face_values, solu, b->n_faces*sizeof(cs_real_t));

  /* Take into account Dirichlet BCs */
  if (eqp->enforcement == CS_PARAM_BC_ENFORCE_STRONG)
    for (i = 0; i < dir_faces->n_nhmg_elts; i++) // interior then border faces
      b->face_values[quant->n_i_faces + dir_faces->elt_ids[i]]
        = b->dir_val[i];

  /* /\* Compute now the value at each cell center *\/ */
  /* cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect, false, h_info); */

  /* Build the remaining discrete operators */
  for (c_id = 0; c_id < b->n_cells; c_id++) {

    int  shft = connect->c2f->idx[c_id];
    double _wf_val = 0.0, dsum = 0.0, rowsum = 0.0;

    /* Build a local discrete Hodge operator */
    /* cs_locmat_t  *_h = cs_hodge_build_local(c_id, connect, quant, hb); */

    /* Compute dsum: the sum of all the entries of the local discrete Hodge
       operator */
    /* for (i = 0, l=shft; i < _h->n_ent; i++, l++) { */
    /*   rowsum = 0; */
    /*   f_id = connect->c2f->col_id[l]; */
    /*   for (j = 0; j < _h->n_ent; j++) */
    /*     rowsum += _h->val[i*_h->n_ent+j]; */
    /*   dsum += rowsum; */
    /*   _wf_val += b->face_values[f_id] * rowsum; */
    /* } */

    field_val[c_id] = 1/dsum*(b->source_terms[c_id] + _wf_val);

  } // loop on cells

  /* /\* Free memory *\/ */
  /* hb = cs_hodge_builder_free(hb); */
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
