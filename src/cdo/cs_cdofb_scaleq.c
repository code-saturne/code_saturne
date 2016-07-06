/*============================================================================
 * Build an algebraic CDO face-based system for convection/diffusion equation
 * with source terms
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_cdo_bc.h"
#include "cs_equation_common.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_post.h"
#include "cs_quadrature.h"
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

#define CDOFB_SCALEQ_DBG 0

/* Algebraic system for CDO face-based discretization */

struct  _cs_cdofb_scaleq_t {

  /* Pointer to a cs_equation_param_t structure shared with a cs_equation_t
     structure.  */

  const cs_equation_param_t  *eqp;

  /* Reduced system (known boundary entities may be removed --> Dirichlet) */

  cs_lnum_t  n_cells;
  cs_lnum_t  n_faces;
  cs_lnum_t  n_dof_faces; /* Number of interior faces + Neumann faces */

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

  cs_param_bc_enforce_t  enforce; // type of enforcement of BCs
  cs_cdo_bc_t           *face_bc; // list of faces sorted by type of BCs
  cs_real_t             *dir_val; // size: face_bc->dir->n_nhmg_elts

  /* Indirection between zipped numbering (without BC) and initial numbering
     Allocated only if the boundary conditions are strongly enforced.
  */
  cs_lnum_t          *f_z2i_ids;  // Mapping n_dof_faces -> n_faces
  cs_lnum_t          *f_i2z_ids;  // Mapping n_faces     -> n_dof_faces

  /* Work buffer */
  cs_real_t  *source_terms;  /* size: n_cells (sum of the contribution in each
                                cell of all the volumic source terms) */
  cs_real_t  *face_values;   /* DoF unknowns (x) + BCs */

};

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_time_step_t  *cs_shared_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

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
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect, h_info);

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
    _h = cs_hodge_build_local(c_id, connect, quant, hb);

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
  hb = cs_hodge_builder_free(hb);

  /* Take into account Dirichlet BCs to update RHS */
  if (dir_faces->n_nhmg_elts > 0) {

    cs_flag_t  dof_flag = CS_FLAG_FACE | CS_FLAG_PRIMAL | CS_FLAG_SCAL;

    cs_cdo_bc_dirichlet_set(dof_flag,
                            time_step,
                            quant,
                            eqp->bc,
                            dir_faces,
                            builder->dir_val);

  } // Dirichlet BCs with non-homogeneous values

  /* Modify the system according to the type of boundary enforcement */
  switch (builder->enforce) {

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
  cs_lnum_t  i;

  /* Sanity checks */
  assert(eqp != NULL);
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOFB);
  assert(eqp->var_type == CS_PARAM_VAR_SCAL);

  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_b_faces = mesh->n_b_faces;

  cs_cdofb_scaleq_t  *builder = NULL;

  BFT_MALLOC(builder, 1, cs_cdofb_scaleq_t);

  /* Shared pointers */
  builder->eqp = eqp;

  /* Dimensions: By default, we set number of DoFs as if there is no
     strong enforcement of the BCs */
  builder->n_cells = n_cells;
  builder->n_faces = n_faces;
  builder->n_dof_faces = n_faces;

  /* Set members and structures related to the management of the BCs */
  const cs_param_bc_t  *bc_param = eqp->bc;

  /* Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.
  */
  builder->face_bc = cs_cdo_bc_init(bc_param, n_b_faces);

  /* Strong enforcement means that we need an indirection list between the
     compress (or zip) and initial numbering of vertices */
  builder->enforce = bc_param->enforcement;

  if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_PENA)
    bft_error(__FILE__, __LINE__, 0,
              " CDO face-based schemes and weak enforcement by a strong"
              " penalization are not compatible yet.\n"
              " Please modify your settings.");

  builder->f_z2i_ids = NULL; // zipped --> initial ids
  builder->f_i2z_ids = NULL; // initial --> zipped ids

  cs_cdo_bc_list_t  *dir_faces = builder->face_bc->dir;

  BFT_MALLOC(builder->dir_val, dir_faces->n_nhmg_elts, cs_real_t);
  for (i = 0; i < dir_faces->n_nhmg_elts; i++)
    builder->dir_val[i] = 0.0;

  if (builder->enforce == CS_PARAM_BC_ENFORCE_STRONG &&
      dir_faces->n_elts > 0) {

    cs_lnum_t  cur_id = 0;
    _Bool  *is_kept = NULL;

    builder->n_dof_faces = builder->n_faces - dir_faces->n_elts;

    BFT_MALLOC(is_kept, builder->n_faces, _Bool);
    for (i = 0; i < builder->n_faces; i++)
      is_kept[i] = true;
    for (i = 0; i < dir_faces->n_elts; i++) // i_faces then b_faces
      is_kept[n_i_faces + dir_faces->elt_ids[i]] = false;

    /* Build builder->v_z2i_ids and builder->i2i_ids */
    BFT_MALLOC(builder->f_z2i_ids, builder->n_dof_faces, cs_lnum_t);
    BFT_MALLOC(builder->f_i2z_ids, builder->n_faces, cs_lnum_t);

    for (i = 0; i < builder->n_faces; i++) {
      /* by default, we consider that it's removed */
      builder->f_i2z_ids[i] = -1;
      if (is_kept[i]) {
        builder->f_i2z_ids[i] = cur_id;
        builder->f_z2i_ids[cur_id++] = i;
      }
    }
    assert(cur_id == builder->n_dof_faces);

    BFT_FREE(is_kept);

  } /* Strong enforcement of BCs */

  /* Contribution in each cell of the source terms */
  BFT_MALLOC(builder->source_terms, builder->n_cells, cs_real_t);
  for (i = 0; i < builder->n_cells; i++)
    builder->source_terms[i] = 0;

  /* Values at each face (interior and border) i.e. take into account BCs */
  BFT_MALLOC(builder->face_values, builder->n_faces, cs_real_t);
  for (i = 0; i < builder->n_faces; i++)
    builder->face_values[i] = 0;

  return builder;
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

  /* Free temporary buffers */
  BFT_FREE(_builder->source_terms);
  BFT_FREE(_builder->face_values);

  BFT_FREE(_builder);

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
cs_cdofb_scaleq_free_sysmat(void              *builder,
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
 * \param[in, out] builder     pointer to a cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_compute_source(void            *builder)
{
  cs_cdofb_scaleq_t  *b = (cs_cdofb_scaleq_t *)builder;

  for (cs_lnum_t i = 0; i < b->n_cells; i++)
    b->source_terms[i] = 0;

  const cs_equation_param_t  *eqp = b->eqp;

  if (eqp->n_source_terms == 0)
    return;

  double  *contrib = cs_equation_get_tmpbuf();
  cs_desc_t  desc = {.location = CS_FLAG_SCAL | cs_cdo_primal_cell,
                     .state = CS_FLAG_STATE_DENSITY};

  for (int  st_id = 0; st_id < eqp->n_source_terms; st_id++) {

    const cs_source_term_t  *st = eqp->source_terms[st_id];

    cs_source_term_compute(desc, st, &contrib); // updated inside this function

    /* Update source term array */
    for (cs_lnum_t i = 0; i < b->n_cells; i++)
      b->source_terms[i] += contrib[i];

  } // Loop on source terms

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO face-based scheme.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to cs_cdofb_scaleq_t structure
 * \param[in, out] rhs        pointer to a right-hand side array pointer
 * \param[in, out] sla_mat    pointer to cs_sla_matrix_t structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_build_system(const cs_mesh_t        *mesh,
                             const cs_real_t        *field_val,
                             double                  dt_cur,
                             void                   *builder,
                             cs_real_t             **rhs,
                             cs_sla_matrix_t       **sla_mat)
{
  // To be removed (avoid compilation warnings)
  CS_UNUSED(field_val);
  CS_UNUSED(dt_cur);

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

  if (*rhs == NULL)
    BFT_MALLOC(*rhs, b->n_dof_faces, cs_real_t);

  /* Build diffusion system: stiffness matrix */
  if (eqp->flag & CS_EQUATION_DIFFUSION)
    diffusion_mat = _build_diffusion_system(mesh, *rhs, b);

  *sla_mat = diffusion_mat;

  /* Build convection system */
  // TODO

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
  if (b->enforce == CS_PARAM_BC_ENFORCE_STRONG)
    for (i = 0; i < dir_faces->n_nhmg_elts; i++) // interior then border faces
      b->face_values[quant->n_i_faces + dir_faces->elt_ids[i]]
        = b->dir_val[i];

  /* Compute now the value at each cell center */
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect, h_info);

  /* Build the remaining discrete operators */
  for (c_id = 0; c_id < b->n_cells; c_id++) {

    int  shft = connect->c2f->idx[c_id];
    double _wf_val = 0.0, dsum = 0.0, rowsum = 0.0;

    /* Build a local discrete Hodge operator */
    cs_locmat_t  *_h = cs_hodge_build_local(c_id, connect, quant, hb);

    /* Compute dsum: the sum of all the entries of the local discrete Hodge
       operator */
    for (i = 0, l=shft; i < _h->n_ent; i++, l++) {
      rowsum = 0;
      f_id = connect->c2f->col_id[l];
      for (j = 0; j < _h->n_ent; j++)
        rowsum += _h->val[i*_h->n_ent+j];
      dsum += rowsum;
      _wf_val += b->face_values[f_id] * rowsum;
    }

    field_val[c_id] = 1/dsum*(b->source_terms[c_id] + _wf_val);

  } // loop on cells

  /* Free memory */
  hb = cs_hodge_builder_free(hb);
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
cs_cdofb_scaleq_extra_op(const char            *eqname,
                         const cs_field_t      *field,
                         void                  *builder)
{
  CS_UNUSED(eqname); // avoid a compilation warning

  char *postlabel = NULL;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_i_faces = connect->f_info->n_i_elts;
  const cs_real_t  *face_pdi = cs_cdofb_scaleq_get_face_values(builder);

  /* Field post-processing */
  int  len = strlen(field->name) + 8 + 1;
  BFT_MALLOC(postlabel, len, char);
  sprintf(postlabel, "%s.Border", field->name);

  cs_post_write_var(-2,                    // id du maillage de post
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
