/*============================================================================
 * Routines to handle the evaluation of boundary conditions when building the
 * algebraic system in CDO/HHO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_equation_common.h"
#include "cs_evaluate.h"
#include "cs_xdef.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_bc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

#define CS_EQUATION_BC_DBG  0

/*============================================================================
 * Local private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the Dirichlet BC values to the face vertices.
 *          Case of vertex-based schemes
 *
 * \param[in]      dim          number of values to assign to each vertex
 * \param[in]      n_vf         number of vertices in a face
 * \param[in]      lst          list of vertex numbering
 * \param[in]      eval         result of the evaluation to set
 * \param[in]      is_constant  same value for all vertices ?
 * \param[in, out] vvals        vertex values to update
 * \param[in, out] counter      counter to update
 */
/*----------------------------------------------------------------------------*/

static inline void
_assign_vb_dirichlet_values(int                dim,
                            int                n_vf,
                            const cs_lnum_t   *lst,
                            const cs_real_t   *eval,
                            bool               is_constant,
                            cs_real_t         *vvals,
                            int                counter[])
{
  if (dim == 1) {

    for (short int v = 0; v < n_vf; v++) {
      const cs_lnum_t  v_id = lst[v];
      const short int  _v = is_constant ? 0 : v;
      counter[v_id] += 1;
      vvals[v_id] += eval[_v];
    }

  }
  else { /* Not a scalar-valued quantity */

    for (short int v = 0; v < n_vf; v++) {
      const cs_lnum_t  v_id = lst[v];
      const short int  _v = is_constant ? 0 : v;
      counter[v_id] += 1;
      for (int k = 0; k < dim; k++)
        vvals[dim*v_id + k] += eval[dim*_v + k];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set members of the cs_cell_sys_t structure related to the boundary
 *         conditions. Only the generic part is done here. The remaining part
 *         is performed specifically for each scheme
 *
 * \param[in]      face_bc   pointer to a cs_cdo_bc_face_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in, out] csys      pointer to a cs_cell_system_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_init_cell_sys_bc(const cs_cdo_bc_face_t     *face_bc,
                  const cs_cell_mesh_t       *cm,
                  cs_cell_sys_t              *csys)
{
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_lnum_t  bf_id = cm->f_ids[f] - cm->bface_shift;

    csys->bf_ids[f] = bf_id;

    if (bf_id > -1) { /* This is a boundary face */
      csys->bf_flag[f] = face_bc->flag[bf_id];
      csys->_f_ids[csys->n_bc_faces] = f;
      csys->n_bc_faces++;
    }

  } /* Loop on cell faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the boundary definitions related to the enforcement of
 *         a circulation along a boundary edge
 *
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2e_idx   index array  to define
 * \param[in, out]  def2e_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

static void
_sync_circulation_def_at_edges(const cs_cdo_connect_t    *connect,
                               int                        n_defs,
                               cs_xdef_t                **defs,
                               cs_lnum_t                  def2e_idx[],
                               cs_lnum_t                  def2e_ids[])
{
  if (n_defs == 0)
    return;

  const cs_lnum_t  n_edges = connect->n_edges;
  const cs_adjacency_t  *f2e = connect->f2e;

  int  *e2def_ids = NULL;
  BFT_MALLOC(e2def_ids, n_edges, int);
# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t e = 0; e < n_edges; e++)
    e2def_ids[e] = -1; /* default: not associated to a definition */

  const cs_lnum_t  face_shift = connect->n_faces[CS_INT_FACES];

  for (int def_id = 0; def_id < n_defs; def_id++) {

    /* Get and then set the definition of the initial condition */
    const cs_xdef_t  *def = defs[def_id];
    assert(def->support == CS_XDEF_SUPPORT_BOUNDARY);

    if (def->meta & CS_CDO_BC_TANGENTIAL_DIRICHLET) {

      const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);

      for (cs_lnum_t i = 0; i < z->n_elts; i++) { /* Loop on selected cells */
        const cs_lnum_t  f_id = face_shift + z->elt_ids[i];
        for (cs_lnum_t j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++)
          e2def_ids[f2e->ids[j]] = def_id;
      }

    } /* Enforcement of the tangential component */

  } /* Loop on definitions */

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    assert(connect->interfaces[CS_CDO_CONNECT_EDGE_SCAL] != NULL);
    /* Last definition is used if there is a conflict between several
       definitions */
    cs_interface_set_max(connect->interfaces[CS_CDO_CONNECT_EDGE_SCAL],
                         n_edges,
                         1,             /* stride */
                         false,         /* interlace (not useful here) */
                         CS_INT_TYPE,   /* int */
                         e2def_ids);

  }

  /* 0. Initialization */
  cs_lnum_t  *count = NULL;
  BFT_MALLOC(count, n_defs, cs_lnum_t);
  memset(count, 0, n_defs*sizeof(cs_lnum_t));
  memset(def2e_idx, 0, (n_defs+1)*sizeof(cs_lnum_t));

  /* 1. Count the number of edges related to each definition */
  for (cs_lnum_t e = 0; e < n_edges; e++)
    if (e2def_ids[e] > -1)
      def2e_idx[e2def_ids[e]+1] += 1;

  /* 2. Build the index */
  for (int def_id = 0; def_id < n_defs; def_id++)
    def2e_idx[def_id+1] += def2e_idx[def_id];

  /* 3. Build the list */
  for (cs_lnum_t e = 0; e < n_edges; e++) {
    const int def_id = e2def_ids[e];
    if (def_id > -1) {
      def2e_ids[def2e_idx[def_id] + count[def_id]] = e;
      count[def_id] += 1;
    }
  }

  /* Free memory */
  BFT_FREE(e2def_ids);
  BFT_FREE(count);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the values for the normal boundary flux stemming from the
 *         Neumann boundary conditions (zero is left where a Dirichlet is
 *         set. This can be updated later on)
 *
 * \param[in]       t_eval   time at which one performs the evaluation
 * \param[in]       cdoq     pointer to a cs_cdo_quantities_t structure
 * \param[in]       eqp      pointer to a cs_equation_param_t structure
 * \param[in, out]  values   pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_boundary_flux_from_bc(cs_real_t                    t_eval,
                                       const cs_cdo_quantities_t   *cdoq,
                                       const cs_equation_param_t   *eqp,
                                       cs_real_t                   *values)
{
  /* We assume a homogeneous Neumann boundary condition as a default */
  memset(values, 0, sizeof(cs_real_t)*cdoq->n_b_faces);

  for (int def_id = 0; def_id < eqp->n_bc_defs; def_id++) {

    const cs_xdef_t  *def = eqp->bc_defs[def_id];
    const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);
    assert(def->support == CS_XDEF_SUPPORT_BOUNDARY);

    if (cs_flag_test(def->meta, CS_CDO_BC_NEUMANN)) {

      switch (def->type) {

      case CS_XDEF_BY_VALUE:
        {
          const cs_real_t  *constant_val = (cs_real_t *)def->input;

          switch (eqp->dim) {

          case 1: /* scalar-valued equation */
#           pragma omp parallel for if (z->n_elts > CS_THR_MIN)
            for (cs_lnum_t i = 0; i < z->n_elts; i++) {
              const cs_lnum_t  elt_id =
                (z->elt_ids != NULL) ? z->elt_ids[i] : i;
              values[elt_id] = constant_val[0];
            }
            break;

          default:
#           pragma omp parallel for if (z->n_elts > CS_THR_MIN)
            for (cs_lnum_t i = 0; i < z->n_elts; i++) {
              const cs_lnum_t  elt_id =
                (z->elt_ids != NULL) ? z->elt_ids[i] : i;
              for (int k = 0; k < eqp->dim; k++)
                values[eqp->dim*elt_id + k] = constant_val[k];
            }
            break;

          } /* switch on dimension */

        }
        break; /* definition by value */

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        {
          cs_xdef_analytic_input_t  *anai =
            (cs_xdef_analytic_input_t *)def->input;

          anai->func(t_eval,
                     z->n_elts, z->elt_ids, cdoq->b_face_center,
                     false,       /* compacted output ? */
                     anai->input,
                     values);
        }
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.", __func__);
      }

    } /* Neumann boundary conditions */

  } /* Loop on boundary conditions applied to the Richards equation */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the BC into a cellwise view of the current system.
 *          Case of vertex-based schemes
 *
 * \param[in]      cm           pointer to a cellwise view of the mesh
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      face_bc      pointer to a cs_cdo_bc_face_t structure
 * \param[in]      vtx_bc_flag  BC flags associated to vertices
 * \param[in]      dir_values   Dirichlet values associated to each vertex
 * \param[in]      t_eval       time at which one performs the evaluation
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_vb_set_cell_bc(const cs_cell_mesh_t         *cm,
                           const cs_equation_param_t    *eqp,
                           const cs_cdo_bc_face_t       *face_bc,
                           const cs_flag_t               vtx_bc_flag[],
                           const cs_real_t               dir_values[],
                           cs_real_t                     t_eval,
                           cs_cell_sys_t                *csys,
                           cs_cell_builder_t            *cb)
{
  CS_UNUSED(cb);

  /* Sanity check */
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_EV | CS_FLAG_COMP_FE));

  /* Initialize the common part */
  _init_cell_sys_bc(face_bc, cm, csys);

  const int  d = eqp->dim;

  /* First pass: Set the DoF flag and the Dirichlet values */
  for (short int v = 0; v < cm->n_vc; v++) {

    const cs_flag_t bc_flag = vtx_bc_flag[cm->v_ids[v]];
    for (int k = 0; k < d; k++)
      csys->dof_flag[d*v+k] = bc_flag;

    if (cs_cdo_bc_is_dirichlet(bc_flag)) {
      csys->has_dirichlet = true;
      if (bc_flag & CS_CDO_BC_HMG_DIRICHLET)
        continue; /* Nothing else to do for this vertex */
      else {
        assert(bc_flag & CS_CDO_BC_DIRICHLET);
        for (int k = 0; k < d; k++)
          csys->dir_values[d*v+k] = dir_values[d*cm->v_ids[v]+k];
      }
    }

  } /* Loop on cell vertices */

  /* Second pass: BC related to faces. Neumann or Robin or sliding condition
     (for vector-valued equations). Identify which face is a boundary face */
  for (short int f = 0; f < cm->n_fc; f++) {
    if (csys->bf_ids[f] > -1) { /* This a boundary face */

      switch(csys->bf_flag[f]) {

      case CS_CDO_BC_NEUMANN:
        csys->has_nhmg_neumann = true;
        cs_equation_compute_neumann_sv(t_eval,
                                       face_bc->def_ids[csys->bf_ids[f]],
                                       f,
                                       eqp,
                                       cm,
                                       csys->neu_values);
        break;

      case CS_CDO_BC_ROBIN:
        csys->has_robin = true;

        /* The values which define the Robin BC are stored for each boundary
           face. This is different from Neumann and Dirichlet where the values
           are defined at each vertices. Be aware of that when computing */
        cs_equation_compute_robin(t_eval,
                                  face_bc->def_ids[csys->bf_ids[f]],
                                  f,
                                  eqp,
                                  cm,
                                  csys->rob_values);
        break;

      case CS_CDO_BC_SLIDING:
        csys->has_sliding = true;
        break;

      default:   /* Nothing to do for */
        /* case CS_CDO_BC_HMG_DIRICHLET: */
        /* case CS_CDO_BC_DIRICHLET: */
        /* case CS_CDO_BC_HMG_NEUMANN: */
        break;

      } /* End of switch */

    } /* Boundary face */
  } /* Loop on cell faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the BC into a cellwise view of the current system.
 *          Case of Face-based schemes
 *
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      face_bc     pointer to a cs_cdo_bc_face_t structure
 * \param[in]      dir_values  Dirichlet values associated to each vertex
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_fb_set_cell_bc(const cs_cell_mesh_t         *cm,
                           const cs_equation_param_t    *eqp,
                           const cs_cdo_bc_face_t       *face_bc,
                           const cs_real_t               dir_values[],
                           cs_real_t                     t_eval,
                           cs_cell_sys_t                *csys,
                           cs_cell_builder_t            *cb)
{
  CS_UNUSED(cb);

  /* Initialize the common part */
  _init_cell_sys_bc(face_bc, cm, csys);

  const int  d = eqp->dim;

  /* Identify which face is a boundary face */
  for (short int f = 0; f < cm->n_fc; f++) {
    if (csys->bf_ids[f] > -1) { /* This a boundary face */

      switch(csys->bf_flag[f]) {

      case CS_CDO_BC_HMG_DIRICHLET:
        csys->has_dirichlet = true;
        for (int k = 0; k < d; k++)
          csys->dof_flag[d*f + k] |= CS_CDO_BC_HMG_DIRICHLET;
        break;

      case CS_CDO_BC_DIRICHLET:
        csys->has_dirichlet = true;
        for (int k = 0; k < d; k++) {
          csys->dof_flag[d*f + k] |= CS_CDO_BC_DIRICHLET;
          csys->dir_values[d*f + k] = dir_values[d*csys->bf_ids[f] + k];
        }
        break;

      case CS_CDO_BC_NEUMANN:
        csys->has_nhmg_neumann = true;
        for (int k = 0; k < d; k++)
          csys->dof_flag[d*f + k] |= CS_CDO_BC_NEUMANN;

        cs_equation_compute_neumann_fb(t_eval,
                                       face_bc->def_ids[csys->bf_ids[f]],
                                       f,
                                       eqp,
                                       cm,
                                       csys->neu_values);
        break;

      case CS_CDO_BC_ROBIN:
        csys->has_robin = true;
        for (int k = 0; k < d; k++)
          csys->dof_flag[d*f + k] |= CS_CDO_BC_ROBIN;

        cs_equation_compute_robin(t_eval,
                                  face_bc->def_ids[csys->bf_ids[f]],
                                  f,
                                  eqp,
                                  cm,
                                  csys->rob_values);
        break;

      case CS_CDO_BC_SLIDING:
        csys->has_sliding = true;
        break;

      default:
        /* Nothing to do for instance in case of homogeneous Neumann */
        break;

      } /* End of switch */

    } /* Boundary face */
  } /* Loop on cell faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs when DoFs are attached to
 *          vertices
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      face_bc     pointer to a cs_cdo_bc_face_t structure
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 * \param[in, out] bcflag      pointer to an array storing type of BC
 * \param[in, out] values      pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_dirichlet_vb(cs_real_t                   t_eval,
                                 const cs_mesh_t            *mesh,
                                 const cs_cdo_quantities_t  *quant,
                                 const cs_cdo_connect_t     *connect,
                                 const cs_equation_param_t  *eqp,
                                 const cs_cdo_bc_face_t     *face_bc,
                                 cs_cell_builder_t          *cb,
                                 cs_flag_t                  *bcflag,
                                 cs_real_t                  *values)
{
  assert(face_bc != NULL && bcflag != NULL && values != NULL);

  const cs_lnum_t  *bf2v_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t  *bf2v_lst = mesh->b_face_vtx_lst;

  /* Initialization */
  cs_real_t  *bcvals = cs_equation_get_tmpbuf();
  memset(bcvals, 0, eqp->dim*quant->n_vertices*sizeof(cs_real_t));

  /* Number of faces with a Dir. related to a vertex */
  int  *counter = NULL;
  BFT_MALLOC(counter, quant->n_vertices, int);
  memset(counter, 0, quant->n_vertices*sizeof(int));

  if (face_bc->is_steady == false) /* Update bcflag if needed */
    cs_equation_set_vertex_bc_flag(connect, face_bc, bcflag);

  /* Define array storing the Dirichlet values */
  for (cs_lnum_t i = 0; i < face_bc->n_nhmg_dir_faces; i++) {

    const cs_lnum_t  bf_id = face_bc->nhmg_dir_ids[i];
    const short int  def_id = face_bc->def_ids[bf_id];
    const cs_xdef_t  *def = eqp->bc_defs[def_id];
    const cs_lnum_t  *idx = bf2v_idx + bf_id;
    const cs_lnum_t  *lst = bf2v_lst + idx[0];
    const int  n_vf = idx[1] - idx[0];

    switch(def->type) {

    case CS_XDEF_BY_VALUE:
      _assign_vb_dirichlet_values(eqp->dim, n_vf, lst,
                                  (const cs_real_t *)def->input,
                                  true, /* is constant for all vertices ? */
                                  bcvals,
                                  counter);
      break;

    case CS_XDEF_BY_ARRAY:
      {
        cs_real_t  *eval = cb->values;

        /* Evaluate the boundary condition at each boundary vertex */
        cs_xdef_eval_at_vertices_by_array(n_vf,
                                          lst,
                                          true, /* compact ouput */
                                          mesh,
                                          connect,
                                          quant,
                                          t_eval,
                                          def->input,
                                          eval);

        _assign_vb_dirichlet_values(eqp->dim, n_vf, lst,
                                    eval,
                                    false, /* is constant for all vertices ? */
                                    bcvals,
                                    counter);
      }
      break; /* By array */

    case CS_XDEF_BY_ANALYTIC_FUNCTION:
      {
        cs_real_t  *eval = cb->values;

        /* Evaluate the boundary condition at each boundary vertex */
        cs_xdef_eval_at_vertices_by_analytic(n_vf,
                                             lst,
                                             true, /* compact output */
                                             mesh,
                                             connect,
                                             quant,
                                             t_eval,
                                             def->input,
                                             eval);

        _assign_vb_dirichlet_values(eqp->dim, n_vf, lst,
                                    eval,
                                    false, /* is constant for all vertices ? */
                                    bcvals,
                                    counter);
      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid type of definition.\n"
                  " Stop computing the Dirichlet value.\n"), __func__);

    } /* End of switch: def_type */

  } /* Loop on faces with a non-homogeneous Dirichlet BC */

  cs_equation_sync_vertex_mean_values(connect, eqp->dim, counter, bcvals);

  /* Homogeneous Dirichlet are always enforced (even in case of multiple BCs).
     If multi-valued Dirichlet BCs are set, a weighted sum is used to set the
     Dirichlet value at each corresponding vertex */
  if (eqp->dim == 1) {

#   pragma omp parallel if (quant->n_vertices > CS_THR_MIN)
    {
#     pragma omp for
      for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {
        if (bcflag[v_id] & CS_CDO_BC_HMG_DIRICHLET)
          bcvals[v_id] = 0.;
      }

      /* BC value overwrites the initial value */
#     pragma omp for
      for (cs_lnum_t v = 0; v < quant->n_vertices; v++)
        if (cs_cdo_bc_is_dirichlet(bcflag[v]))
          values[v] = bcvals[v];
    }

  }
  else { /* eqp->dim > 1 */

#   pragma omp parallel if (quant->n_vertices > CS_THR_MIN)
    {
#     pragma omp for
      for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {
        if (bcflag[v_id] & CS_CDO_BC_HMG_DIRICHLET)
          memset(bcvals + eqp->dim*v_id, 0, eqp->dim*sizeof(cs_real_t));
      }

      /* BC value overwrites the initial value */
#     pragma omp for
      for (cs_lnum_t v = 0; v < quant->n_vertices; v++) {
        if (cs_cdo_bc_is_dirichlet(bcflag[v])) {
          for (int k = 0; k < 3; k++)
            values[3*v+k] = bcvals[3*v+k];
        }
      }

    }

  } /* eqp->dim ? */

  /* Free temporary buffers */
  BFT_FREE(counter);

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_BC_DBG > 1
  cs_dbg_darray_to_listing("DIRICHLET_VALUES",
                           eqp->dim*quant->n_vertices, bcvals, 6*eqp->dim);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs when DoFs are attached to
 *          CDO face-based schemes
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      eqp        pointer to a cs_equation_param_t
 * \param[in]      face_bc    pointer to a cs_cdo_bc_face_t structure
 * \param[in]      t_eval     time at which one evaluates the boundary cond.
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_dirichlet_fb(const cs_mesh_t            *mesh,
                                 const cs_cdo_quantities_t  *quant,
                                 const cs_cdo_connect_t     *connect,
                                 const cs_equation_param_t  *eqp,
                                 const cs_cdo_bc_face_t     *face_bc,
                                 cs_real_t                   t_eval,
                                 cs_cell_builder_t          *cb,
                                 cs_real_t                  *values)
{
  assert(face_bc != NULL);

  CS_UNUSED(cb);

  /* Define the array storing the Dirichlet values */
  for (int def_id = 0; def_id < eqp->n_bc_defs; def_id++) {

    const cs_xdef_t  *def = eqp->bc_defs[def_id];

    if (def->meta & CS_CDO_BC_DIRICHLET) {

      assert(eqp->dim == def->dim);

      const cs_zone_t  *bz = cs_boundary_zone_by_id(def->z_id);
      const cs_lnum_t  *elt_ids = bz->elt_ids;

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        {
          const cs_real_t  *constant_val = (cs_real_t *)def->input;

          if (def->dim ==  1) {

#           pragma omp parallel for if (bz->n_elts > CS_THR_MIN)
            for (cs_lnum_t i = 0; i < bz->n_elts; i++) {
              const cs_lnum_t  elt_id = (elt_ids == NULL) ? i : elt_ids[i];
              values[elt_id] = constant_val[0];
            }

          }
          else {

#           pragma omp parallel for if (bz->n_elts > CS_THR_MIN)
            for (cs_lnum_t i = 0; i < bz->n_elts; i++) {
              const cs_lnum_t  elt_id = (elt_ids == NULL) ? i : elt_ids[i];
              for (int k = 0; k < def->dim; k++)
                values[def->dim*elt_id+k] = constant_val[k];
            }

          }
        }
        break;

      case CS_XDEF_BY_ARRAY:
        {
          cs_xdef_array_input_t  *array_input
            = (cs_xdef_array_input_t *)def->input;

          assert(eqp->n_bc_defs == 1); /* Only one definition allowed */
          assert(bz->n_elts == quant->n_b_faces);
          assert(array_input->stride == eqp->dim);
          assert(cs_flag_test(array_input->loc, cs_flag_primal_face));

          memcpy(values, array_input->values,
                 sizeof(cs_real_t)*bz->n_elts*eqp->dim);

        }
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        /* Evaluate the boundary condition at each boundary face */
        switch(eqp->dof_reduction) {

        case CS_PARAM_REDUCTION_DERHAM:
          cs_xdef_eval_at_b_faces_by_analytic(bz->n_elts,
                                              bz->elt_ids,
                                              false, /* compact output */
                                              mesh,
                                              connect,
                                              quant,
                                              t_eval,
                                              def->input,
                                              values);
          break;

        case CS_PARAM_REDUCTION_AVERAGE:
          cs_xdef_eval_avg_at_b_faces_by_analytic(bz->n_elts,
                                                  bz->elt_ids,
                                                  false, /* compact output */
                                                  mesh,
                                                  connect,
                                                  quant,
                                                  t_eval,
                                                  def->input,
                                                  def->qtype,
                                                  def->dim,
                                                  values);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    _(" %s: Invalid type of reduction.\n"
                      " Stop computing the Dirichlet value.\n"), __func__);

        } /* switch on reduction */
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: Invalid type of definition.\n"
                    " Stop computing the Dirichlet value.\n"), __func__);

      } /* switch on def_type */

    } /* Definition based on Dirichlet BC */
  } /* Loop on definitions */

  /* Set the values to zero for all faces attached to a homogeneous
     Diriclet BC */
# pragma omp parallel for if (quant->n_b_faces > CS_THR_MIN)
  for (cs_lnum_t bf_id = 0; bf_id < quant->n_b_faces; bf_id++)
    if (face_bc->flag[bf_id] & CS_CDO_BC_HMG_DIRICHLET)
      for (int k = 0; k < eqp->dim; k++)
        values[eqp->dim*bf_id + k] = 0;

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_BC_DBG > 1
  cs_dbg_darray_to_listing("DIRICHLET_VALUES",
                           eqp->dim*quant->n_b_faces, values, 9);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define an array of flags for each vertex collecting the flags
 *          of associated boundary faces
 *
 * \param[in]      connect   pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]      face_bc   pointer to a structure collecting boundary
 *                           conditions applied to faces
 * \param[in, out] vflag     BC flag on vertices to define
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_vertex_bc_flag(const cs_cdo_connect_t     *connect,
                               const cs_cdo_bc_face_t     *face_bc,
                               cs_flag_t                  *vflag)
{
  if (vflag == NULL)
    return;

  assert(connect->bf2v != NULL);

  const cs_adjacency_t  *bf2v = connect->bf2v;
  const cs_lnum_t  n_vertices = connect->n_vertices;
  const cs_lnum_t  n_b_faces = connect->n_faces[CS_BND_FACES];

  /* Initialization */
  memset(vflag, 0, n_vertices*sizeof(cs_flag_t));

  for (cs_lnum_t bf_id = 0; bf_id < n_b_faces; bf_id++) {

    const cs_flag_t  bc_flag = face_bc->flag[bf_id];
    for (cs_lnum_t j = bf2v->idx[bf_id]; j < bf2v->idx[bf_id+1]; j++) {
      vflag[bf2v->ids[j]] |= bc_flag;
    }

  } /* Loop on border faces */

#if defined(DEBUG) && !defined(NDEBUG)
  for (cs_lnum_t bf_id = 0; bf_id < n_b_faces; bf_id++) {
    for (cs_lnum_t j = bf2v->idx[bf_id]; j < bf2v->idx[bf_id+1]; j++) {
      const cs_lnum_t v_id = bf2v->ids[j];
      if (vflag[v_id] == 0)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Border vertices %d without any boundary conditions.",
                  __func__, v_id);
    }
  } /* Loop on border faces */
#endif

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    assert(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL] != NULL);

    cs_interface_set_max(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                         n_vertices,
                         1,             /* stride */
                         false,         /* interlace (not useful here) */
                         CS_FLAG_TYPE,  /* unsigned short int */
                         vflag);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are scalar-valued
 *          and attached to vertices.
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing the Neumann values
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_sv(cs_real_t                   t_eval,
                               short int                   def_id,
                               short int                   f,
                               const cs_equation_param_t  *eqp,
                               const cs_cell_mesh_t       *cm,
                               double                     *neu_values)
{
  assert(neu_values != NULL && cm != NULL && eqp != NULL);
  assert(def_id > -1);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_EV | CS_FLAG_COMP_FE | CS_FLAG_COMP_FV));

  const cs_xdef_t  *def = eqp->bc_defs[def_id];

  assert(def->dim == 3); /* flux is a vector in the scalar case */
  assert(def->meta & CS_CDO_BC_NEUMANN); /* Neuman BC */

  /* Evaluate the boundary condition at each boundary vertex */
  switch(def->type) {

  case CS_XDEF_BY_VALUE:
    cs_xdef_cw_eval_flux_at_vtx_by_val(cm, f, t_eval, def->input, neu_values);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_xdef_cw_eval_flux_at_vtx_by_analytic(cm,
                                            f,
                                            t_eval,
                                            def->input,
                                            def->qtype,
                                            neu_values);
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *array_input
        = (cs_xdef_array_input_t *)def->input;

      assert(eqp->n_bc_defs == 1); /* Only one definition allowed */
      assert(array_input->stride == 3);

      cs_lnum_t  bf_id = cm->f_ids[f] - cm->bface_shift;
      assert(bf_id > -1);

      if (cs_flag_test(array_input->loc, cs_flag_primal_face))
        cs_xdef_cw_eval_flux_at_vtx_by_val(cm, f, t_eval,
                                           array_input->values + 3*bf_id,
                                           neu_values);

      else if (cs_flag_test(array_input->loc, cs_flag_dual_closure_byf)) {

        assert(array_input->index != NULL);

        /* Retrieve the bf2v->idx stored in the cs_cdo_connect_t structure */
        cs_lnum_t  shift = array_input->index[bf_id];
        for (short int iv = cm->f2v_idx[f]; iv < cm->f2v_idx[f+1]; iv++)
          neu_values[cm->f2v_ids[iv]] = array_input->values[shift++];

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid array location.", __func__);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of definition.\n"
                " Stop computing the Neumann value.\n"));

  } /* switch def_type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are attached to
 *          faces.
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] neu_values  array storing Neumann values to use
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_fb(cs_real_t                    t_eval,
                               short int                    def_id,
                               short int                    f,
                               const cs_equation_param_t   *eqp,
                               const cs_cell_mesh_t        *cm,
                               double                      *neu_values)
{
  assert(neu_values != NULL && cm != NULL && eqp != NULL);
  assert(def_id > -1);
  assert(eqp->dim == 1 || eqp->dim == 3);

  const cs_xdef_t  *def = eqp->bc_defs[def_id];

  /* Flux is a vector in the scalar-valued case and a tensor in the
     vector-valued case */
  assert(def->meta & CS_CDO_BC_NEUMANN); /* Neuman BC */

  /* Evaluate the boundary condition at each boundary face */
  switch(def->type) {

  case CS_XDEF_BY_VALUE:
    if (eqp->dim == 1)
      cs_xdef_cw_eval_flux_by_val(cm, f, t_eval, def->input, neu_values);
    else if (eqp->dim == 3)
      cs_xdef_cw_eval_tensor_flux_by_val(cm, f, t_eval, def->input, neu_values);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    if (eqp->dim == 1)
      cs_xdef_cw_eval_flux_by_analytic(cm,
                                       f,
                                       t_eval,
                                       def->input,
                                       def->qtype,
                                       neu_values);
    else if (eqp->dim == 3)
      cs_xdef_cw_eval_tensor_flux_by_analytic(cm,
                                              f,
                                              t_eval,
                                              def->input,
                                              def->qtype,
                                              neu_values);
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *array_input
        = (cs_xdef_array_input_t *)def->input;

      assert(eqp->n_bc_defs == 1); /* Only one definition allowed */
      assert(array_input->stride == 3);
      assert(cs_flag_test(array_input->loc, cs_flag_primal_face));

      cs_lnum_t  bf_id = cm->f_ids[f] - cm->bface_shift;
      assert(bf_id > -1);

      cs_real_t  *face_val = array_input->values + 3*bf_id;

      cs_xdef_cw_eval_flux_by_val(cm, f, t_eval, face_val, neu_values);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of definition.\n"
                " Stop computing the Neumann value.\n"));

  } /* switch def_type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Robin BCs
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in, out] rob_values  array storing Robin values to use
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_robin(cs_real_t                    t_eval,
                          short int                    def_id,
                          short int                    f,
                          const cs_equation_param_t   *eqp,
                          const cs_cell_mesh_t        *cm,
                          double                      *rob_values)
{
  assert(rob_values != NULL && cm != NULL && eqp != NULL);
  assert(def_id > -1);
  assert(eqp->dim == 1);

  const cs_xdef_t  *def = eqp->bc_defs[def_id];

  /* Flux is a vector in the scalar-valued case and a tensor in the
     vector-valued case */
  assert(def->meta & CS_CDO_BC_ROBIN); /* Robin BC */

  /* Evaluate the boundary condition at each boundary face */
  switch(def->type) {

  case CS_XDEF_BY_VALUE:
    {
      const cs_real_t  *parameters = (cs_real_t *)def->input;

      rob_values[3*f  ] = parameters[0];
      rob_values[3*f+1] = parameters[1];
      rob_values[3*f+2] = parameters[2];
    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_real_3_t  parameters = {0, 0, 0};
      cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)def->input;

      anai->func(t_eval, 1, NULL,
                 cm->face[f].center,
                 true, /* compacted output ? */
                 anai->input,
                 (cs_real_t *)parameters);

      rob_values[3*f  ] = parameters[0];
      rob_values[3*f+1] = parameters[1];
      rob_values[3*f+2] = parameters[2];
    }
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *array_input
        = (cs_xdef_array_input_t *)def->input;

      assert(eqp->n_bc_defs == 1); /* Only one definition allowed */
      assert(array_input->stride == 3);
      assert(cs_flag_test(array_input->loc, cs_flag_primal_face));

      cs_lnum_t  bf_id = cm->f_ids[f] - cm->bface_shift;
      assert(bf_id > -1);
      cs_real_t  *parameters = array_input->values + 3*bf_id;

      rob_values[3*f  ] = parameters[0];
      rob_values[3*f+1] = parameters[1];
      rob_values[3*f+2] = parameters[2];
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of definition.\n"
                " Stop computing the Neumann value.\n"));

  } /* switch def_type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the tangential component lying on the domain
 *          boundary. Kind of BCs used when DoFs are attached to CDO (primal)
 *          edge-based schemes. One sets the values of the circulation.
 *
 * \param[in]      t_eval     time at which one evaluates the boundary cond.
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      eqp        pointer to a cs_equation_param_t
 * \param[in, out] values     pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_circulation_eb(cs_real_t                    t_eval,
                                   const cs_mesh_t             *mesh,
                                   const cs_cdo_quantities_t   *quant,
                                   const cs_cdo_connect_t      *connect,
                                   const cs_equation_param_t   *eqp,
                                   cs_real_t                   *values)
{
  assert(values != NULL);

  /* Synchronization of the definition of the circulation if needed */
  cs_lnum_t  *def2e_ids = (cs_lnum_t *)cs_equation_get_tmpbuf();
  cs_lnum_t  *def2e_idx = NULL;
  BFT_MALLOC(def2e_idx, eqp->n_bc_defs + 1, cs_lnum_t);

  _sync_circulation_def_at_edges(connect,
                                 eqp->n_bc_defs,
                                 eqp->bc_defs,
                                 def2e_idx,
                                 def2e_ids);

  /* Define the array storing the circulation values */
  for (int def_id = 0; def_id < eqp->n_bc_defs; def_id++) {

    const cs_xdef_t  *def = eqp->bc_defs[def_id];

    if (def->meta & CS_CDO_BC_TANGENTIAL_DIRICHLET) {

      const cs_lnum_t  n_elts = def2e_idx[def_id+1] - def2e_idx[def_id];
      const cs_lnum_t  *elt_ids = def2e_ids + def2e_idx[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_circulation_along_edges_by_value(def,
                                                     n_elts,
                                                     elt_ids,
                                                     values);
        break;

      case CS_XDEF_BY_ARRAY:
        cs_evaluate_circulation_along_edges_by_array(def,
                                                     n_elts,
                                                     elt_ids,
                                                     values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_evaluate_circulation_along_edges_by_analytic(def,
                                                        t_eval,
                                                        n_elts,
                                                        elt_ids,
                                                        values);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: Invalid type of definition.\n"
                    " Stop computing the circulation.\n"), __func__);

      } /* Switch on def_type */

    } /* Definition related to a circulation */

  } /* Loop on definitions */

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_BC_DBG > 1
  cs_dbg_darray_to_listing("CIRCULATION_VALUES", quant->n_edges, values, 9);
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
