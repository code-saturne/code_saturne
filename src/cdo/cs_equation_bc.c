/*============================================================================
 * Routines to handle the evaluation of boundary conditions when building the
 * algebraic system in CDO/HHO schemes
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
#include "cs_xdef.h"

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
 * \param[in, out] flag         flag to update
 * \param[in, out] counter      counter to update
 */
/*----------------------------------------------------------------------------*/

static inline void
_assign_vb_dirichlet_values(int                dim,
                            int                n_vf,
                            const cs_lnum_t   *lst,
                            const cs_real_t   *eval,
                            _Bool              is_constant,
                            cs_real_t         *vvals,
                            cs_flag_t          flag[],
                            cs_flag_t          counter[])
{
  switch (dim) {

  case 1:
    for (short int v = 0; v < n_vf; v++) {

      const cs_lnum_t  v_id = lst[v];
      const short int  _v = is_constant ? 0 : v;

      flag[v_id] |= CS_CDO_BC_DIRICHLET;
      counter[v_id] += 1;
      vvals[v_id] += eval[_v];
    }
    break;

  default:
    for (short int v = 0; v < n_vf; v++) {

      const cs_lnum_t  v_id = lst[v];
      const short int  _v = is_constant ? 0 : v;

      flag[v_id] |= CS_CDO_BC_DIRICHLET;
      counter[v_id] += 1;
      for (int k = 0; k < dim; k++)
        vvals[dim*v_id + k] += eval[dim*_v + k];
    }
    break;

  } /* End of switch */
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the values for the normal boundary flux stemming from the
 *         Neumann boundary conditions (zero is left where a Dirichlet is
 *         set. This can be updated later one)
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
 * \param[in]      bf_id       id of the border face
 * \param[in]      f           id of the current face in a cellwise numbering
 * \param[in]      face_flag   metadata about the current face
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      dir_values  Dirichlet values associated to each vertex
 * \param[in]      neu_tags    definition id related to each Neumann face
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_vb_set_cell_bc(cs_lnum_t                     bf_id,
                           short int                     f,
                           cs_flag_t                     face_flag,
                           const cs_cell_mesh_t         *cm,
                           const cs_cdo_connect_t       *connect,
                           const cs_cdo_quantities_t    *quant,
                           const cs_equation_param_t    *eqp,
                           const cs_real_t               dir_values[],
                           const short int               neu_tags[],
                           cs_real_t                     t_eval,
                           cs_cell_sys_t                *csys,
                           cs_cell_builder_t            *cb)
{
  CS_UNUSED(connect);

  /* Sanity check */
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE));

  /* Identify which face is a boundary face */
  short int  n_vf;

  csys->bf_flag[csys->n_bc_faces] = face_flag;
  csys->bf_ids[csys->n_bc_faces] = bf_id;
  csys->_f_ids[csys->n_bc_faces++] = f;

  cs_cell_mesh_get_f2v(f, cm, &n_vf, cb->ids);
  assert(n_vf == cm->f2e_idx[f+1] - cm->f2e_idx[f]);

  if (face_flag & CS_CDO_BC_HMG_DIRICHLET) {

    csys->has_dirichlet = true;
    for (short int i = 0; i < n_vf; i++)
      csys->dof_flag[cb->ids[i]] |= CS_CDO_BC_HMG_DIRICHLET;

  }
  else if (face_flag & CS_CDO_BC_DIRICHLET) {

    csys->has_dirichlet = true;
    for (short int i = 0; i < n_vf; i++) {
      short int  v = cb->ids[i];
      csys->dir_values[v] = dir_values[cm->v_ids[v]];
      csys->dof_flag[v] |= CS_CDO_BC_DIRICHLET;
    }

  }
  else if (face_flag & CS_CDO_BC_NEUMANN) {

    csys->has_nhmg_neumann = true;
    for (short int i = 0; i < n_vf; i++)
      csys->dof_flag[cb->ids[i]] |= CS_CDO_BC_NEUMANN;

    cs_equation_compute_neumann_sv(neu_tags[bf_id],
                                   f,
                                   quant,
                                   eqp,
                                   cm,
                                   cb->ids,
                                   t_eval,
                                   csys->neu_values);

  }
  else if (face_flag & CS_CDO_BC_ROBIN) {
    csys->has_robin = true;
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the BC into a cellwise view of the current system.
 *          Case of Face-based schemes
 *
 * \param[in]      bf_id       id of the border face
 * \param[in]      f           id of the current face in a cellwise numbering
 * \param[in]      face_flag   metadata about the current face
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      dir_values  Dirichlet values associated to each vertex
 * \param[in]      neu_tags    definition id related to each Neumann face
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_fb_set_cell_bc(cs_lnum_t                     bf_id,
                           short int                     f,
                           cs_flag_t                     face_flag,
                           const cs_cell_mesh_t         *cm,
                           const cs_cdo_connect_t       *connect,
                           const cs_cdo_quantities_t    *quant,
                           const cs_equation_param_t    *eqp,
                           const cs_real_t               dir_values[],
                           const short int               neu_tags[],
                           cs_real_t                     t_eval,
                           cs_cell_sys_t                *csys,
                           cs_cell_builder_t            *cb)
{
  CS_UNUSED(connect);
  CS_UNUSED(cb);

  csys->bf_flag[csys->n_bc_faces] = face_flag;
  csys->bf_ids[csys->n_bc_faces] = bf_id;
  csys->_f_ids[csys->n_bc_faces++] = f;

  if (face_flag & CS_CDO_BC_HMG_DIRICHLET) {

    csys->has_dirichlet = true;
    for (int k = 0; k < eqp->dim; k++)
      csys->dof_flag[eqp->dim*f + k] |= CS_CDO_BC_HMG_DIRICHLET;

  }
  else if (face_flag & CS_CDO_BC_DIRICHLET) {

    csys->has_dirichlet = true;
    for (int k = 0; k < eqp->dim; k++) {
      csys->dof_flag[eqp->dim*f + k] |= CS_CDO_BC_DIRICHLET;
      csys->dir_values[eqp->dim*f + k] = dir_values[eqp->dim*bf_id + k];
    }

  }
  else if (face_flag & CS_CDO_BC_NEUMANN) {

    csys->has_nhmg_neumann = true;
    for (int k = 0; k < eqp->dim; k++)
      csys->dof_flag[eqp->dim*f + k] |= CS_CDO_BC_NEUMANN;

    cs_equation_compute_neumann_fb(neu_tags[bf_id], f,
                                   quant,
                                   eqp, cm,
                                   t_eval,
                                   csys->neu_values);

  }
  else if (face_flag & CS_CDO_BC_ROBIN) {
    csys->has_robin = true;
    /* TODO */
    bft_error(__FILE__, __LINE__, 0, "%s: TODO", __func__);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs when DoFs are attached to
 *          vertices
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      face_bc     pointer to a cs_cdo_bc_t structure
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 * \param[in, out] values      pointer to the array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_dirichlet_vb(const cs_mesh_t            *mesh,
                                 const cs_cdo_quantities_t  *quant,
                                 const cs_cdo_connect_t     *connect,
                                 const cs_equation_param_t  *eqp,
                                 const cs_cdo_bc_t          *face_bc,
                                 cs_real_t                   t_eval,
                                 cs_cell_builder_t          *cb,
                                 cs_real_t                  *values)
{
  assert(face_bc != NULL);
  assert(face_bc->n_elts == quant->n_b_faces);

  const cs_cdo_bc_list_t  *dir = face_bc->dir;
  const cs_lnum_t  *bf2v_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t  *bf2v_lst = mesh->b_face_vtx_lst;

  cs_flag_t  *flag = NULL, *counter = NULL;

  /* Initialization */
  BFT_MALLOC(counter, quant->n_vertices, cs_flag_t);
  BFT_MALLOC(flag, quant->n_vertices, cs_flag_t);

# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {
    flag[v_id] = 0;    /* No flag by default */
    counter[v_id] = 0; /* Number of faces with a Dir. related to a vertex */
  }

  /* Set the values to zero for all vertices attached to a Diriclet BC */
# pragma omp parallel for if (quant->n_b_faces > CS_THR_MIN)
  for (cs_lnum_t bf_id = 0; bf_id < quant->n_b_faces; bf_id++) {
    if ((face_bc->flag[bf_id] & CS_CDO_BC_DIRICHLET) ||
        (face_bc->flag[bf_id] & CS_CDO_BC_HMG_DIRICHLET)) {
      for (cs_lnum_t i = bf2v_idx[bf_id]; i < bf2v_idx[bf_id+1]; i++) {

        const cs_lnum_t v_id = bf2v_lst[i];
        for (int k = 0; k < eqp->dim; k++)
          values[eqp->dim*v_id + k] = 0;

      }
    } /* face with a Dirichlet */
  } /* Loop on boundary faces */

  /* Define array storing the Dirichlet values */
  for (cs_lnum_t i = 0; i < dir->n_nhmg_elts; i++) {

    const cs_lnum_t  bf_id = dir->elt_ids[i];
    const cs_lnum_t  *idx = bf2v_idx + bf_id;
    const cs_lnum_t  *lst = bf2v_lst + idx[0];
    const int  n_vf = idx[1] - idx[0];
    const short int  def_id = dir->def_ids[i];
    const cs_xdef_t  *def = eqp->bc_defs[def_id];

    switch(def->type) {

    case CS_XDEF_BY_VALUE:
      _assign_vb_dirichlet_values(eqp->dim, n_vf, lst,
                                  (const cs_real_t *)def->input,
                                  true, /* is constant for all vertices ? */
                                  values, flag, counter);
      break;

    case CS_XDEF_BY_ARRAY:
      {
        cs_real_t  *eval = cb->values;

        /* Evaluate the boundary condition at each boundary vertex */
        cs_xdef_eval_at_vertices_by_array(n_vf,
                                          lst,
                                          true, // compact ouput
                                          mesh,
                                          connect,
                                          quant,
                                          t_eval,
                                          def->input,
                                          eval);

        _assign_vb_dirichlet_values(eqp->dim, n_vf, lst,
                                    eval,
                                    false, /* is constant for all vertices ? */
                                    values, flag, counter);

      }
      break; /* By array */

    case CS_XDEF_BY_ANALYTIC_FUNCTION:
      {
        cs_real_t  *eval = cb->values;

        /* Evaluate the boundary condition at each boundary vertex */
        cs_xdef_eval_at_vertices_by_analytic(n_vf,
                                             lst,
                                             true, // compact output
                                             mesh,
                                             connect,
                                             quant,
                                             t_eval,
                                             def->input,
                                             eval);

        _assign_vb_dirichlet_values(eqp->dim, n_vf, lst,
                                    eval,
                                    false, /* is constant for all vertices ? */
                                    values, flag, counter);
      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid type of definition.\n"
                  " Stop computing the Dirichlet value.\n"), __func__);

    } /* End of switch: def_type */

  } /* Loop on faces with a non-homogeneous Dirichlet BC */

  /* Define array storing the Dirichlet values */
  for (cs_lnum_t i = dir->n_nhmg_elts; i < dir->n_elts; i++) {

    const cs_lnum_t  f_id = dir->elt_ids[i];
    const cs_lnum_t  *idx = bf2v_idx + f_id;
    const int  n_vf = idx[1] - idx[0];
    const cs_lnum_t  *lst = bf2v_lst + idx[0];

    for (short int v = 0; v < n_vf; v++)
      flag[lst[v]] |= CS_CDO_BC_HMG_DIRICHLET;

  } /* Loop on faces with a non-homogeneous Dirichlet BC */

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    cs_interface_set_max(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                         quant->n_vertices,
                         1,            // stride
                         false,        // interlace (not useful here)
                         CS_FLAG_TYPE, // unsigned short int
                         flag);

    cs_interface_set_sum(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                         quant->n_vertices,
                         1,            // stride
                         false,        // interlace (not useful here)
                         CS_FLAG_TYPE, // unsigned short int
                         counter);

    cs_interface_set_sum(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                         quant->n_vertices,
                         eqp->dim,     // stride
                         false,        // interlace (not useful here)
                         CS_REAL_TYPE,
                         values);

  }

  /* Homogeneous Dirichlet are always enforced (even in case of multiple BCs).
     If multi-valued Dirichlet BCs are set, a weighted sum is used to set the
     Dirichlet value at each corresponding vertex */
  if (eqp->dim == 1) {

#   pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {

      if (flag[v_id] & CS_CDO_BC_HMG_DIRICHLET)
        values[v_id] = 0.;
      else if (flag[v_id] & CS_CDO_BC_DIRICHLET) {
        assert(counter[v_id] > 0);
        if (counter[v_id] > 1)
          values[v_id] /= counter[v_id];
      }

    } /* Loop on vertices */

  }
  else { /* eqp->dim > 1 */

#   pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {

      if (flag[v_id] & CS_CDO_BC_HMG_DIRICHLET) {
        for (int j = 0; j < eqp->dim; j++)
          values[eqp->dim*v_id + j] = 0.;
      }
      else if (flag[v_id] & CS_CDO_BC_DIRICHLET) {
        assert(counter[v_id] > 0);
        if (counter[v_id] > 1) {
          const cs_real_t  inv_count = counter[v_id];
          for (int j = 0; j < eqp->dim; j++)
            values[eqp->dim*v_id + j] *= inv_count;
        }
      }

    } /* Loop on vertices */

  } /* eqp->dim ? */

  /* Free temporary buffers */
  BFT_FREE(counter);
  BFT_FREE(flag);

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_BC_DBG > 1
  cs_dbg_darray_to_listing("DIRICHLET_VALUES",
                           eqp->dim*quant->n_vertices, values, 6*eqp->dim);
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
 * \param[in]      face_bc    pointer to a cs_cdo_bc_t structure
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
                                 const cs_cdo_bc_t          *face_bc,
                                 cs_real_t                   t_eval,
                                 cs_cell_builder_t          *cb,
                                 cs_real_t                  *values)
{
  assert(face_bc != NULL);
  assert(face_bc->n_elts == quant->n_b_faces);

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
          cs_xdef_array_input_t  *array_input =
            (cs_xdef_array_input_t *)def->input;

          assert(eqp->n_bc_defs == 1); // Only one definition allowed
          assert(bz->n_elts == quant->n_b_faces);
          assert(array_input->stride == eqp->dim);
          assert(cs_flag_test(array_input->loc, cs_flag_primal_face));

          memcpy(values, array_input->values,
                 sizeof(cs_real_t)*bz->n_elts*eqp->dim);

        }
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        /* Evaluate the boundary condition at each boundary vertex */
        cs_xdef_eval_at_b_faces_by_analytic(bz->n_elts,
                                            bz->elt_ids,
                                            false, // compact output
                                            mesh,
                                            connect,
                                            quant,
                                            t_eval,
                                            def->input,
                                            values);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Invalid type of definition.\n"
                    " Stop computing the Dirichlet value.\n"));

      } // switch def_type

    } // Definition based on Dirichlet BC
  } // Loop on definitions

  /* Set the values to zero for all vertices attached to a homogeneous
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
 * \brief   Tag each face related to a Neumann BC with its definition id.
 *          Default tag is -1 (not a Neumann face)
 *
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp        pointer to a cs_equation_param_t

 * \return an array with prescribed tags
 */
/*----------------------------------------------------------------------------*/

short int *
cs_equation_tag_neumann_face(const cs_cdo_quantities_t    *quant,
                             const cs_equation_param_t    *eqp)
{
  short int  *face_tag = NULL;

  /* Initialization */
  BFT_MALLOC(face_tag, quant->n_b_faces, short int);

# pragma omp parallel for if (quant->n_b_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < quant->n_b_faces; f_id++)
    face_tag[f_id] = -1;

  /* Tag faces with Neumann BCs */
  for (int def_id = 0; def_id < eqp->n_bc_defs; def_id++) {

    const cs_xdef_t  *def = eqp->bc_defs[def_id];
    if (def->meta & CS_CDO_BC_NEUMANN) {

      const cs_zone_t  *bz = cs_boundary_zone_by_id(def->z_id);

#     pragma omp parallel for if (bz->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < bz->n_elts; i++)
        face_tag[bz->elt_ids[i]] = def_id;

    }
  }

  return face_tag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are scalar-valued
 *          and attached to vertices.
 *
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      f2v_ids     vertex ids of the face in the cell numbering
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] neu_values  array storing the Neumann values
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_sv(short int                   def_id,
                               short int                   f,
                               const cs_cdo_quantities_t  *quant,
                               const cs_equation_param_t  *eqp,
                               const cs_cell_mesh_t       *cm,
                               const short int            *f2v_ids,
                               cs_real_t                   t_eval,
                               double                     *neu_values)
{
  assert(neu_values != NULL && cm != NULL && eqp != NULL);
  assert(def_id > -1);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE));

  const cs_xdef_t  *def = eqp->bc_defs[def_id];

  assert(def->dim == 3);                 // flux is a vector in the scalar case
  assert(def->meta & CS_CDO_BC_NEUMANN); // Neuman BC

  /* Evaluate the boundary condition at each boundary vertex */
  switch(def->type) {

  case CS_XDEF_BY_VALUE:
    cs_xdef_eval_cw_at_vtx_flux_by_val(cm, f, t_eval, def->input, neu_values);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_xdef_eval_cw_at_vtx_flux_by_analytic(cm,
                                            f,
                                            t_eval,
                                            def->input,
                                            def->qtype,
                                            neu_values);
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *array_input =
        (cs_xdef_array_input_t *)def->input;

      assert(eqp->n_bc_defs == 1); // Only one definition allowed
      assert(array_input->stride == 3);

      cs_lnum_t  bf_id = cm->f_ids[f] - quant->n_i_faces;
      assert(bf_id > -1);

      if (cs_flag_test(array_input->loc, cs_flag_primal_face))
        cs_xdef_eval_cw_at_vtx_flux_by_val(cm, f, t_eval,
                                           array_input->values + 3*bf_id,
                                           neu_values);

      else if (cs_flag_test(array_input->loc, cs_flag_dual_face_byc)) {

        assert(array_input->index != NULL);

        const short int  n_vf = cm->f2e_idx[f+1] - cm->f2e_idx[f];
        /* Retrieve the bf2v->idx stored in the cs_cdo_connect_t structure */
        const cs_lnum_t  shift = array_input->index[bf_id];
        for (short int v = 0; v < n_vf; v++)
          neu_values[f2v_ids[v]] = array_input->values[shift + v];

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

  } // switch def_type

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are attached to
 *          faces.
 *
 * \param[in]      def_id      id of the definition for setting the Neumann BC
 * \param[in]      f           local face number in the cs_cell_mesh_t
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] neu_values  array storing Neumann values to use
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_fb(short int                    def_id,
                               short int                    f,
                               const cs_cdo_quantities_t   *quant,
                               const cs_equation_param_t   *eqp,
                               const cs_cell_mesh_t        *cm,
                               cs_real_t                    t_eval,
                               double                      *neu_values)
{
  assert(neu_values != NULL && cm != NULL && eqp != NULL);
  assert(def_id > -1);
  assert(eqp->dim == 1 || eqp->dim == 3);

  const cs_xdef_t  *def = eqp->bc_defs[def_id];

  /* Flux is a vector in the scalar-valued case and a tensor in the
     vector-valued case */
  assert(def->meta & CS_CDO_BC_NEUMANN); // Neuman BC

  /* Evaluate the boundary condition at each boundary vertex */
  switch(def->type) {

  case CS_XDEF_BY_VALUE:
    if (eqp->dim == 1)
      cs_xdef_eval_cw_flux_by_val(cm, f, t_eval, def->input, neu_values);
    else if (eqp->dim == 3)
      cs_xdef_eval_cw_tensor_flux_by_val(cm, f, t_eval, def->input, neu_values);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    if (eqp->dim == 1)
      cs_xdef_eval_cw_flux_by_analytic(cm,
                                       f,
                                       t_eval,
                                       def->input,
                                       def->qtype,
                                       neu_values);
    else if (eqp->dim == 3)
      cs_xdef_eval_cw_tensor_flux_by_analytic(cm,
                                              f,
                                              t_eval,
                                              def->input,
                                              def->qtype,
                                              neu_values);
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *array_input =
        (cs_xdef_array_input_t *)def->input;

      assert(eqp->n_bc_defs == 1); // Only one definition allowed
      assert(array_input->stride == 3);
      assert(cs_flag_test(array_input->loc, cs_flag_primal_face));

      cs_lnum_t  bf_id = cm->f_ids[f] - quant->n_i_faces;
      assert(bf_id > -1);

      cs_real_t  *face_val = array_input->values + 3*bf_id;

      cs_xdef_eval_cw_flux_by_val(cm, f, t_eval, face_val, neu_values);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of definition.\n"
                " Stop computing the Neumann value.\n"));

  } // switch def_type

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
