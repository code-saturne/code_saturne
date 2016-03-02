/*============================================================================
 * Manage boundary conditions
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

#include <errno.h>
#include <locale.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_bc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/* Tag to define the underpinning behavior of boundary conditions */
#define CS_BC_FLAG_HMG    (1 <<  0)  /*     1: homogeneous (i.e. 0) */
#define CS_BC_FLAG_CONST  (1 <<  1)  /*     2: constant value */
#define CS_BC_FLAG_DIRI   (1 <<  2)  /*     4: Dirichlet BC*/
#define CS_BC_FLAG_NEUM   (1 <<  3)  /*     8: Neumann BC */
#define CS_BC_FLAG_ROBIN  (1 <<  4)  /*    16: Robin BC*/
#define CS_BC_FLAG_SCAL   (1 <<  5)  /*    32: scalar-valued */
#define CS_BC_FLAG_VECT   (1 <<  6)  /*    64: vector-valued */
#define CS_BC_FLAG_TENS   (1 <<  7)  /*   128: tensor-valued */
#define CS_BC_FLAG_TANG   (1 <<  8)  /*   256: tangential component */
#define CS_BC_FLAG_NORM   (1 <<  9)  /*   512: normal component */
#define CS_BC_FLAG_VTX    (1 << 10)  /*  1024: on vertices */
#define CS_BC_FLAG_EDGE   (1 << 11)  /*  2048: on edges */
#define CS_BC_FLAG_FACE   (1 << 12)  /*  4096: on faces */
#define CS_BC_FLAG_CELL   (1 << 13)  /*  8192: on cells */
#define CS_BC_FLAG_PRIM   (1 << 14)  /* 16384: on primal mesh */
#define CS_BC_FLAG_DUAL   (1 << 15)  /* 32768: on dual mesh */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_cdo_bc_list_t structure
 *
 * \param[in]  n_elts      number of entries of the list
 * \param[in]  n_nhmg_elts number of elements attached to a homogeneous BC
 *
 * \return  a new allocated pointer to a cs_cdo_bc_list_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_list_t *
cs_cdo_bc_list_create(cs_lnum_t   n_elts,
                      cs_lnum_t   n_nhmg_elts)
{
  cs_lnum_t  i;

  cs_cdo_bc_list_t  *bcl = NULL;

  /* Sanity check */
  assert(n_elts >= n_nhmg_elts && n_elts > -1);

  BFT_MALLOC(bcl, 1, cs_cdo_bc_list_t);

  bcl->n_elts = n_elts;
  bcl->n_nhmg_elts = n_nhmg_elts;
  bcl->elt_ids = NULL;
  bcl->def_ids = NULL;

  if (n_elts > 0) {  /* Allocate and initialize by default */

    BFT_MALLOC(bcl->elt_ids, n_elts, cs_lnum_t);
    for (i = 0; i < n_elts; i++)
      bcl->elt_ids[i] = -1;

    BFT_MALLOC(bcl->def_ids, n_nhmg_elts, short int);
    for (i = 0; i < n_nhmg_elts; i++)
      bcl->def_ids[i] = -1;

  }

  return bcl;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_cdo_bc_list_t structure
 *
 * \param[in]  bcl     pointer to the cs_cdo_bc_list_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_list_t *
cs_cdo_bc_list_free(cs_cdo_bc_list_t   *bcl)
{
  if (bcl == NULL)
    return bcl;

  if (bcl->n_elts > 0) {
    BFT_FREE(bcl->def_ids);
    BFT_FREE(bcl->elt_ids);
  }
  BFT_FREE(bcl);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Prepare the treatment of the boundary conditions.
 *         Compile the information detailed in a cs_param_bc_t structure
 *         into the structure cs_cdo_bc_t (based on border faces).
 *         This is a primilary step to be ready to set the values of the BC
 *
 * \param[in] param_bc    pointer to the parameters related to BCs
 * \param[in] n_b_faces   number of border faces
 *
 * \return a pointer to a new allocated cs_cdo_bc_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_t *
cs_cdo_bc_init(const cs_param_bc_t  *param_bc,
               cs_lnum_t             n_b_faces)
{
  cs_lnum_t  i, id, shift, n_ents;
  cs_lnum_t  count[CS_PARAM_N_BC_TYPES];
  cs_param_bc_type_t  type;

  cs_param_bc_type_t  *bc_types = NULL;
  cs_cdo_bc_t  *bc = NULL;

  /* Sanity check */
  assert(param_bc != NULL);

  /* Create and initialize a new structure */
  BFT_MALLOC(bc, 1, cs_cdo_bc_t);

  bc->n_b_faces = n_b_faces;

  bc->dir = NULL;
  bc->neu = NULL;
  bc->rob = NULL;

  if (param_bc->default_bc != CS_PARAM_BC_HMG_DIRICHLET &&
      param_bc->default_bc != CS_PARAM_BC_HMG_NEUMANN)
    bft_error(__FILE__, __LINE__, 0,
              _(" Incompatible type of boundary condition by default.\n"
                " Please modify your settings.\n"));

  if (n_b_faces > 0) {

    // 1) build bc_types
    BFT_MALLOC(bc_types, n_b_faces, cs_param_bc_type_t);

    for (i = 0; i < n_b_faces; i++)
      bc_types[i] = param_bc->default_bc;

    /* Loop on the definition of each boundary condition */
    for (id = 0; id < param_bc->n_defs; id++) {

      const cs_param_bc_def_t  *def = param_bc->defs + id;
      const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(def->loc_id);
      const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(def->loc_id);

      if (elt_ids == NULL) {
        assert(n_elts[0] == n_b_faces);
        for (i = 0; i < n_elts[0]; i++)
          bc_types[i] = def->bc_type;
      }
      else
        for (i = 0; i < n_elts[0]; i++)
          bc_types[elt_ids[i]] = def->bc_type;

    } // Loop on boundary conditions

    // 2) Define bc->type_shift and bc->type_size
    for (i = 0; i < CS_PARAM_N_BC_TYPES; i++)
      count[i] = 0;
    for (i = 0; i < n_b_faces; i++)
      count[bc_types[i]] += 1;

    // 3) Allocate cs_cdo_bc_list_t structures
    n_ents = count[CS_PARAM_BC_DIRICHLET] + count[CS_PARAM_BC_HMG_DIRICHLET];
    bc->dir = cs_cdo_bc_list_create(n_ents, count[CS_PARAM_BC_DIRICHLET]);

    n_ents = count[CS_PARAM_BC_NEUMANN] + count[CS_PARAM_BC_HMG_NEUMANN];
    bc->neu = cs_cdo_bc_list_create(n_ents, count[CS_PARAM_BC_NEUMANN]);

    n_ents = count[CS_PARAM_BC_ROBIN];
    bc->rob = cs_cdo_bc_list_create(n_ents, n_ents);

    /* Sanity checks */
    cs_lnum_t  n_counted_faces = 0;
    for (i = 0; i < CS_PARAM_N_BC_TYPES; i++)
      n_counted_faces += count[i];
    assert(n_counted_faces == n_b_faces);

    // 4) Define each cs_cdo_bc_list_t
    for (i = 0; i < CS_PARAM_N_BC_TYPES; i++)
      count[i] = 0;

    /* Loop on the definition of each boundary condition */
    for (id = 0; id < param_bc->n_defs; id++) {

      const cs_param_bc_def_t  *def = param_bc->defs + id;
      const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(def->loc_id);
      const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(def->loc_id);

      type = def->bc_type;
      switch (type) {

      case CS_PARAM_BC_DIRICHLET:
        if (elt_ids == NULL) { /* Full selection */
          for (i = 0; i < n_elts[0]; i++) {
            bc->dir->elt_ids[i] = i;
            bc->dir->def_ids[i] = id;
          }
        }
        else { /* Partial selection */

          shift = count[type];
          for (i = 0; i < n_elts[0]; i++) {
            bc->dir->elt_ids[shift + i] = elt_ids[i];
            bc->dir->def_ids[shift + i] = id;
          }
          count[type] += n_elts[0];

        }
        break;

      case CS_PARAM_BC_HMG_DIRICHLET:
        if (elt_ids == NULL) /* Full selection */
          for (i = 0; i < n_elts[0]; i++)
            bc->dir->elt_ids[i] = i;

        else { /* Partial selection */

          shift = count[type] + bc->dir->n_nhmg_elts;
          for (i = 0; i < n_elts[0]; i++)
            bc->dir->elt_ids[shift + i] = elt_ids[i];
          count[type] += n_elts[0];

        }
        break;

      case CS_PARAM_BC_NEUMANN:
        if (elt_ids == NULL) { /* Full selection */
          for (i = 0; i < n_elts[0]; i++) {
            bc->neu->elt_ids[i] = i;
            bc->neu->def_ids[i] = id;
          }
        }
        else { /* Partial selection */

          shift = count[type];
          for (i = 0; i < n_elts[0]; i++) {
            bc->neu->elt_ids[shift + i] = elt_ids[i];
            bc->neu->def_ids[shift + i] = id;
          }
          count[type] += n_elts[0];

        }
        break;

      case CS_PARAM_BC_HMG_NEUMANN:
        if (elt_ids == NULL) /* Full selection */
          for (i = 0; i < n_elts[0]; i++)
            bc->neu->elt_ids[i] = i;

        else { /* Partial selection */

          shift = count[type] + bc->neu->n_nhmg_elts;
          for (i = 0; i < n_elts[0]; i++)
            bc->neu->elt_ids[shift + i] = elt_ids[i];
          count[type] += n_elts[0];

        }
        break;

      case CS_PARAM_BC_ROBIN:
        if (elt_ids == NULL) { /* Full selection */
          for (i = 0; i < n_elts[0]; i++) {
            bc->rob->elt_ids[i] = i;
            bc->rob->def_ids[i] = id;
          }
        }
        else { /* Partial selection */

          shift = count[type];
          for (i = 0; i < n_elts[0]; i++) {
            bc->rob->elt_ids[shift + i] = elt_ids[i];
            bc->rob->def_ids[shift + i] = id;
          }
          count[type] += n_elts[0];

        }
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Invalid type of boundary condition.\n"
                    " Stop generating the boundary condition structure."));
      }

    } // Loop on boundary conditions

    BFT_FREE(bc_types);

  } /* n_b_faces > 0 */

  return bc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cdo_bc_t structure
 *
 * \param[in, out]  face_bc   pointer to a cs_cdo_bc_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_t *
cs_cdo_bc_free(cs_cdo_bc_t   *face_bc)
{
  if (face_bc == NULL)
    return face_bc;

  face_bc->dir = cs_cdo_bc_list_free(face_bc->dir);
  face_bc->neu = cs_cdo_bc_list_free(face_bc->neu);
  face_bc->rob = cs_cdo_bc_list_free(face_bc->rob);

  BFT_FREE(face_bc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build cs_cdo_bc_list_t structures for Dirichlet BC on primal
 *         vertices.
 *         When there is a choice between homogeneous or non homogeneous BCs,
 *         we always set the homogeneous condition for Dirichlet BCs.
 *
 * \param[in] m         pointer to a cs_mesh_t structure
 * \param[in] face_bc   pointer to a cs_cdo_bc_t structure
 *
 * \return a pointer to a new allocated cs_cdo_bc_list_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_bc_list_t *
cs_cdo_bc_vtx_dir_create(const cs_mesh_t    *m,
                         const cs_cdo_bc_t  *face_bc)
{
  cs_lnum_t  i, j, f_id, v_id, def_id;

  cs_lnum_t  n_nhmg_vertices = 0, n_hmg_vertices = 0;
  cs_param_bc_type_t  *vtx_type = NULL;
  short int *vtx_def = NULL;

  const cs_cdo_bc_list_t  *face_dir = face_bc->dir;
  const cs_lnum_t  *f2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t  *f2v_lst = m->b_face_vtx_lst;

  /* Initialization */
  BFT_MALLOC(vtx_type, m->n_vertices, cs_param_bc_type_t);
  BFT_MALLOC(vtx_def, m->n_vertices, short int);
  for (i = 0; i < m->n_vertices; i++) {
    vtx_type[i] = CS_PARAM_N_BC_TYPES;
    vtx_def[i] = -1;
  }

  /* 1) Tag vertices attached to a face where a Dirichlet BC is defined */
  for (i = 0; i < face_dir->n_nhmg_elts; i++) {
    f_id = face_dir->elt_ids[i];
    def_id = face_dir->def_ids[i];
    for (j = f2v_idx[f_id]; j < f2v_idx[f_id+1]; j++) {
      v_id = f2v_lst[j];
      vtx_type[v_id] = CS_PARAM_BC_DIRICHLET;
      vtx_def[v_id] = def_id;
    }
  }

  /* 2) Tag vertices attached to a face where a hmg Dirichlet BC is defined */
  for (i = face_dir->n_nhmg_elts; i < face_dir->n_elts; i++) {
    f_id = face_dir->elt_ids[i];
    for (j = f2v_idx[f_id]; j < f2v_idx[f_id+1]; j++)
      vtx_type[f2v_lst[j]] = CS_PARAM_BC_HMG_DIRICHLET;
  }

  /* 3) In parallel, we need to synchronize the type/def associated to a vertex
     belonging to the parallel frontier (TODO-MPI) */

  /* 4) Count the number of hmg and nhmg Dirichlet vertices */
  for (v_id = 0; v_id < m->n_vertices; v_id++) {

    if (vtx_type[v_id] == CS_PARAM_BC_HMG_DIRICHLET)
      n_hmg_vertices++;
    else if (vtx_type[v_id] == CS_PARAM_BC_DIRICHLET)
      n_nhmg_vertices++;

  }

  /* Define the structure to return */
  cs_lnum_t  n_dir_vertices = n_hmg_vertices + n_nhmg_vertices;
  cs_cdo_bc_list_t  *vtx_dir = cs_cdo_bc_list_create(n_dir_vertices,
                                                     n_nhmg_vertices);

  if (n_dir_vertices > 0) { /* Fill elt_ids */

    /* Non-homogeneous BC are first listed */
    cs_lnum_t  hmg_count = 0, nhmg_count = 0;
    for (v_id = 0; v_id < m->n_vertices; v_id++) {

      if (vtx_type[v_id] == CS_PARAM_BC_DIRICHLET) {
        vtx_dir->elt_ids[nhmg_count] = v_id;
        vtx_dir->def_ids[nhmg_count++] = vtx_def[v_id];
      }
      else if (vtx_type[v_id] == CS_PARAM_BC_HMG_DIRICHLET) {
        vtx_dir->elt_ids[ n_nhmg_vertices + hmg_count] = v_id;
        hmg_count++;
      }

    }

    /* Sanity checks */
    assert(hmg_count == n_hmg_vertices);
    assert(nhmg_count == n_nhmg_vertices);

  } /* n_dir_vertices > 0 */

  BFT_FREE(vtx_type);
  BFT_FREE(vtx_def);

  return vtx_dir;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the Dirichlet values to enforce on the corresponding entities
 *
 * \param[in]      dof_flag  information about the corresponding DoF to treat
 * \param[in]      time_step  pointer to a time step structure
 * \param[in]      geom      structure storing geometric information
 * \param[in]      bc        pointer to a cs_param_bc_t structure
 * \param[in]      ent_dir   pointer to a cs_cdo_bc_list_t
 * \param[in, out] dir_val   array used to store Dirichlet values
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_bc_dirichlet_set(cs_flag_t                dof_flag,
                        const cs_time_step_t    *time_step,
                        const void              *geom,
                        const cs_param_bc_t     *bc,
                        const cs_cdo_bc_list_t  *ent_dir,
                        double                  *dir_val)
{
  cs_real_3_t  xyz;
  cs_get_t  get;

  const double  tcur = time_step->t_cur;

  if (ent_dir->n_nhmg_elts == 0) /* Nothing to compute */
    return;

  /* Sanity check */
  assert(dir_val != NULL);

  for (cs_lnum_t i = 0; i < ent_dir->n_nhmg_elts; i++) {

    cs_lnum_t  id = ent_dir->elt_ids[i];
    cs_param_bc_def_t  *bc_def = bc->defs + ent_dir->def_ids[i];

    if (bc_def->var_type != CS_PARAM_VAR_SCAL)
      bft_error(__FILE__, __LINE__, 0,
                _(" This situation is not handled yet."));

    switch(bc_def->def_type) {

    case CS_PARAM_DEF_BY_VALUE:
      dir_val[i] = bc_def->def_coef1.get.val;
      break;

    case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:

      if (dof_flag & CS_FLAG_VERTEX) { // Values at vertices

        const cs_mesh_t *m = (const cs_mesh_t *)geom;

        for (int k = 0; k < 3; k++) xyz[k] = m->vtx_coord[3*id+k];
        bc_def->def_coef1.analytic(tcur, xyz, &get);
        dir_val[i] = get.val; // stride = 1
        break;

      }
      else if (dof_flag & CS_FLAG_FACE) { // Values at face centers

        const cs_cdo_quantities_t *q = (const cs_cdo_quantities_t *)geom;

        cs_lnum_t  f_id = q->n_i_faces + id;
        for (int k = 0; k < 3; k++) xyz[k] = q->face[f_id].center[k];
        bc_def->def_coef1.analytic(tcur, xyz, &get);
        dir_val[i] = get.val; // stride = 1
        break;

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(" Invalid couple (definition type, degrees of freedom).\n"
                    " Stop computing the Dirichlet values.\n"));
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of definition.\n"
                  " Stop computing the Dirichlet value.\n"));

    } // switch def_type

  } // loop on entities

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
