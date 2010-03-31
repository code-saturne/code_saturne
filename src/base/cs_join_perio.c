/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 2008-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *===========================================================================*/

/*============================================================================
 * Management of conforming and non-conforming joining in case of periodicity
 *===========================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *---------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *---------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *---------------------------------------------------------------------------*/

#include <fvm_parall.h>
#include <fvm_periodicity.h>
#include <fvm_interface.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "cs_search.h"
#include "cs_join_mesh.h"
#include "cs_join_post.h"
#include "cs_join_set.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_perio.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Structure Definitions
 *===========================================================================*/

typedef struct _cs_join_perio_builder_t {

  int   n_perio;   /* Number of initial periodicity
                      No composition of periodicity taken into account */

  fvm_periodicity_t   *periodicity;   /* Structure keeping periodicity
                                         information like matrix of
                                         the transformation */

  int   *n_perio_couples;             /* Local number of periodic face
                                         couples for each periodicity */
  fvm_gnum_t  **perio_couples;        /* List of global numbering of
                                         periodic faces. */

} cs_join_perio_builder_t;

/*============================================================================
 * Global variables
 *===========================================================================*/

int  cs_glob_n_join_perio = 0; /* Number of periodicity defined through
                                  a joining operation */

/*============================================================================
 * Static global variables
 *===========================================================================*/

static cs_join_perio_builder_t  *cs_glob_join_perio_builder = NULL;

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Create and initialize a cs_join_perio_builder_t structure.
 *
 * returns:
 *  a pointer to a new allocated cs_join_perio_builder_t
 *---------------------------------------------------------------------------*/

static cs_join_perio_builder_t *
_create_perio_builder(void)
{
  cs_join_perio_builder_t  *perio_builder = NULL;

  BFT_MALLOC(perio_builder, 1, cs_join_perio_builder_t);

  perio_builder->n_perio = 1;

  perio_builder->periodicity = fvm_periodicity_create(0.001);

  BFT_MALLOC(perio_builder->n_perio_couples, 1, int);
  BFT_MALLOC(perio_builder->perio_couples, 1, fvm_gnum_t *);

  perio_builder->n_perio_couples[0] = 0;
  perio_builder->perio_couples[0] = NULL;

  return perio_builder;
}

/*----------------------------------------------------------------------------
 * Free a cs_join_perio_builder_t structure.
 * builder->periodicity is not freed here because it's transfered to a
 * cs_mesh_t structure.
 *
 * parameter:
 *   builder     <-> a pointer to a cs_join_perio_builder_t struct. to free
 *---------------------------------------------------------------------------*/

static void
_delete_perio_builder(cs_join_perio_builder_t   **builder)
{
  int  i;

  cs_join_perio_builder_t  *_bdr = *builder;

  if (_bdr == NULL)
    return;

  for (i = 0; i < _bdr->n_perio; i++)
    BFT_FREE(_bdr->perio_couples[i]);
  BFT_FREE(_bdr->perio_couples);
  BFT_FREE(_bdr->n_perio_couples);

  BFT_FREE(_bdr);

  *builder = NULL;
}

/*----------------------------------------------------------------------------
 * Update a cs_join_perio_builder_t structure. Add a new periodicity.
 *
 * parameter:
 *   builder     <-> a pointer to a cs_join_perio_builder_t structure
 *---------------------------------------------------------------------------*/

static void
_increment_perio_builder(cs_join_perio_builder_t   *builder)
{
  assert(builder != NULL);

  builder->n_perio += 1;

  BFT_REALLOC(builder->n_perio_couples, builder->n_perio, int);
  BFT_REALLOC(builder->perio_couples, builder->n_perio, fvm_gnum_t *);

  builder->n_perio_couples[builder->n_perio - 1] = 0;
  builder->perio_couples[builder->n_perio - 1] = NULL;
}

/*----------------------------------------------------------------------------
 * Add a new periodic cs_join_t structure.
 *
 * parameters:
 *   join_number  <-- number related to the joining operation
 *   sel_criteria <-- boundary face selection criteria
 *   fraction     <-- value of the fraction parameter
 *   plane        <-- value of the plane parameter
 *   perio_num    <-- periodicity number
 *   verbosity    <-- level of verbosity required
 *
 * returns:
 *   a pointer to a new allocated cs_join_t structure
 *---------------------------------------------------------------------------*/

static void
_add_perio_join(int           join_number,
                const char   *criteria,
                float         fraction,
                float         plane,
                int           perio_num,
                int           verbosity)
{
  /* Check parameters value */

  if (fraction < 0.0 || fraction >= 1.0)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the fraction parameter.\n"
                "  It must be between [0.0, 1.0[ and is here: %f\n"),
              fraction);

  if (plane < 0.0 || plane >= 90.0)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the plane parameter.\n"
                "  It must be between [0, 90] and is here: %f\n"),
              plane);

  if (perio_num < 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the periodicity number.\n"
                "  It must be between > 0 and is here: %d\n"),
              perio_num);

   /* Allocate and initialize a cs_join_t structure */

  BFT_REALLOC(cs_glob_join_array, cs_glob_n_joinings + 1, cs_join_t *);

  cs_glob_join_array[cs_glob_n_joinings] = cs_join_create(join_number,
                                                          criteria,
                                                          fraction,
                                                          plane,
                                                          perio_num,
                                                          verbosity);

  cs_glob_n_joinings++;
}

/*----------------------------------------------------------------------------
 * Update work_mesh by redistributing local join mesh.
 *
 * parameters:
 *   param             <--  set of user-defined parameter
 *   gnum_rank_index   <--  index on ranks for the old global face numbering
 *   local_mesh        <--  mesh on local selected faces to be joined
 *   p_work_mesh       <->  distributed mesh on faces to join
 *---------------------------------------------------------------------------*/

static void
_redistribute_mesh(cs_join_param_t         param,
                   const fvm_gnum_t        gnum_rank_index[],
                   const cs_join_mesh_t   *local_mesh,
                   cs_join_mesh_t        **p_work_mesh)
{
  cs_join_mesh_t  *work_mesh = *p_work_mesh;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  /* sanity checks */

  assert(local_mesh != NULL);
  assert(work_mesh != NULL);

  if (n_ranks == 1)
    cs_join_mesh_copy(&work_mesh, local_mesh);

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Parallel mode */

    int  i, n_work_faces;

    char  *mesh_name = NULL;
    fvm_gnum_t  *work_faces = NULL;

    n_work_faces = work_mesh->n_faces;
    BFT_MALLOC(work_faces, n_work_faces, fvm_gnum_t);

    for (i = 0; i < n_work_faces; i++)
      work_faces[i] = work_mesh->face_gnum[i];

    /* Replace current work mesh */

    cs_join_mesh_destroy(&work_mesh);

    BFT_MALLOC(mesh_name, strlen("WorkMesh_j_n") + 2 + 5 + 1, char);
    sprintf(mesh_name,"%s%02d%s%05d",
            "WorkMesh_j", param.num, "_n", local_rank);

    work_mesh = cs_join_mesh_create_from_glob_sel(mesh_name,
                                                  n_work_faces,
                                                  work_faces,
                                                  gnum_rank_index,
                                                  local_mesh);

    BFT_FREE(mesh_name);
    BFT_FREE(work_faces);

  }
#endif

  /* Return pointers */

  *p_work_mesh = work_mesh;

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Complete mesh builder definition for periodic faces.
 * Algorithm similar to the one in the file cs_preprocessor_data.c.
 * Face relation are defined in local numbering.
 *
 * parameters:
 *   n_init_perio <-- number of initial periodicities defined
 *   mb           <-> mesh builder structure to update
 *   face_ifs     <-- face interface set structure
 *---------------------------------------------------------------------------*/

static void
_extract_perio_couples(int                   n_init_perio,
                       cs_mesh_builder_t    *mb,
                       fvm_interface_set_t  *face_ifs)
{
  int  i, j, k, l;

  fvm_lnum_t  *per_face_count = NULL;
  fvm_lnum_t  *if_index = NULL;
  fvm_lnum_t  *send_num = NULL, *recv_num = NULL;

  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const int  n_interfaces = fvm_interface_set_size(face_ifs);
  const fvm_lnum_t  tr_index_size = n_init_perio*2 + 2;

  /* Define interface index only for the periodic data */

  BFT_MALLOC(if_index, n_interfaces + 1, fvm_lnum_t);

  if_index[0] = 0;
  for (j = 0; j < n_interfaces; j++) {
    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);
    if_index[j+1] = if_index[j] + tr_index[tr_index_size - 1]-tr_index[1];
  }

  /* Define send_num and recv_num */

  BFT_MALLOC(send_num, if_index[n_interfaces], fvm_lnum_t);
  BFT_MALLOC(recv_num, if_index[n_interfaces], fvm_lnum_t);

  for (j = 0; j < n_interfaces; j++) {

    const fvm_lnum_t start_id = if_index[j];
    const fvm_lnum_t end_id = if_index[j+1];
    const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
    const fvm_lnum_t *tr_index = fvm_interface_get_tr_index(face_if);
    const fvm_lnum_t *loc_num = fvm_interface_get_local_num(face_if);
    const int distant_rank = fvm_interface_rank(face_if);

    for (k = start_id, l = tr_index[1]; k < end_id; k++, l++)
      send_num[k] = loc_num[l];

    if (distant_rank == local_rank) {
      const fvm_lnum_t *dist_num = fvm_interface_get_distant_num(face_if);
      for (k = start_id, l = tr_index[1]; k < end_id; k++, l++)
        recv_num[k] = dist_num[l];
    }
  }

  {   /* Exchange local face numbers */

    MPI_Request  *request = NULL;
    MPI_Status  *status  = NULL;

    int request_count = 0;

    BFT_MALLOC(request, n_interfaces*2, MPI_Request);
    BFT_MALLOC(status, n_interfaces*2, MPI_Status);

    for (j = 0; j < n_interfaces; j++) {
      const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
      int distant_rank = fvm_interface_rank(face_if);

      if (distant_rank != local_rank)
        MPI_Irecv(recv_num + if_index[j],
                  if_index[j+1] - if_index[j],
                  FVM_MPI_LNUM,
                  distant_rank,
                  distant_rank,
                  cs_glob_mpi_comm,
                  &(request[request_count++]));
    }

    for (j = 0; j < n_interfaces; j++) {
      const fvm_interface_t *face_if = fvm_interface_set_get(face_ifs, j);
      int distant_rank = fvm_interface_rank(face_if);

      if (distant_rank != local_rank)
        MPI_Isend(send_num + if_index[j],
                  if_index[j+1] - if_index[j],
                  FVM_MPI_LNUM,
                  distant_rank,
                  (int)local_rank,
                  cs_glob_mpi_comm,
                  &(request[request_count++]));
    }

    MPI_Waitall(request_count, request, status);

    BFT_FREE(request);
    BFT_FREE(status);
  }

  /* Copy interface information to mesh builder */

  BFT_MALLOC(per_face_count, n_init_perio, fvm_lnum_t);
  for (i = 0; i < n_init_perio; i++)
    per_face_count[i] = 0;

  for (j = 0; j < n_interfaces; j++) {

    fvm_lnum_t  tr_shift = 0;
    const fvm_interface_t  *face_if = fvm_interface_set_get(face_ifs, j);
    const int  distant_rank = fvm_interface_rank(face_if);
    const fvm_lnum_t  *tr_index = fvm_interface_get_tr_index(face_if);

    for (i = 1; i < tr_index_size - 1; i++) {

      fvm_lnum_t n_elts = tr_index[i+1] - tr_index[i];

      if ((distant_rank != local_rank) || (i%2 == 1)) {

        fvm_lnum_t  send_shift, recv_shift;

        int  perio_id = (i-1)/2;
        int  perio_sgn = (i%2)*2 - 1; /* 1 for odd, -1 for even */
        fvm_lnum_t  n_dir_elts = tr_index[2*perio_id+2]-tr_index[2*perio_id+1];
        fvm_lnum_t  n_rev_elts = tr_index[2*perio_id+3]-tr_index[2*perio_id+2];

        send_shift = if_index[j] + tr_shift;
        if (distant_rank != local_rank) {
          if (perio_sgn > 0)
            recv_shift = if_index[j] + n_rev_elts + tr_shift;
          else
            recv_shift = if_index[j] - n_dir_elts + tr_shift;
        }
        else /* if (i%2 == 1) */
          recv_shift = send_shift;

        for (k = 0; k < n_elts; k++) {
          l = mb->per_face_idx[perio_id] + per_face_count[perio_id];
          mb->per_face_lst[l*2]     = send_num[send_shift + k]*perio_sgn;
          mb->per_face_lst[l*2 + 1] = recv_num[recv_shift + k];
          mb->per_rank_lst[l] = distant_rank + 1;
          per_face_count[perio_id] += 1;
        }
      }

      tr_shift += n_elts;

    } /* End of loop on tr_index */

  } /* End of loop on interfaces */


  BFT_FREE(per_face_count);
  BFT_FREE(recv_num);
  BFT_FREE(send_num);
  BFT_FREE(if_index);

}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Delete interior faces temporary added for periodicity operation.
 * Only done if n_ranks > 1 because in serial, these periodic faces are
 * merged with initially present faces
 *
 * parameters:
 *   param <-- set of parameters for the joining operation
 *   mesh  <-> cs_mesh_t structure to update
 *---------------------------------------------------------------------------*/

static void
_perio_face_clean(cs_join_param_t      param,
                  cs_mesh_t           *mesh)
{
  int  i, j, k, shift;

  int  n_ii_faces = mesh->n_i_faces;
  int  n_fi_faces = 0;
  cs_int_t  *new_f2v_idx = NULL;
  int  *tag = NULL;

  assert(cs_glob_n_ranks > 1);

  BFT_MALLOC(tag, n_ii_faces, int);

  for (i = 0; i < n_ii_faces; i++) {

    if (mesh->i_face_cells[2*i] == 0 && mesh->i_face_cells[2*i+1] == 0)
      tag[i] = -1;
    else {
      mesh->i_face_cells[2*n_fi_faces] = mesh->i_face_cells[2*i];
      mesh->i_face_cells[2*n_fi_faces+1] = mesh->i_face_cells[2*i+1];
      n_fi_faces++;
      tag[i] = n_fi_faces;
    }

  }

  if (param.verbosity > 2)
    bft_printf(_("\n  Delete %d interior periodic faces locally\n"),
               n_ii_faces - n_fi_faces);

  mesh->n_i_faces = n_fi_faces;
  BFT_REALLOC(mesh->i_face_cells, 2*mesh->n_i_faces, cs_int_t);
  BFT_MALLOC(new_f2v_idx, n_fi_faces + 1, cs_int_t);

  n_fi_faces = 0;
  for (i = 0; i < n_ii_faces; i++) {
    if (tag[i] > 0) {
      mesh->global_i_face_num[n_fi_faces] = mesh->global_i_face_num[i];
      mesh->i_face_family[n_fi_faces] = mesh->i_face_family[i];
      new_f2v_idx[n_fi_faces + 1] =  mesh->i_face_vtx_idx[i+1]
                                   - mesh->i_face_vtx_idx[i];
      n_fi_faces++;
    }
  }

  BFT_REALLOC(mesh->global_i_face_num, mesh->n_i_faces, fvm_gnum_t);
  BFT_REALLOC(mesh->i_face_family, mesh->n_i_faces, cs_int_t);

  /* Update interior face connectivity */

  new_f2v_idx[0] = 1;
  for (i = 0; i < n_fi_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  n_fi_faces = 0;
  for (i = 0; i < n_ii_faces; i++) {
    if (tag[i] > 0) {
      shift = new_f2v_idx[n_fi_faces] - 1;
      for (k = 0, j = mesh->i_face_vtx_idx[i]-1;
           j < mesh->i_face_vtx_idx[i+1]-1; j++, k++)
        mesh->i_face_vtx_lst[shift+k] = mesh->i_face_vtx_lst[j];
      n_fi_faces++;
    }
  }

  BFT_REALLOC(mesh->i_face_vtx_lst, new_f2v_idx[n_fi_faces]-1, cs_int_t);
  BFT_FREE(mesh->i_face_vtx_idx);

  mesh->i_face_vtx_idx = new_f2v_idx;

  /* There is no need to define a new glbal interior face numbering
     because the excluded faces are always defined on an another rank */

  BFT_FREE(tag);
}

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get the number of periodic transformations already defined
 *
 * Fortran Interface:
 *
 * SUBROUTINE NUMPER
 * *****************
 *
 * INTEGER   numper      : --> : number of periodicities  op. already defined
 *----------------------------------------------------------------------------*/

void CS_PROCF(numper, NUMPER)
(
 cs_int_t    *numper
)
{
  *numper = cs_glob_n_join_perio;

  return;
}

/*----------------------------------------------------------------------------
 * Define a translation
 *
 * Fortran Interface:
 *
 * SUBROUTINE DEFPT1
 * *****************
 *
 * INTEGER        numper    : <-- : number related to the periodic op.
 * CHARACTER*     criteria  : <-- : boundary face selection criteria
 * REAL           fraction  : <-- : parameter for merging vertices
 * REAL           plane     : <-- : parameter for splitting faces
 * INTEGER        verbosity : <-- : verbosity level
 * REAL           tx        : <-- : X coordinate of the translation vector
 * REAL           ty        : <-- : Y coordinate of the translation vector
 * REAL           tz        : <-- : Z coordinate of the translation vector
 * INTEGER        crit_len  : <-- : length of criteria
 *----------------------------------------------------------------------------*/

void CS_PROCF(defpt1, DEFPT1)
(
 cs_int_t    *numper,
 const char  *criteria,
 cs_real_t   *fraction,
 cs_real_t   *plane,
 cs_int_t    *verbosity,
 cs_real_t   *tx,
 cs_real_t   *ty,
 cs_real_t   *tz,
 cs_int_t    *crit_len
 CS_ARGF_SUPP_CHAINE
)
{
  double  trans[3];

  char *_criteria = NULL;

  assert(CS_ABS(*tx) > 0.0 || CS_ABS(*ty) > 0.0 || CS_ABS(*tz) > 0.0);

  trans[0] = *tx;
  trans[1] = *ty;
  trans[2] = *tz;

  if (criteria != NULL && *crit_len > 0)
    _criteria = cs_base_string_f_to_c_create(criteria, *crit_len);
  if (_criteria != NULL && strlen(_criteria) == 0)
    cs_base_string_f_to_c_free(&_criteria);

  bft_printf(_("  Adding periodicity %d "
               "(translation [%10.4e, %10.4e, %10.4e]).\n"),
             cs_glob_n_join_perio, trans[0], trans[1], trans[2]);

  cs_join_perio_add_translation(*numper,
                                _criteria,
                                *fraction,
                                *plane,
                                *verbosity,
                                trans);

  if (_criteria != NULL)
    cs_base_string_f_to_c_free(&_criteria);
}

/*----------------------------------------------------------------------------
 * Define a rotation
 *
 * Fortran Interface:
 *
 * SUBROUTINE DEFPR1
 * *****************
 *
 * INTEGER        numper    : <-- : number related to the periodic op.
 * CHARACTER*     criteria  : <-- : boundary face selection criteria
 * REAL           fraction  : <-- : parameter for merging vertices
 * REAL           plane     : <-- : parameter for splitting faces
 * INTEGER        verbosity : <-- : verbosity level
 * REAL           ax        : <-- : X coordinate of the rotation axis
 * REAL           ay        : <-- : Y coordinate of the rotation axis
 * REAL           az        : <-- : Z coordinate of the rotation axis
 * REAL           theta     : <-- : angle of the rotation (radian)
 * REAL           ix        : <-- : X coordinate of the invariant point
 * REAL           iy        : <-- : Y coordinate of the invariant point
 * REAL           iz        : <-- : Z coordinate of the invariant point
 * INTEGER        crit_len  : <-- : length of criteria string
 *----------------------------------------------------------------------------*/

void CS_PROCF(defpr1, DEFPR1)
(
 cs_int_t    *numper,
 const char  *criteria,
 cs_real_t   *fraction,
 cs_real_t   *plane,
 cs_int_t    *verbosity,
 cs_real_t   *ax,
 cs_real_t   *ay,
 cs_real_t   *az,
 cs_real_t   *theta,
 cs_real_t   *ix,
 cs_real_t   *iy,
 cs_real_t   *iz,
 cs_int_t    *crit_len
 CS_ARGF_SUPP_CHAINE
)
{
  double  axis[3], inv[3];
  char  *_criteria = NULL;

  assert(CS_ABS(*ax) > 0.0 || CS_ABS(*ay) > 0.0 || CS_ABS(*az) > 0.0);
  assert(CS_ABS(*theta) > 0.0);

  axis[0] = *ax;
  axis[1] = *ay;
  axis[2] = *az;

  inv[0] = *ix;
  inv[1] = *iy;
  inv[2] = *iz;

  if (criteria != NULL && *crit_len > 0)
    _criteria = cs_base_string_f_to_c_create(criteria, *crit_len);
  if (_criteria != NULL && strlen(_criteria) == 0)
    cs_base_string_f_to_c_free(&_criteria);

  bft_printf(_("  Adding periodicity %d (rotation).\n"), cs_glob_n_join_perio);

  cs_join_perio_add_rotation(*numper,
                             _criteria,
                             *fraction,
                             *plane,
                             *verbosity,
                             *theta,
                             axis,
                             inv);

  if (_criteria != NULL)
    cs_base_string_f_to_c_free(&_criteria);
}

/*----------------------------------------------------------------------------
 * Define a general transformation through a homogeneous matrix (4x4)
 *     _               _
 *    | r11 r12 r13 tx  |  t(x,y,z) : translation vector
 *    | r21 r22 r23 ty  |  r(i,j)   : rotation matrix
 *    | r31 r32 r33 tz  |
 *    |_  0   0   0  1 _|
 *
 * Fortran Interface:
 *
 * SUBROUTINE DEFPG1
 * *****************
 *
 * INTEGER        numper    : <-- : number related to the periodic op.
 * CHARACTER*     criteria  : <-- : boundary face selection criteria
 * REAL           fraction  : <-- : parameter for merging vertices
 * REAL           plane     : <-- : parameter for splitting faces
 * INTEGER        verbosity : <-- : verbosity level
 * REAL           r11       : <-- : coef. (1,1) of the homogeneous matrix
 * REAL           r12       : <-- : coef. (1,2) of the homogeneous matrix
 * REAL           r13       : <-- : coef. (1,3) of the homogeneous matrix
 * REAL           tx        : <-- : coef. (1,4) of the homogeneous matrix
 * REAL           r21       : <-- : coef. (2,1) of the homogeneous matrix
 * REAL           r22       : <-- : coef. (2,2) of the homogeneous matrix
 * REAL           r23       : <-- : coef. (2,3) of the homogeneous matrix
 * REAL           ty        : <-- : coef. (2,4) of the homogeneous matrix
 * REAL           r31       : <-- : coef. (3,1) of the homogeneous matrix
 * REAL           r32       : <-- : coef. (3,2) of the homogeneous matrix
 * REAL           r33       : <-- : coef. (3,3) of the homogeneous matrix
 * REAL           tz        : <-- : coef. (3,4) of the homogeneous matrix
 * INTEGER        crit_len  : <-- : length of criteria string
 *----------------------------------------------------------------------------*/

void CS_PROCF(defpg1, DEFPG1)
(
 cs_int_t    *numper,
 const char  *criteria,
 cs_real_t   *fraction,
 cs_real_t   *plane,
 cs_int_t    *verbosity,
 cs_real_t   *r11,
 cs_real_t   *r12,
 cs_real_t   *r13,
 cs_real_t   *tx,
 cs_real_t   *r21,
 cs_real_t   *r22,
 cs_real_t   *r23,
 cs_real_t   *ty,
 cs_real_t   *r31,
 cs_real_t   *r32,
 cs_real_t   *r33,
 cs_real_t   *tz,
 cs_int_t    *crit_len
 CS_ARGF_SUPP_CHAINE
)
{
  double  matrix[3][4];
  char  *_criteria = NULL;

  /* Build the matrix */

  matrix[0][0] = *r11;
  matrix[0][1] = *r12;
  matrix[0][2] = *r13;
  matrix[0][3] = *tx;

  matrix[1][0] = *r21;
  matrix[1][1] = *r22;
  matrix[1][2] = *r23;
  matrix[1][3] = *ty;

  matrix[2][0] = *r31;
  matrix[2][1] = *r32;
  matrix[2][2] = *r33;
  matrix[2][3] = *tz;

  if (criteria != NULL && *crit_len > 0)
    _criteria = cs_base_string_f_to_c_create(criteria, *crit_len);
  if (_criteria != NULL && strlen(_criteria) == 0)
    cs_base_string_f_to_c_free(&_criteria);

  bft_printf(_("  Adding periodicity %d (general formulation).\n"),
             cs_glob_n_join_perio);

  cs_join_perio_add_mixed(*numper,
                          _criteria,
                          *fraction,
                          *plane,
                          *verbosity,
                          (double (*)[4])matrix);

  if (_criteria != NULL)
    cs_base_string_f_to_c_free(&_criteria);
}

/*----------------------------------------------------------------------------
 * Set advanced parameters for the joining algorithm in case of periodicity
 *
 * Fortran Interface:
 *
 * SUBROUTINE SETAPP
 * *****************
 *
 * INTEGER      perio_num         : <-- : perio number
 * REAL         mtf               : <-- : merge tolerance coefficient
 * REAL         pmf               : <-- : pre-merge factor
 * INTEGER      tcm               : <-- : tolerance computation mode
 * INTEGER      icm               : <-- : intersection computation mode
 * INTEGER      maxbrk            : <-- : max number of equiv. breaks
 * INTEGER      max_sub_faces     : <-- : max. possible number of sub-faces
 *                                        by splitting a selected face
 * INTEGER      tml               : <-- : tree max level
 * INTEGER      tmb               : <-- : tree max boxes
 * REAL         tmr               : <-- : tree max ratio
 *----------------------------------------------------------------------------*/

void CS_PROCF(setapp, SETAPP)
(
 cs_int_t    *perio_num,
 cs_real_t   *mtf,
 cs_real_t   *pmf,
 cs_int_t    *tcm,
 cs_int_t    *icm,
 cs_int_t    *maxbrk,
 cs_int_t    *max_sub_faces,
 cs_int_t    *tml,
 cs_int_t    *tmb,
 cs_real_t   *tmr
)
{
  int  i, join_id = -1;

  cs_join_t  *join = NULL;

  assert(*perio_num > 0);

  /* Look for the joining structure related to "perio_num" */

  for (i = 0; i < cs_glob_n_joinings; i++) {

    join = cs_glob_join_array[i];
    if (*perio_num == join->param.perio_num) {
      join_id = i;
      break;
    }

  }

  if (join_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" Periodicity number %d is not defined\n"
                " %d periodicities are defined\n"),
              *perio_num, cs_glob_n_join_perio);

  assert(join != NULL);

  cs_join_set_advanced_param(join,
                             *mtf,
                             *pmf,
                             *tcm,
                             *icm,
                             *maxbrk,
                             *max_sub_faces,
                             *tml,
                             *tmb,
                             *tmr);

}

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Define a translational periodicity
 *
 * parameters:
 *   perio_num    <-- number related to the periodicity
 *   sel_criteria <-- boundary face selection criteria
 *   fraction     <-- value of the fraction parameter
 *   plane        <-- value of the plane parameter
 *   verbosity    <-- level of verbosity required
 *   trans        <-- translation vector
 *----------------------------------------------------------------------------*/

void
cs_join_perio_add_translation(int            perio_num,
                              const char    *sel_criteria,
                              double         fraction,
                              double         plane,
                              int            verbosity,
                              const double   trans[3])
{
  cs_int_t  tr_id;

  assert((trans[0]*trans[0] + trans[1]*trans[1] + trans[2]*trans[2]) > 0.0);

  if (cs_glob_n_join_perio == 0) {
    assert(cs_glob_join_perio_builder == NULL);
    cs_glob_join_perio_builder = _create_perio_builder();
  }
  else
    _increment_perio_builder(cs_glob_join_perio_builder);

  cs_glob_n_join_perio++;

  tr_id =
    fvm_periodicity_add_translation(cs_glob_join_perio_builder->periodicity,
                                    cs_glob_n_join_perio,
                                    trans);

  _add_perio_join(perio_num + cs_glob_join_count,
                  sel_criteria,
                  fraction,
                  plane,
                  cs_glob_n_join_perio,
                  verbosity);
}

/*----------------------------------------------------------------------------
 * Define a rotational periodicity
 *
 * parameters:
 *   perio_num    <-- number related to the periodicity
 *   sel_criteria <-- boundary face selection criteria
 *   fraction     <-- value of the fraction parameter
 *   plane        <-- value of the plane parameter
 *   verbosity    <-- level of verbosity required
 *   theta        <-- rotation angle (in degrees)
 *   axis         <-- axis vector
 *   invariant    <-- invariant point coordinates
 *----------------------------------------------------------------------------*/

void
cs_join_perio_add_rotation(int            perio_num,
                           const char    *sel_criteria,
                           double         fraction,
                           double         plane,
                           int            verbosity,
                           double         theta,
                           const double   axis[3],
                           const double   invariant[3])
{
  cs_int_t  tr_id;

  assert(theta*theta > 0.0);
  assert((axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]) > 0.0);

  if (cs_glob_n_join_perio == 0) {
    assert(cs_glob_join_perio_builder == NULL);
    cs_glob_join_perio_builder = _create_perio_builder();
  }
  else
    _increment_perio_builder(cs_glob_join_perio_builder);

  cs_glob_n_join_perio++;

  tr_id =
    fvm_periodicity_add_rotation(cs_glob_join_perio_builder->periodicity,
                                 cs_glob_n_join_perio,
                                 theta,
                                 axis,
                                 invariant);

  _add_perio_join(perio_num + cs_glob_join_count,
                  sel_criteria,
                  fraction,
                  plane,
                  cs_glob_n_join_perio,
                  verbosity);
}

/*----------------------------------------------------------------------------
 * Define a periodicity using a matrix
 *
 * parameters:
 *   perio_num    <-- number related to the periodicity
 *   sel_criteria <-- boundary face selection criteria
 *   fraction     <-- value of the fraction parameter
 *   plane        <-- value of the plane parameter
 *   verbosity    <-- level of verbosity required
 *   matrix       <-- transformation matrix
 *----------------------------------------------------------------------------*/

void
cs_join_perio_add_mixed(int            perio_num,
                        const char    *sel_criteria,
                        double         fraction,
                        double         plane,
                        int            verbosity,
                        double         matrix[3][4])
{
  cs_int_t  tr_id;

  if (cs_glob_n_join_perio == 0) {
    assert(cs_glob_join_perio_builder == NULL);
    cs_glob_join_perio_builder = _create_perio_builder();
  }
  else
    _increment_perio_builder(cs_glob_join_perio_builder);

  cs_glob_n_join_perio++;

  tr_id =
    fvm_periodicity_add_by_matrix(cs_glob_join_perio_builder->periodicity,
                                  cs_glob_n_join_perio,
                                  FVM_PERIODICITY_MIXED,
                                  matrix);

  _add_perio_join(perio_num + cs_glob_join_count,
                  sel_criteria,
                  fraction,
                  plane,
                  cs_glob_n_join_perio,
                  verbosity);
}

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Modify compact_face_gnum to take into account
 * new periodic faces and create a periodic vertex couple list.
 *
 * parameters:
 *   this_join <-- high level join structure
 *   jmesh     <-> local join mesh struct. to duplicate and transform
 *   mesh      <-- pointer to a cs_mesh_t struct.
 *---------------------------------------------------------------------------*/

void
cs_join_perio_apply(cs_join_t          *this_join,
                    cs_join_mesh_t     *jmesh,
                    const cs_mesh_t    *mesh)
{
  cs_int_t  i, j, k, shift;
  cs_real_t  matrix[3][4], xyz[4];

  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *select = this_join->selection;
  fvm_periodicity_t  *periodicity = cs_glob_join_perio_builder->periodicity;

  const int  n_ranks = cs_glob_n_ranks;
  const int  perio_id = param.perio_num - 1;
  const int  n_init_vertices = jmesh->n_vertices;
  const int  n_init_faces = jmesh->n_faces;

  assert(perio_id > -1);

  /* Retrieve related transformation */

  fvm_periodicity_get_matrix(periodicity, 2*perio_id, matrix);

  /* Duplicate and transform vertices */

  jmesh->n_vertices *= 2;
  jmesh->n_g_vertices *= 2;

  BFT_REALLOC(jmesh->vertices, jmesh->n_vertices, cs_join_vertex_t);

  shift = n_init_vertices;
  for (i = 0; i < n_init_vertices; i++) {

    /* Copy tolerance, coord, state and gnum to initialize new_vtx */

    cs_join_vertex_t  new_vtx = jmesh->vertices[i];

    for (j = 0; j < 3; j++) {
      xyz[j] = new_vtx.coord[j];
      new_vtx.coord[j] = 0.0;
    }
    xyz[3] = 1;

    for (j = 0; j < 3; j++)
      for (k = 0; k < 4; k++)
        new_vtx.coord[j] += matrix[j][k]*xyz[k];

    new_vtx.state = CS_JOIN_STATE_PERIO;
    jmesh->vertices[shift++] = new_vtx;

  }

  /* Add a periodic vertex couple list */

  select->n_couples = n_init_vertices;
  BFT_MALLOC(select->per_v_couples, 2*n_init_vertices, fvm_gnum_t);

  if (n_ranks > 1) { /* Global numbering update */

    fvm_gnum_t  *gnum = NULL;
    fvm_io_num_t  *io_num = NULL;
    const fvm_gnum_t  *io_gnum = NULL;

    BFT_MALLOC(gnum, n_init_vertices, fvm_gnum_t);

    for (i = 0, shift = n_init_vertices; i < n_init_vertices; i++, shift++)
      gnum[i] = jmesh->vertices[shift].gnum;

    io_num = fvm_io_num_create(NULL, gnum, n_init_vertices, 0);
    io_gnum = fvm_io_num_get_global_num(io_num);

    for (i = 0, shift = n_init_vertices; i < n_init_vertices; i++, shift++) {
      jmesh->vertices[shift].gnum = io_gnum[i] + mesh->n_g_vertices;
      select->per_v_couples[2*i] = jmesh->vertices[i].gnum;
      select->per_v_couples[2*i+1] = jmesh->vertices[shift].gnum;
    }

    fvm_io_num_destroy(io_num);
    BFT_FREE(gnum);

  }
  else { /* Serial run */

    for (i = 0, shift = n_init_vertices; i < n_init_vertices; i++, shift++) {
      jmesh->vertices[shift].gnum = i + 1 + mesh->n_g_vertices;
      select->per_v_couples[2*i] = jmesh->vertices[i].gnum;
      select->per_v_couples[2*i+1] = jmesh->vertices[shift].gnum;
    }

  }

  /* Duplicate and transform faces */

  jmesh->n_faces *= 2;
  jmesh->n_g_faces *= 2;

  BFT_REALLOC(jmesh->face_vtx_idx, jmesh->n_faces + 1, cs_int_t);
  BFT_REALLOC(jmesh->face_gnum, jmesh->n_faces, fvm_gnum_t);
  BFT_REALLOC(jmesh->face_vtx_lst,
              2*(jmesh->face_vtx_idx[n_init_faces]-1), cs_int_t);

  for (i = 0; i < n_init_faces; i++) {

    int  pfid = n_init_faces + i;
    int  s = jmesh->face_vtx_idx[i] - 1;
    int  e = jmesh->face_vtx_idx[i+1] - 1;
    int  ps = jmesh->face_vtx_idx[pfid] - 1;
    int  pe = jmesh->face_vtx_idx[pfid] - 1 + e - s;
    fvm_gnum_t  new_gnum = 2*jmesh->face_gnum[i];

    jmesh->face_gnum[i] = new_gnum - 1;
    jmesh->face_gnum[pfid] = new_gnum;

    for (j = s, shift = ps; j < e; j++, shift++)
      jmesh->face_vtx_lst[shift] = jmesh->face_vtx_lst[j] + n_init_vertices;
    jmesh->face_vtx_idx[pfid+1] =  pe + 1;

  }

  /* Modify cs_join_select_t structure */

  for (i = 0; i < n_ranks + 1; i++)
    select->compact_rank_index[i] *= 2;

  for (i = 0; i < select->n_faces; i++)
    select->compact_face_gnum[i] = 2*select->compact_face_gnum[i] - 1;

  /* Order faces in mesh struct. by increasing global face number.
     We have to be order in this way in order to keep an exact
     redistribution (work_mesh) */

  cs_join_mesh_face_order(jmesh);

  if (param.verbosity > 1)
    bft_printf(_("  Apply periodicity to the local join mesh structure\n"
                 "  New number of faces to treat locally: %8d\n"),
               jmesh->n_faces);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n Periodic vertex couples:\n");
  for (i = 0; i < n_init_vertices; i++)
    bft_printf(" %6d  (%9u, %9u)\n", i+1,
               select->per_v_couples[2*i],  select->per_v_couples[2*i+1]);
  bft_printf_flush();
#endif

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(_("\n  Selected faces for the joining operation:\n"));
  for (i = 0; i < select->n_faces; i++)
    bft_printf(" %9d | %9d | %10u | %10u\n",
               i, select->faces[i], select->compact_face_gnum[i],
               select->cell_gnum[i]);
  bft_printf("\n");
#endif

}

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Modify compact_face_gnum to take into account
 * new periodic faces and create a periodic vertex couple list.
 *
 * parameters:
 *   this_join          <-- pointer to a high level join structure
 *   jmesh              <-> local join mesh struct. to duplicate and transform
 *   p_work_jmesh       <-> distributed join mesh struct. on which operations
 *                          take place
 *   p_work_edges       <-> join edges struct. related to work_jmesh
 *   init_max_vtx_gnum  <-- initial max. global numbering for vertices
 *   n_g_new_vertices   <-- global number of vertices created during the
 *                          intersection of edges
 *---------------------------------------------------------------------------*/

void
cs_join_perio_merge_back(cs_join_t          *this_join,
                         cs_join_mesh_t     *jmesh,
                         cs_join_mesh_t    **p_work_jmesh,
                         cs_join_edges_t   **p_work_edges,
                         fvm_gnum_t          init_max_vtx_gnum,
                         fvm_gnum_t          n_g_new_vertices)
{
  cs_int_t  i, j, k, shift, vid, start, end, perio_start, perio_end;
  cs_int_t  n_new_vertices, n_init_faces;
  cs_real_t  matrix[3][4], xyz[4];
  cs_bool_t  is_modified;
  fvm_gnum_t  new_gnum;
  cs_join_state_t  state;

  cs_int_t  *new_f2v_idx = NULL, *new_f2v_lst = NULL, *vtag = NULL;
  cs_int_t  *linked_id = NULL;
  fvm_gnum_t  *gnum = NULL;
  cs_bool_t  *f_state = NULL;
  cs_join_mesh_t  *work_jmesh = *p_work_jmesh;
  cs_join_edges_t  *work_edges = *p_work_edges;
  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *select = this_join->selection;
  cs_join_perio_builder_t  *builder = cs_glob_join_perio_builder;

  const int  n_ranks = cs_glob_n_ranks;
  const int  perio_id = param.perio_num - 1;

  assert(perio_id > -1);
  assert(builder != NULL);

  /* Retrieve related back transformation */

  fvm_periodicity_get_matrix(builder->periodicity, 2*perio_id+1, matrix);

  BFT_MALLOC(linked_id, jmesh->n_vertices, cs_int_t);
  BFT_MALLOC(gnum, jmesh->n_vertices, fvm_gnum_t);

  for (i = 0; i < jmesh->n_vertices; i++) {
    linked_id[i] = -1; /* Default: no link */
    gnum[i] = jmesh->vertices[i].gnum; /* ordered list */
  }

  /* Linked vertex id by inverse transformation */

  for (i = 0; i < select->n_couples; i++) {

    j = cs_search_g_binary(jmesh->n_vertices,
                           select->per_v_couples[2*i],
                           gnum);

    k = cs_search_g_binary(jmesh->n_vertices,
                           select->per_v_couples[2*i+1],
                           gnum);

    linked_id[k] = j;

  }

  BFT_FREE(gnum);

  /* Scan faces to detect new vertices in order to define and build a
     new face->vertex index */

  n_init_faces = jmesh->n_faces/2;

  BFT_MALLOC(f_state, jmesh->n_faces, cs_bool_t);
  BFT_MALLOC(new_f2v_idx, jmesh->n_faces + 1, cs_int_t);
  BFT_MALLOC(vtag, jmesh->n_vertices, cs_int_t);

  for (i = 0; i < jmesh->n_vertices; i++)
    vtag[i] = 0;

  for (i = 0; i < n_init_faces; i++) {

    is_modified = false;
    start = jmesh->face_vtx_idx[2*i]-1;
    end = jmesh->face_vtx_idx[2*i+1]-1;
    perio_start = end;
    perio_end = jmesh->face_vtx_idx[2*i+2]-1;

    for (j = perio_start; j < perio_end; j++) {

      vid = jmesh->face_vtx_lst[j] - 1;
      state = jmesh->vertices[vid].state;

      if (state == CS_JOIN_STATE_PERIO_MERGE) {

        is_modified = true;
        vtag[vid] = -1;
        assert(linked_id[vid] > -1);

      }
      else if (state == CS_JOIN_STATE_NEW) {

        is_modified = true;
        vtag[vid] = 1;
        assert(linked_id[vid] == -1);

      }
      else if (state == CS_JOIN_STATE_MERGE) {

        is_modified = true;
        if (linked_id[vid] > -1) /* Update is enough */
          vtag[vid] = -2;

        else /* New vertex for the periodic face but not for
                  the original intersected face.
                  Add a new vertex for the related original face */

          vtag[vid] = 2;

      }

    }

    if (is_modified == true)
      new_f2v_idx[2*i+1] = perio_end - perio_start;
    else
      new_f2v_idx[2*i+1] = end - start;
    new_f2v_idx[2*i+2] = perio_end - perio_start;

    f_state[2*i] = is_modified;
    f_state[2*i+1] = false;

  }

  n_new_vertices = 0;
  for (i = 0; i < jmesh->n_vertices; i++)
    if (vtag[i] > 0)
      n_new_vertices++;

  BFT_REALLOC(jmesh->vertices,
              jmesh->n_vertices + n_new_vertices,
              cs_join_vertex_t);

  /* Transform back new periodic vertices */

  n_new_vertices = 0;

  for (i = 0; i < jmesh->n_vertices; i++) {

    if (vtag[i] > 0) { /* New vertex to transform back */

      cs_join_vertex_t  new_vtx = jmesh->vertices[i];

      for (j = 0; j < 3; j++) {
        xyz[j] = new_vtx.coord[j];
        new_vtx.coord[j] = 0.0;
      }
      xyz[3] = 1;

      for (j = 0; j < 3; j++)
        for (k = 0; k < 4; k++)
          new_vtx.coord[j] += matrix[j][k]*xyz[k];

      jmesh->vertices[jmesh->n_vertices + n_new_vertices] = new_vtx;
      n_new_vertices++;
      vtag[i] = jmesh->n_vertices + n_new_vertices;

    }
    else if (vtag[i] < 0) { /* Existing vertex to update */

      cs_join_vertex_t  new_vtx = jmesh->vertices[linked_id[i]];

      assert(   new_vtx.state != CS_JOIN_STATE_MERGE
             || new_vtx.state != CS_JOIN_STATE_PERIO_MERGE);

      for (j = 0; j < 3; j++) {
        xyz[j] = jmesh->vertices[i].coord[j];
        new_vtx.coord[j] = 0.0;
      }
      xyz[3] = 1;

      for (j = 0; j < 3; j++)
        for (k = 0; k < 4; k++)
          new_vtx.coord[j] += matrix[j][k]*xyz[k];

      new_vtx.state = CS_JOIN_STATE_MERGE;
      jmesh->vertices[linked_id[i]] = new_vtx;

    }

  }

  if (n_ranks > 1) { /* Global numbering update */

    fvm_io_num_t  *io_num = NULL;
    const fvm_gnum_t  *io_gnum = NULL;

    BFT_MALLOC(gnum, n_new_vertices, fvm_gnum_t);

    for (i = 0, shift = jmesh->n_vertices; i < n_new_vertices; i++, shift++)
      gnum[i] = jmesh->vertices[shift].gnum;

    io_num = fvm_io_num_create(NULL, gnum, n_new_vertices, 0);
    io_gnum = fvm_io_num_get_global_num(io_num);

    for (i = 0, shift = jmesh->n_vertices; i < n_new_vertices; i++, shift++) {
      new_gnum = io_gnum[i] + init_max_vtx_gnum + n_g_new_vertices;
      jmesh->vertices[shift].gnum = new_gnum;
    }

    jmesh->n_g_vertices += fvm_io_num_get_global_count(io_num);

    fvm_io_num_destroy(io_num);
    BFT_FREE(gnum);

  }
  else { /* Serial run */

    for (i = 0, shift = jmesh->n_vertices; i < n_new_vertices; i++, shift++) {
      new_gnum = i + 1 + init_max_vtx_gnum + n_g_new_vertices;
      jmesh->vertices[shift].gnum = new_gnum;
    }

    jmesh->n_g_vertices += n_new_vertices;

  }

  /* Update face->vertex connectivity for original faces if needed
     Copy connectivity for periodic faces */

  new_f2v_idx[0] = 1;
  for (i = 0; i < jmesh->n_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  BFT_MALLOC(new_f2v_lst, new_f2v_idx[jmesh->n_faces] - 1, cs_int_t);

  for (i = 0; i < n_init_faces; i++) {

    start = jmesh->face_vtx_idx[2*i]-1;
    end = jmesh->face_vtx_idx[2*i+1]-1;
    perio_start = end;
    perio_end = jmesh->face_vtx_idx[2*i+2]-1;

    if (f_state[2*i] == false) { /* No modification to apply */

      for (j = start, shift = new_f2v_idx[2*i] - 1;
           j < perio_start; j++, shift++)
        new_f2v_lst[shift] = jmesh->face_vtx_lst[j];

    }
    else { /* Modification to apply from the periodic face */

      for (j = perio_start, shift = new_f2v_idx[2*i] - 1;
           j < perio_end; j++, shift++) {

        vid = jmesh->face_vtx_lst[j] - 1;
        state = jmesh->vertices[vid].state;

        if (   state == CS_JOIN_STATE_PERIO_MERGE
            || state == CS_JOIN_STATE_PERIO)
          new_f2v_lst[shift] = linked_id[vid] + 1;
        else if (state == CS_JOIN_STATE_MERGE) {
          if (linked_id[vid] > -1)
            new_f2v_lst[shift] = linked_id[vid] + 1;
          else
            new_f2v_lst[shift] = vtag[vid];
        }
        else if (state == CS_JOIN_STATE_NEW)
          new_f2v_lst[shift] = vtag[vid];
        else
          bft_error(__FILE__, __LINE__, 0,
                    _("  Vertex state (%d) is not consistent.\n"
                      "  Can not apply changes from periodic faces to"
                      " original face.\n"
                      "  Check your periodicity parameters.\n"), state);
      }

    } /* End of test on face state */

    /* Copy periodic face connectivity */

    for (j = perio_start, shift = new_f2v_idx[2*i+1] - 1;
         j < perio_end; j++, shift++)
      new_f2v_lst[shift] = jmesh->face_vtx_lst[j];

  } /* End of loop on initial faces */

  /* Add new vertices to the periodic vertex list */

  if (param.verbosity > 2)
  bft_printf("  Add locally %d new vertices for periodicity\n",
             n_new_vertices);

  shift = select->n_couples;
  select->n_couples += n_new_vertices;
  BFT_REALLOC(select->per_v_couples, 2*select->n_couples, fvm_gnum_t);

  for (i = 0; i < jmesh->n_vertices; i++) {
    if (vtag[i] > 0) {
      vid = vtag[i] - 1;
      select->per_v_couples[2*shift] = jmesh->vertices[vid].gnum;
      select->per_v_couples[2*shift+1] = jmesh->vertices[i].gnum;
      shift++;
    }
  }

  BFT_FREE(vtag);
  BFT_FREE(linked_id);
  BFT_FREE(f_state);

  /* Reshape join_mesh structure */

  BFT_FREE(jmesh->face_vtx_idx);
  BFT_FREE(jmesh->face_vtx_lst);

  jmesh->face_vtx_idx = new_f2v_idx;
  jmesh->face_vtx_lst = new_f2v_lst;
  jmesh->n_vertices += n_new_vertices;

  /* Update now work_jmesh by exchanging jmesh over the ranks */

  _redistribute_mesh(param,
                     select->compact_rank_index,
                     jmesh,
                     &work_jmesh);

  /* Define a new cs_join_edges_t structure related to work_jmesh */

  cs_join_mesh_destroy_edges(&work_edges);
  work_edges = cs_join_mesh_define_edges(work_jmesh);

  /* Reshape join_mesh structure by deleting periodic faces */

  shift = 0;

  for (i = 0; i < n_init_faces; i++) {

    jmesh->face_gnum[i] = jmesh->face_gnum[2*i];

    for (j = jmesh->face_vtx_idx[2*i]-1; j < jmesh->face_vtx_idx[2*i+1]-1; j++)
      jmesh->face_vtx_lst[shift++] = jmesh->face_vtx_lst[j];

    jmesh->face_vtx_idx[i+1] = shift + 1;

  }

  BFT_REALLOC(jmesh->face_gnum, n_init_faces, fvm_gnum_t);
  BFT_REALLOC(jmesh->face_vtx_idx, n_init_faces + 1, cs_int_t);
  BFT_REALLOC(jmesh->face_vtx_lst, shift, cs_int_t);

  jmesh->n_faces = n_init_faces;
  jmesh->n_g_faces /= 2;

  /* Return pointer */

  *p_work_jmesh = work_jmesh;
  *p_work_edges = work_edges;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n Periodic vertex couples:\n");
  for (i = 0; i < select->n_couples; i++)
    bft_printf(" %6d  (%9u, %9u)\n", i+1,
               select->per_v_couples[2*i], select->per_v_couples[2*i+1]);
  bft_printf_flush();
#endif

}

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Update jmesh structure.
 * Define a new n2o_hist.
 *
 * parameters:
 *   this_join  <-- pointer to a high level join structure
 *   jmesh      <-> local join mesh struct. to duplicate and transform
 *   mesh       <-- pointer to a cs_mesh_t structure
 *   o2n_hist   <-- old global face -> new local face numbering
 *   p_n2o_hist <-- new global face -> old local face numbering
 *---------------------------------------------------------------------------*/

void
cs_join_perio_split_back(cs_join_t          *this_join,
                         cs_join_mesh_t     *jmesh,
                         cs_mesh_t          *mesh,
                         cs_join_gset_t     *o2n_hist,
                         cs_join_gset_t    **p_n2o_hist)
{
  int  i, j, k, shift, vid, fid, start, end, perio_start, perio_end;
  int  n_final_faces, n1_faces, n2_faces;
  int  shift1, shift2, shift3, shift4;
  int  n_sub_ori, n_sub_per, n_contrib, n_couples;
  fvm_gnum_t  n2_g_faces;

  int  n_vertices_to_add = 0, n_g_vertices_to_add = 0;
  cs_int_t  *new_f2v_idx = NULL, *new_f2v_lst = NULL;
  cs_int_t  *linked_id = NULL, *f_tag = NULL;
  fvm_gnum_t  *gnum = NULL, *f2_gnum = NULL, *new_fgnum = NULL;
  cs_join_gset_t  *new_history = NULL, *n2o_hist = *p_n2o_hist;

  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *select = this_join->selection;
  cs_join_perio_builder_t  *builder = cs_glob_join_perio_builder;

  const int  n_ranks = cs_glob_n_ranks;
  const int  perio_id = param.perio_num - 1;

  assert(perio_id > -1);
  assert(builder != NULL);

  /* Detect periodic face to delete and associate a tag for each new face */

  BFT_MALLOC(f_tag, jmesh->n_faces, cs_int_t);

  assert(n2o_hist->n_elts == jmesh->n_faces);
  n_couples = 0;

  for (i = 0; i < n2o_hist->n_elts; i++) {

    start = n2o_hist->index[i];
    end = n2o_hist->index[i+1];
    n_contrib = end - start;

    if (n_contrib == 1) { /* New face remains a border face */

      if (n2o_hist->g_list[start]%2 == 0)
        f_tag[i] = 0; /* To delete without transferring back face connect. */
      else
        f_tag[i] = 1; /* To keep. Original face */

    }
    else if (n_contrib == 2) { /* New face becomes an interior face */

      if (n2o_hist->g_list[start]%2 == 0) { /* Old periodic face */

        assert(n2o_hist->g_list[start+1]%2 == 1); /* Associated face should
                                                     be an old original face */

        f_tag[i] = -1; /* To keep with transferring back face connect. */
        n_couples++;

      }
      else { /* Old original face */

        assert(n2o_hist->g_list[start]%2 == 1);
        assert(n2o_hist->g_list[start+1]%2 == 0); /* Associated face should
                                                     be an old periodic face */

        f_tag[i] = -1; /* To keep with transferring back face connect. */
        n_couples++;

      }

    }
    else { /* n_contrib > 2 => remains a border face and keep the new face
              if not all the old faces are periodic */

      cs_bool_t  have_perio = true;

      f_tag[i] = 0; /* Initialize as if we want to delete the new face */

      for (j = start; j < end; j++)
        if (n2o_hist->g_list[j]%2 == 1) /* Old original face are taking part
                                           in the composition of the new
                                           face */
          have_perio = false;

      if (have_perio == false)
        f_tag[i] = 1; /* To keep without transferring back face connect. */
    }

  } /* End of loop on new faces */

  /* Loop over old faces and change the current tag for the new faces
     according to the number of subdivisions applied to the original
     faces */

  for (i = 0; i < select->n_faces; i++) {

    start = o2n_hist->index[2*i];
    end = o2n_hist->index[2*i+1];
    perio_start = end;
    perio_end = o2n_hist->index[2*i+2];

    n_sub_ori = end - start;
    n_sub_per = perio_end - perio_start;

    assert(n_sub_per > 0 && n_sub_ori > 0);

    if (n_sub_ori == 1 && n_sub_per > 1) {

      fid = cs_search_g_binary(jmesh->n_faces,
                               o2n_hist->g_list[start],
                               jmesh->face_gnum);

      f_tag[fid] = 2; /* Original face to delete and to replace by
                         a set of new sub-faces */

      for (j = perio_start; j < perio_end; j++) {

        fid = cs_search_g_binary(jmesh->n_faces,
                                 o2n_hist->g_list[j],
                                 jmesh->face_gnum);

        assert(fid > -1);

        if (f_tag[fid] == 0)
          f_tag[fid] = -2; /* Switch tag:
                              Will be deleted after periodicity application */

      }

    } /* n_sub_ori == 1 && n_sub_per > 1 */

    else if (n_sub_ori == 1 && n_sub_per == 1) {

      fid = cs_search_g_binary(jmesh->n_faces,
                               o2n_hist->g_list[perio_start],
                               jmesh->face_gnum);

      if (n2o_hist->index[fid+1] - n2o_hist->index[fid] == 2) {

        assert(f_tag[fid] == -1); /* New face to keep with transferring
                                     back face connectivity */

        fid = cs_search_g_binary(jmesh->n_faces,
                                 o2n_hist->g_list[start],
                                 jmesh->face_gnum);


        assert(fid > -1);
        f_tag[fid] = 2; /* Original face to delete and to replace by
                           a set the new sub-face */

      }

    }

  } /* End of loop on old faces */

  /* Count the final number of new subfaces */

  n_final_faces = 0;
  for (i = 0; i < jmesh->n_faces; i++) {

    if (f_tag[i] == 1 || f_tag[i] == -2)
      n_final_faces++;
    else if (f_tag[i] == -1)
      n_final_faces += 2;

  }

  /* Define a global number for each face to transfer back by periodicity */

  n2_faces = 0;
  for (i = 0; i < jmesh->n_faces; i++)
    if (f_tag[i] < 0)
      n2_faces++;
  n2_g_faces = n2_faces;

  BFT_MALLOC(f2_gnum, n2_faces, fvm_gnum_t);

  if (n_ranks > 1) { /* Parallel run */

    fvm_io_num_t  *io_num = NULL;
    const fvm_gnum_t  *io_gnum = NULL;

    n2_faces = 0;
    for (i = 0; i < jmesh->n_faces; i++)
      if (f_tag[i] < 0)
        f2_gnum[n2_faces++] = jmesh->face_gnum[i];

    io_num = fvm_io_num_create(NULL, f2_gnum, n2_faces, 0);
    io_gnum = fvm_io_num_get_global_num(io_num);

    for (i = 0; i < n2_faces; i++)
      f2_gnum[i] = jmesh->n_g_faces + io_gnum[i];

    n2_g_faces = fvm_io_num_get_global_count(io_num);

    fvm_io_num_destroy(io_num);

  }
  else { /* Serial run */

    n2_faces = 0;
    for (i = 0; i < jmesh->n_faces; i++) {
      if (f_tag[i] < 0) {
        f2_gnum[n2_faces] = jmesh->n_faces + n2_faces + 1;
        n2_faces++;
      }
    }

  }

  /* Define the new face -> vertex index and the new global face numbering */

  BFT_MALLOC(new_f2v_idx, n_final_faces + 1, cs_int_t);
  BFT_MALLOC(new_fgnum, n_final_faces, fvm_gnum_t);

  new_history = cs_join_gset_create(n_final_faces);
  new_history->n_g_elts = new_history->n_elts;

  n1_faces = 0;

  for (i = 0; i < jmesh->n_faces; i++) {
    if (f_tag[i] == 1 || f_tag[i] == -1) {

      new_f2v_idx[n1_faces+1] =
        jmesh->face_vtx_idx[i+1] - jmesh->face_vtx_idx[i];
      new_fgnum[n1_faces] = jmesh->face_gnum[i];

      new_history->index[n1_faces+1] =
        n2o_hist->index[i+1] - n2o_hist->index[i];

      n1_faces++;

    }
  }

  shift = n1_faces;
  n2_faces = 0;
  for (i = 0; i < jmesh->n_faces; i++) {
    if (f_tag[i] < 0) {

      assert(f_tag[i] == -1 || f_tag[i] == -2);

      new_f2v_idx[shift + 1] =
        jmesh->face_vtx_idx[i+1] - jmesh->face_vtx_idx[i];
      new_fgnum[shift] = f2_gnum[n2_faces];

      new_history->index[shift+1] =
        n2o_hist->index[i+1] - n2o_hist->index[i];

      n2_faces++;
      shift++;

    }
  }

  assert(n1_faces + n2_faces == n_final_faces);

  /* Detect if there are new vertices to add to the jmesh definition */

  BFT_MALLOC(gnum, jmesh->n_vertices, fvm_gnum_t);
  BFT_MALLOC(linked_id, jmesh->n_vertices, cs_int_t);

  for (i = 0; i < jmesh->n_vertices; i++) {
    linked_id[i] = -1; /* Default: no link */
    gnum[i] = jmesh->vertices[i].gnum; /* ordered list */
  }

  for (i = 0; i < select->n_couples; i++) {

    j = cs_search_g_binary(jmesh->n_vertices,
                           select->per_v_couples[2*i],
                           gnum);

    k = cs_search_g_binary(jmesh->n_vertices,
                           select->per_v_couples[2*i+1],
                           gnum);

    linked_id[k] = j;

  }

  BFT_FREE(f2_gnum);
  BFT_FREE(gnum);

  for (i = 0; i < jmesh->n_faces; i++) {
    if (f_tag[i] < 0) { /* Periodicity to transfer back */

      for (j = jmesh->face_vtx_idx[i]-1; j < jmesh->face_vtx_idx[i+1]-1; j++) {
        vid = jmesh->face_vtx_lst[j] - 1;

        if (linked_id[vid] == -1)
          linked_id[vid] = -2; /* Have to get the related id. Add a new
                                  vertex */

      }

    }
  } /* End of loop on new faces */

  n_vertices_to_add = 0;
  for (i = 0; i < jmesh->n_vertices; i++)
    if (linked_id[i] == -2)
      n_vertices_to_add += 1;

  if (n_vertices_to_add > 0) {

    cs_int_t  i1, i2;
    cs_real_t  matrix[3][4], xyz[4];

    BFT_REALLOC(jmesh->vertices, jmesh->n_vertices + n_vertices_to_add,
                cs_join_vertex_t);
    BFT_REALLOC(linked_id,  jmesh->n_vertices + n_vertices_to_add, cs_int_t);

    /* Retrieve related back transformation */

    fvm_periodicity_get_matrix(builder->periodicity, 2*perio_id+1, matrix);

    n_vertices_to_add = 0;
    for (i1 = 0; i1 < jmesh->n_vertices; i1++) {
      if (linked_id[i1] == -2) {

        cs_join_vertex_t  new_vtx = jmesh->vertices[i1];

        i2 = jmesh->n_vertices + n_vertices_to_add;
        linked_id[i1] = i2;

        for (k = 0; k < 3; k++) {
          xyz[k] = new_vtx.coord[k];
          new_vtx.coord[k] = 0.0;
        }
        xyz[3] = 1;

        for (j = 0; j < 3; j++)
          for (k = 0; k < 4; k++)
            new_vtx.coord[j] += matrix[j][k]*xyz[k];

        jmesh->vertices[i2] = new_vtx;
        n_vertices_to_add += 1;

      }
    }

  } /* n_vertices_to_add > 0 */

  /* Define a global vertex num for the new added vertices */

#if defined(HAVE_MPI)
  if (n_ranks > 1) {

    fvm_io_num_t  *io_num = NULL;
    const fvm_gnum_t  *io_gnum = NULL;

    MPI_Allreduce(&n_vertices_to_add, &n_g_vertices_to_add, 1, MPI_INT,
                  MPI_SUM, cs_glob_mpi_comm);

    /* Define global numbering for vertices to add
       and keep relation between periodic couples */

    if (n_g_vertices_to_add > 0) {

      BFT_MALLOC(gnum, n_vertices_to_add, fvm_gnum_t);

      for (i = 0, shift = jmesh->n_vertices; i < n_vertices_to_add;
           i++, shift++)
        gnum[i] = jmesh->vertices[shift].gnum;

      io_num = fvm_io_num_create(NULL, gnum, n_vertices_to_add, 0);
      io_gnum = fvm_io_num_get_global_num(io_num);

      for (i = 0, shift = jmesh->n_vertices; i < n_vertices_to_add;
           i++, shift++)
        jmesh->vertices[shift].gnum = io_gnum[i] + mesh->n_g_vertices;

      n_g_vertices_to_add = fvm_io_num_get_global_count(io_num);
      jmesh->n_g_vertices += n_g_vertices_to_add;

      fvm_io_num_destroy(io_num);
      BFT_FREE(gnum);

    }

  }
#endif

  if (n_ranks == 1 && n_vertices_to_add > 0) {

    for (i = 0, shift = jmesh->n_vertices; i < n_vertices_to_add; i++, shift++)
      jmesh->vertices[shift].gnum = i + 1 + mesh->n_g_vertices;

    jmesh->n_g_vertices += n_vertices_to_add;

  }

  jmesh->n_vertices += n_vertices_to_add;

  new_f2v_idx[0] = 1;
  new_history->index[0] = 0;

  for (i = 0; i < n_final_faces; i++) {
    new_f2v_idx[i+1] += new_f2v_idx[i];
    new_history->index[i+1] += new_history->index[i];
  }

  BFT_MALLOC(new_f2v_lst, new_f2v_idx[n_final_faces] - 1, cs_int_t);
  BFT_MALLOC(new_history->g_list, new_history->index[new_history->n_elts],
             fvm_gnum_t);

  /* Define a new face connectivity and a new face history */

  shift1 = 0; /*kept faces in face connect. */
  shift2 = new_f2v_idx[n1_faces] - 1; /* transfered faces in face connect. */

  shift3 = 0; /* kept faces in face history */
  shift4 = new_history->index[n1_faces]; /* transfered faces in face history */

  /* Store periodic couples. First in local join numbering.
     Will move next to a global numbering (after interior face add) */

  builder->n_perio_couples[perio_id] = n_couples;
  BFT_MALLOC(builder->perio_couples[perio_id], 2*n_couples, fvm_gnum_t);

  n2_faces = n1_faces;
  n1_faces = 0;
  n_couples = 0;

  for (i = 0; i < jmesh->n_faces; i++) {

    start = jmesh->face_vtx_idx[i]-1;
    end = jmesh->face_vtx_idx[i+1]-1;

    if (f_tag[i] == 1) { /* Original face to keep */

      for (j = start; j < end; j++)
        new_f2v_lst[shift1++] = jmesh->face_vtx_lst[j];

      for (j = n2o_hist->index[i]; j < n2o_hist->index[i+1]; j++)
        new_history->g_list[shift3++] = n2o_hist->g_list[j];

      n1_faces++;

    }
    else if (f_tag[i] == -1) { /* Periodic face to keep and to transfer back */

      /* Store periodic face couple */

      n1_faces++;
      n2_faces++;
      builder->perio_couples[perio_id][2*n_couples] = n2_faces;
      builder->perio_couples[perio_id][2*n_couples+1] = n1_faces;
      n_couples++;

      /* Define face connectivity */

      for (j = start; j < end; j++) {
        vid = jmesh->face_vtx_lst[j] - 1;
        assert(linked_id[vid] > -1);
        new_f2v_lst[shift1++] = vid + 1;
        new_f2v_lst[shift2++] = linked_id[vid] + 1;
      }

      /* Define new old->new face history */

      for (j = n2o_hist->index[i]; j < n2o_hist->index[i+1]; j++) {

        if (n2o_hist->g_list[j]%2 == 0) { /* Periodic face */
          new_history->g_list[shift3++] = n2o_hist->g_list[j];
          /* Original face related to this periodic one */
          new_history->g_list[shift4++] = n2o_hist->g_list[j] - 1;
        }
        else { /* Original face */
          new_history->g_list[shift3++] = n2o_hist->g_list[j];
          /* Periodic face related to this original one */
          new_history->g_list[shift4++] = n2o_hist->g_list[j] + 1;
        }

      }

    }
    else if (f_tag[i] == -2) { /* Periodic face to transfer back */

      n2_faces++;

      for (j = start; j < end; j++) {
        vid = jmesh->face_vtx_lst[j] - 1;
        assert(linked_id[vid] > -1);
        new_f2v_lst[shift2++] = linked_id[vid] + 1;
      }

      for (j = n2o_hist->index[i]; j < n2o_hist->index[i+1]; j++) {

        assert(n2o_hist->g_list[j]%2 == 0); /* Periodic face */

        /* Original face related to this periodic one */
        new_history->g_list[shift4++] = n2o_hist->g_list[j] - 1;

      }

    }

  } /* End of loop on faces */

  BFT_FREE(linked_id);
  BFT_FREE(f_tag);

  /* Reshape join_mesh structure */

  jmesh->n_faces = n_final_faces;
  jmesh->n_g_faces = jmesh->n_faces;

  BFT_FREE(jmesh->face_gnum);
  BFT_FREE(jmesh->face_vtx_idx);
  BFT_FREE(jmesh->face_vtx_lst);

  jmesh->face_gnum = new_fgnum;
  jmesh->face_vtx_idx = new_f2v_idx;
  jmesh->face_vtx_lst = new_f2v_lst;

  if (n_ranks > 1) { /* Update global face number */

    fvm_io_num_t  *io_num = NULL;
    const fvm_gnum_t  *io_gnum = NULL;

    io_num = fvm_io_num_create(NULL, jmesh->face_gnum, jmesh->n_faces, 0);
    io_gnum = fvm_io_num_get_global_num(io_num);

    for (i = 0; i < jmesh->n_faces; i++) {
      jmesh->face_gnum[i] = io_gnum[i];
      new_history->g_elts[i] = io_gnum[i];
    }

    jmesh->n_g_faces = fvm_io_num_get_global_count(io_num);
    new_history->n_g_elts = jmesh->n_g_faces;

    fvm_io_num_destroy(io_num);
    BFT_FREE(gnum);

  }
  else { /* Serial run */

    for (i = 0; i < jmesh->n_faces; i++) {
      jmesh->face_gnum[i] = i+1;
      new_history->g_elts[i] = i+1;
    }

  }

  /* Return pointer */

  cs_join_gset_destroy(&n2o_hist);

  *p_n2o_hist = new_history;

}

/*----------------------------------------------------------------------------
 * Define a list of coupled faces by periodicty in global numbering.
 * For parallel runs:
 *  - remove isolated periodic faces in the mesh definition
 *  - define a consistent face connectivity in order to prepare the building
 *    of periodic vertex couples
 *
 *
 * parameters:
 *   param      <-- set of parameters for the joining operation
 *   n_ii_faces <-- initial local number of interior faces
 *   face_type  <-- type of faces in join mesh (interior or border ...)
 *   jmesh      <-- pointer on a cs_join_mesh_t struct.
 *   mesh       <-> pointer on a cs_mesh_t struct.
 *---------------------------------------------------------------------------*/

void
cs_join_perio_split_update(cs_join_param_t             param,
                           cs_int_t                    n_ii_faces,
                           const cs_join_face_type_t   face_type[],
                           const cs_join_mesh_t       *jmesh,
                           cs_mesh_t                  *mesh)
{
  int  i, shift;

  fvm_gnum_t  *o2n_num = NULL;
  cs_join_perio_builder_t  *builder = cs_glob_join_perio_builder;

  const int  n_j_faces = jmesh->n_faces;
  const int  n_ranks = cs_glob_n_ranks;
  const int  perio_id = param.perio_num - 1;

  assert(perio_id > -1);
  assert(builder != NULL);

  /* Initialize o2n_num */

  BFT_MALLOC(o2n_num, n_j_faces, fvm_gnum_t);

  for (i = 0; i < n_j_faces; i++)
    o2n_num[i] = 0;

  if (n_ranks == 1) {

    shift = n_ii_faces + 1;
    for (i = 0; i < n_j_faces; i++)
      if (face_type[i] == CS_JOIN_FACE_INTERIOR)
        o2n_num[i] = shift++;

  }
  else { /* n_ranks > 1 */

    shift = n_ii_faces;
    for (i = 0; i < n_j_faces; i++)
      if (face_type[i] == CS_JOIN_FACE_INTERIOR)
        o2n_num[i] = mesh->global_i_face_num[shift++];

  }

  /* Apply new numbering */

  for (i = 0; i < builder->n_perio_couples[perio_id]; i++) {

    fvm_gnum_t  old1 = builder->perio_couples[perio_id][2*i] - 1;
    fvm_gnum_t  old2 = builder->perio_couples[perio_id][2*i+1] - 1;

    assert(o2n_num[old1] > 0);
    assert(o2n_num[old2] > 0);

    builder->perio_couples[perio_id][2*i] = o2n_num[old1];
    builder->perio_couples[perio_id][2*i+1] = o2n_num[old2];

  }

  BFT_FREE(o2n_num);

  if (n_ranks > 1) /* Remove isolated periodic face for the current mesh */
    _perio_face_clean(param, mesh);

}

/*----------------------------------------------------------------------------
 * Use periodic face couples in cs_glob_join_perio_builder to define
 * cs_glob_mesh_builder.
 * Free all elements which can be freed.
 * Transfer data to cs_glob_mesh and cs_glob_mesh_builder.
 *---------------------------------------------------------------------------*/

void
cs_join_perio_transfer_builder(void)
{
  int  i, j, k;

  cs_join_perio_builder_t  *p_builder = cs_glob_join_perio_builder;
  cs_mesh_builder_t  *m_builder = cs_glob_mesh_builder;
  cs_mesh_t  *mesh = cs_glob_mesh;

  const int  n_perio = cs_glob_n_join_perio;
  const int  n_ranks = cs_glob_n_ranks;

  /* Periodicity is either defined by a joining operation, either defined
     by a preprocessor command. One cannot use the way together */

  assert(m_builder != NULL);
  assert(m_builder->per_face_idx == NULL);
  assert(m_builder->per_face_lst == NULL);
  assert(m_builder->per_rank_lst == NULL);
  assert(p_builder != NULL);
  assert(mesh != NULL);
  assert(n_perio == p_builder->n_perio);
  assert(mesh->periodicity == NULL); /* No periodicity already defined by
                                        the preprocessor */

  /* Transfer data to mesh structure */

  mesh->n_init_perio = n_perio;
  mesh->periodicity = p_builder->periodicity;

  /* Define mesh builder */

  BFT_MALLOC(m_builder->per_face_idx, n_perio + 1, cs_int_t);

  m_builder->per_face_idx[0] = 0;
  for (i = 0; i < n_perio; i++)
    m_builder->per_face_idx[i+1] =  m_builder->per_face_idx[i]
                                  + p_builder->n_perio_couples[i];

  BFT_MALLOC(m_builder->per_face_lst, 2*m_builder->per_face_idx[n_perio],
             cs_int_t);

  if (n_ranks == 1) {

    for (i = 0, k = 0; i < n_perio; i++) {
      for (j = 0; j < p_builder->n_perio_couples[i]; j++, k++) {
        m_builder->per_face_lst[2*k] = p_builder->perio_couples[i][2*j];
        m_builder->per_face_lst[2*k+1] = p_builder->perio_couples[i][2*j+1];
      }
    }

  }
  else { /* Parallel run: n_ranks > 1 */

    fvm_interface_set_t *face_ifs = NULL;
    fvm_lnum_t  *periodicity_num = NULL;

    BFT_MALLOC(periodicity_num, n_perio, fvm_lnum_t);
    BFT_MALLOC(m_builder->per_rank_lst, m_builder->per_face_idx[n_perio],
               cs_int_t);

    for (i = 0; i < n_perio; i++)
      periodicity_num[i] = i+1;

    face_ifs = fvm_interface_set_create(mesh->n_i_faces,
                                        NULL,
                                        mesh->global_i_face_num,
                                        mesh->periodicity,
                                        n_perio,
                                        periodicity_num,
                                        p_builder->n_perio_couples,
              (const fvm_gnum_t **const)p_builder->perio_couples);

#if defined(HAVE_MPI) /* Algo. similar to the one in cs_preprocessor_data.c */
    _extract_perio_couples(n_perio, m_builder, face_ifs);
#endif

    /* Keep face interface to build the future i_face_cells connect. */

    m_builder->face_ifs = face_ifs;

    /* Free memory */

    BFT_FREE(periodicity_num);

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
 {
   int  perio_id;

   const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

   bft_printf("\n  Dump periodic data built from the joining algorithm\n");

   for (perio_id = 0; perio_id < n_perio; perio_id++) {

     int  start_id = m_builder->per_face_idx[perio_id];
     int  end_id = m_builder->per_face_idx[perio_id+1];

     bft_printf("\n  Perio id: %4d - Number of elements: %7d "
                "(start: %7d - end: %7d)\n",
                perio_id, end_id-start_id, start_id, end_id);
     bft_printf("   id    | 1st face | 2nd face | associated rank\n");

     for (i = start_id; i < end_id; i++) {
       if (cs_glob_n_ranks > 1) {

         int  f1_id = CS_ABS(m_builder->per_face_lst[2*i]) -1;
         int  f2_id = CS_ABS(m_builder->per_face_lst[2*i+1]) -1;

         bft_printf("%8d | %8d (%9u) | %8d (%9u) | %5d\n",
                    i, m_builder->per_face_lst[2*i],
                    cs_glob_mesh->global_i_face_num[f1_id],
                    m_builder->per_face_lst[2*i+1],
                    (m_builder->per_rank_lst[i]-1 == local_rank ?
                     cs_glob_mesh->global_i_face_num[f2_id] : 0),
                    m_builder->per_rank_lst[i]-1);

       }
       else
         bft_printf("%8d | %10d | %10d | %6d\n",
                    i, m_builder->per_face_lst[2*i],
                    m_builder->per_face_lst[2*i+1],
                    local_rank);
     }
     bft_printf_flush();

   }

 }
#endif

  _delete_perio_builder(&p_builder);

}

/*---------------------------------------------------------------------------*/

END_C_DECLS

