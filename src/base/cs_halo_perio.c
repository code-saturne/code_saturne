/*============================================================================
 * Functions handling periodicity.
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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_mesh.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_halo_perio.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure definitions
 *============================================================================*/

/* Structure used for building mesh structure */
/* ------------------------------------------ */

typedef struct {

  /* Periodic features */

  cs_lnum_t  *per_face_idx;    /* Index on periodicity for per_face_lst */

  cs_lnum_t  *per_face_lst;    /* Periodic faces list. For each couple,
                                  we have the local face number on local rank
                                  and the local face number on distant rank */

  cs_lnum_t  *per_rank_lst;    /* Remote ranks list. For each couple,
                                  we have the distant rank number. Exists
                                  only in case of parallelism. */

} _perio_mesh_builder_t ;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Table giving the Reynolds stress component for [i][j] */

/* Warning: old fashion to store Rij */
static const int _rij[3][3] = {{0, 3, 4},
                               {3, 1, 5},
                               {4, 5, 2}};

static const int _symt[3][3] = {{0, 3, 5},
                                {3, 1, 4},
                                {5, 4, 2}};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a transformation to a given
 * vector.
 *
 * parameters:
 *   matrix[3][4] --> matrix of the transformation in homogeneous coord.
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   in_cell_id   --> cell_id of the parent cell.
 *   out_cell_id  --> cell_id of the generated cell.
 *   xyz          <-> array of coordinates
 *----------------------------------------------------------------------------*/

static void
_apply_vector_transfo(cs_real_t    matrix[3][4],
                      cs_lnum_t    in_cell_id,
                      cs_lnum_t    out_cell_id,
                      cs_real_t   *xyz)
{
  cs_lnum_t  i, j;

  cs_real_t  xyz_a[3 + 1];
  cs_real_t  xyz_b[3];

  /* Define the cell center in homogeneous coordinates before
     transformation */

  for (j = 0; j < 3; j++)
    xyz_a[j] = xyz[in_cell_id*3 + j];
  xyz_a[3] = 1;

  /* Initialize output */

  for (i = 0; i < 3; i++)
    xyz_b[i] = 0.;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 4; j++)
      xyz_b[i] += matrix[i][j]*xyz_a[j];

  /* Store updated cell center */

  for (j = 0; j < 3; j++)
    xyz[out_cell_id*3 + j] = xyz_b[j];

}

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a rotation to a given vector.
 *
 * parameters:
 *   matrix[3][4] --> matrix of the transformation in homogeneous coord.
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   x_in         --> X coord. of the incoming vector
 *   y_in         --> Y coord. of the incoming vector
 *   z_in         --> Z coord. of the incoming vector
 *   x_out        <-- pointer to the X coord. of the output
 *   y_out        <-- pointer to the Y coord. of the output
 *   z_out        <-- pointer to the Z coord. of the output
 *----------------------------------------------------------------------------*/

static void
_apply_vector_rotation_ni(cs_real_t   matrix[3][4],
                       cs_real_t   x_in,
                       cs_real_t   y_in,
                       cs_real_t   z_in,
                       cs_real_t   *x_out,
                       cs_real_t   *y_out,
                       cs_real_t   *z_out)
{
  *x_out = matrix[0][0] * x_in + matrix[0][1] * y_in + matrix[0][2] * z_in;
  *y_out = matrix[1][0] * x_in + matrix[1][1] * y_in + matrix[1][2] * z_in;
  *z_out = matrix[2][0] * x_in + matrix[2][1] * y_in + matrix[2][2] * z_in;
}

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a transformation to a given
 * interleaved vector.
 *
 * parameters:
 *   matrix[3][4] --> matrix of the transformation in homogeneous coord.
 *                    last line = [0; 0; 0; 1] (Not used here)
 *   xyz          <-> array of coordinates
 *----------------------------------------------------------------------------*/

static void
_apply_vector_rotation(cs_real_t    matrix[3][4],
                       cs_real_t   *xyz)
{
  cs_lnum_t  i;

  cs_real_t  t[3];
  for (i = 0; i < 3; i++)
    t[i] = xyz[i];

  /* Initialize output */

  for (i = 0; i < 3; i++)
    xyz[i] = matrix[i][0]*t[0] + matrix[i][1]*t[1] + matrix[i][2]*t[2];

}

/*----------------------------------------------------------------------------
 * Compute a matrix * tensor * Tmatrix product to apply a rotation to a
 * given tensor
 *
 * parameters:
 *   matrix[3][4]        --> transformation matrix in homogeneous coords.
 *                           last line = [0; 0; 0; 1] (Not used here)
 *   in11, in12, in13    --> incoming first line of the tensor
 *   in21, in22, in23    --> incoming second line of the tensor
 *   in31, in32, in33    --> incoming third line of the tensor
 *   out11, out12, out13 <-- pointer to the first line of the output
 *   out21, out22, out23 <-- pointer to the second line of the output
 *   out31, out32, out33 <-- pointer to the third line of the output
 *----------------------------------------------------------------------------*/

static void
_apply_tensor_rotation_ni(cs_real_t   matrix[3][4],
                          cs_real_t   in11,
                          cs_real_t   in12,
                          cs_real_t   in13,
                          cs_real_t   in21,
                          cs_real_t   in22,
                          cs_real_t   in23,
                          cs_real_t   in31,
                          cs_real_t   in32,
                          cs_real_t   in33,
                          cs_real_t   *out11,
                          cs_real_t   *out12,
                          cs_real_t   *out13,
                          cs_real_t   *out21,
                          cs_real_t   *out22,
                          cs_real_t   *out23,
                          cs_real_t   *out31,
                          cs_real_t   *out32,
                          cs_real_t   *out33)
{
  cs_lnum_t  i, j, k;
  cs_real_t  tensorA[3][3], tensorB[3][3];

  tensorA[0][0] = matrix[0][0]*in11 + matrix[0][1]*in12 + matrix[0][2]*in13;
  tensorA[0][1] = matrix[1][0]*in11 + matrix[1][1]*in12 + matrix[1][2]*in13;
  tensorA[0][2] = matrix[2][0]*in11 + matrix[2][1]*in12 + matrix[2][2]*in13;

  tensorA[1][0] = matrix[0][0]*in21 + matrix[0][1]*in22 + matrix[0][2]*in23;
  tensorA[1][1] = matrix[1][0]*in21 + matrix[1][1]*in22 + matrix[1][2]*in23;
  tensorA[1][2] = matrix[2][0]*in21 + matrix[2][1]*in22 + matrix[2][2]*in23;

  tensorA[2][0] = matrix[0][0]*in31 + matrix[0][1]*in32 + matrix[0][2]*in33;
  tensorA[2][1] = matrix[1][0]*in31 + matrix[1][1]*in32 + matrix[1][2]*in33;
  tensorA[2][2] = matrix[2][0]*in31 + matrix[2][1]*in32 + matrix[2][2]*in33;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      tensorB[i][j] = 0.;
      for (k = 0; k < 3; k++)
        tensorB[i][j] += matrix[i][k] * tensorA[k][j];
    }
  }

  *out11 = tensorB[0][0];
  *out22 = tensorB[1][1];
  *out33 = tensorB[2][2];

  if (out12 != NULL) {
    *out12 = tensorB[0][1];
    *out13 = tensorB[0][2];
    *out21 = tensorB[1][0];
    *out23 = tensorB[1][2];
    *out31 = tensorB[2][0];
    *out32 = tensorB[2][1];
  }

}

/*----------------------------------------------------------------------------
 * Compute a matrix * tensor * Tmatrix product to apply a rotation to a
 * given interleaved tensor
 *
 * parameters:
 *   matrix[3][4]        --> transformation matrix in homogeneous coords.
 *                           last line = [0; 0; 0; 1] (Not used here)
 *   tensor              <-> incoming 3x3 tensor
 *----------------------------------------------------------------------------*/

static void
_apply_tensor_rotation(cs_real_t   matrix[3][4],
                       cs_real_t   *tensor)
{
  cs_lnum_t  i, j, k, l;

  cs_real_t  t[3][3];

  for (k = 0; k < 3; k++) {
    for (j = 0; j < 3; j++) {
      t[k][j] = 0.;
      for (l = 0; l < 3; l++)
        t[k][j] += matrix[j][l] * tensor[k*3+l];
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      tensor[i*3+j] = 0.;
      for (k = 0; k < 3; k++)
        tensor[i*3+j] += matrix[i][k] * t[k][j];
    }
  }

}

/*----------------------------------------------------------------------------
 * Compute a matrix * tensor * Tmatrix product to apply a rotation to a
 * given symmetric interleaved tensor
 *
 * parameters:
 *   matrix[3][4]        --> transformation matrix in homogeneous coords.
 *                           last line = [0; 0; 0; 1] (Not used here)
 *   tensor              <-> incoming (6) symmetric tensor
 *----------------------------------------------------------------------------*/

static void
_apply_sym_tensor_rotation(cs_real_t   matrix[3][4],
                           cs_real_t   *tensor)
{
  cs_lnum_t  i, j, k, l;

  cs_real_t  t[3][3];
  cs_real_t  t0[3][3];

  t0[0][0] = tensor[0];
  t0[1][1] = tensor[1];
  t0[2][2] = tensor[2];
  t0[0][1] = tensor[3];
  t0[1][0] = tensor[3];
  t0[1][2] = tensor[4];
  t0[2][1] = tensor[4];
  t0[0][2] = tensor[5];
  t0[2][0] = tensor[5];

  for (k = 0; k < 3; k++) {
    for (j = 0; j < 3; j++) {
      t[k][j] = 0.;
      for (l = 0; l < 3; l++)
        t[k][j] += matrix[j][l] * t0[k][l];
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      t0[i][j] = 0.;
      for (k = 0; k < 3; k++)
        t0[i][j] += matrix[i][k] * t[k][j];
    }
  }

  tensor[0] = t0[0][0];
  tensor[1] = t0[1][1];
  tensor[2] = t0[2][2];
  tensor[3] = t0[0][1];
  tensor[3] = t0[1][0];
  tensor[4] = t0[2][1];
  tensor[5] = t0[2][0];

}

/*----------------------------------------------------------------------------
 * Compute the rotation of a third-order symmetric interleaved tensor
 * (18 components)
 * TENSOR_ijk = M_ip M_jq M_kr TENSOR_pqr
 *
 * WARNING: old fashion to store Rij (11, 22, 33, 12, 13, 23)
 *
 * parameters:
 *   matrix[3][4]        --> transformation matrix in homogeneous coords.
 *                           last line = [0; 0; 0; 1] (Not used here)
 *   tensor              <-> incoming 3x3x3 tensor
 *                           (in fact 3x6 due to symmetry)
 *----------------------------------------------------------------------------*/

static void
_apply_tensor3sym_rotation_old(cs_real_t   matrix[3][4],
                               cs_real_t   *tensor)
{
  cs_lnum_t  i, j, k, p, q, r;

  cs_real_t  t1[3][3][3], t2[3][3][3];

  for (p = 0; p < 3; p++) {
    for (q = 0; q < 3; q++) {
      for (k = 0; k < 3; k++) {
        t1[p][q][k] = 0.;
        for (r = 0; r < 3; r++)
          t1[p][q][k] += matrix[k][r] * tensor[3*_rij[p][q] + r];
      }
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        t2[i][j][k] = 0.;
        for (p = 0; p < 3; p++) {
          for (q = 0; q < 3; q++)
            t2[i][j][k] += matrix[i][p] * matrix[j][q] * t1[p][q][k];
        }
      }
    }
  }

  /* Output */

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++)
        tensor[3*_rij[i][j] + k] = t2[i][j][k];
    }
  }

}

/*----------------------------------------------------------------------------
 * Compute the rotation of a third-order symmetric interleaved tensor
 * (18 components)
 * TENSOR_ijk = M_ip M_jq M_kr TENSOR_pqr
 *
 * Warning: Rij stored as (11, 22, 33, 12, 23, 13)
 *
 * parameters:
 *   matrix[3][4]        --> transformation matrix in homogeneous coords.
 *                           last line = [0; 0; 0; 1] (Not used here)
 *   tensor              <-> incoming 3x3x3 tensor
 *                           (in fact 3x6 due to symmetry)
 *----------------------------------------------------------------------------*/

static void
_apply_tensor3sym_rotation(cs_real_t   matrix[3][4],
                           cs_real_t   *tensor)
{
  cs_lnum_t  i, j, k, p, q, r;

  cs_real_t  t1[3][3][3], t2[3][3][3];

  for (p = 0; p < 3; p++) {
    for (q = 0; q < 3; q++) {
      for (k = 0; k < 3; k++) {
        t1[p][q][k] = 0.;
        for (r = 0; r < 3; r++)
          t1[p][q][k] += matrix[k][r] * tensor[3*_symt[p][q] + r];
      }
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        t2[i][j][k] = 0.;
        for (p = 0; p < 3; p++) {
          for (q = 0; q < 3; q++)
            t2[i][j][k] += matrix[i][p] * matrix[j][q] * t1[p][q][k];
        }
      }
    }
  }

  /* Output */

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++)
        tensor[3*_symt[i][j] + k] = t2[i][j][k];
    }
  }

}

/*----------------------------------------------------------------------------
 * Test if a halo seems compatible with the main mesh's periodic
 * transformations.
 *
 * If a halo is not compatible, abort with an error message.
 *
 * parameters:
 *   halo --> pointer to halo structure
 *----------------------------------------------------------------------------*/

static void
_test_halo_compatibility(const cs_halo_t  *halo)
{
  assert(halo != NULL);

  if (cs_glob_mesh->n_transforms != halo->n_transforms)
    bft_error(__FILE__, __LINE__, 0,
              _("The %d periodic transformations of the halo do not comply\n"
                "with the main mesh transformations (numbering %d).\n"),
              halo->n_transforms, (int)(cs_glob_mesh->n_transforms));
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Rotate tensor values for periodic cells on extended halos.
 *
 * Fortran API:
 *
 * subroutine perrte
 * *****************
 *
 * double precision var11(ncelet) : <-> : component 11 of rank 2 tensor
 * double precision var12(ncelet) : <-> : component 12 of rank 2 tensor
 * double precision var13(ncelet) : <-> : component 13 of rank 2 tensor
 * double precision var21(ncelet) : <-> : component 21 of rank 2 tensor
 * double precision var22(ncelet) : <-> : component 22 of rank 2 tensor
 * double precision var23(ncelet) : <-> : component 23 of rank 2 tensor
 * double precision var31(ncelet) : <-> : component 31 of rank 2 tensor
 * double precision var32(ncelet) : <-> : component 32 of rank 2 tensor
 * double precision var33(ncelet) : <-> : component 33 of rank 2 tensor
 *----------------------------------------------------------------------------*/

void
CS_PROCF (perrte, PERRTE) (cs_real_t  var11[],
                           cs_real_t  var12[],
                           cs_real_t  var13[],
                           cs_real_t  var21[],
                           cs_real_t  var22[],
                           cs_real_t  var23[],
                           cs_real_t  var31[],
                           cs_real_t  var32[],
                           cs_real_t  var33[])
{
  const cs_halo_t *halo = cs_glob_mesh->halo;

  if (halo == NULL)
    return;

  cs_halo_perio_sync_var_tens_ni(halo,
                                 CS_HALO_EXTENDED,
                                 var11, var12, var13,
                                 var21, var22, var23,
                                 var31, var32, var33);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Apply transformation on coordinates.
 *
 * parameters:
 *   halo      <-> halo associated with coordinates to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   coords    <-- coordinates on which transformation have to be done.
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_coords(const cs_halo_t  *halo,
                          cs_halo_type_t    sync_mode,
                          cs_real_t        *coords)
{
  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

  cs_real_t  matrix[3][4];

  const fvm_periodicity_t  *periodicity = cs_glob_mesh->periodicity;
  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts = halo->n_local_elts;

  if (sync_mode == CS_HALO_N_TYPES)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  /* Compute the new cell centers through periodicity */

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    fvm_periodicity_get_matrix(periodicity, t_id, matrix);

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      /* apply transformation for standard halo */

      start_std = halo->perio_lst[shift + 4*rank_id];
      end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

      for (i = start_std; i < end_std; i++)
        _apply_vector_transfo(matrix, n_elts+i, n_elts+i, coords);

      /* apply transformation for extended halo */

      if (sync_mode == CS_HALO_EXTENDED) {

        start_ext = halo->perio_lst[shift + 4*rank_id + 2];
        end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

        for (i = start_ext; i < end_ext; i++)
          _apply_vector_transfo(matrix, n_elts+i, n_elts+i, coords);

      } /* End if extended halo */

    } /* End of loop on ranks */

  } /* End of loop on transformation */
}

/*----------------------------------------------------------------------------
 * Synchronize values for a real vector (interleaved) between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   var       <-> vector to update
 *   incvar    <-- specifies the increment for the elements of var
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_vect(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            cs_real_t         var[],
                            int               incvar)
{
  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

  cs_real_t matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t  *periodicity = cs_glob_mesh->periodicity;
  const int  have_rotation = cs_glob_mesh->have_rotation_perio;

  if (sync_mode == CS_HALO_N_TYPES || have_rotation == 0)
    return;

  assert(halo != NULL);
  assert(incvar == 3);

  _test_halo_compatibility(halo);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = n_elts + halo->perio_lst[shift + 4*rank_id];
        end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

        for (i = start_std; i < end_std; i++)
          _apply_vector_rotation(matrix, var + i*incvar);

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = n_elts + halo->perio_lst[shift + 4*rank_id + 2];
          end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

          for (i = start_ext; i < end_ext; i++)
            _apply_vector_rotation(matrix, var + i*incvar);

        }

      } /* End of loop on ranks */

    } /* End of the treatment of rotation */

  } /* End of loop on transformations */
}

/*----------------------------------------------------------------------------
 * Synchronize values for a real vector between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   var_x     <-> component of the vector to update
 *   var_y     <-> component of the vector to update
 *   var_z     <-> component of the vector to update
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_vect_ni(const cs_halo_t     *halo,
                               cs_halo_type_t       sync_mode,
                               cs_real_t            var_x[],
                               cs_real_t            var_y[],
                               cs_real_t            var_z[])
{
  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;
  cs_real_t  x_in, y_in, z_in;

  cs_real_t matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = cs_glob_mesh->periodicity;
  const int  have_rotation = cs_glob_mesh->have_rotation_perio;

  if (sync_mode == CS_HALO_N_TYPES || have_rotation == 0)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

        for (i = start_std; i < end_std; i++) {

          x_in = var_x[n_elts + i];
          y_in = var_y[n_elts + i];
          z_in = var_z[n_elts + i];

          _apply_vector_rotation_ni(matrix,
                                    x_in,
                                    y_in,
                                    z_in,
                                    &var_x[n_elts+i],
                                    &var_y[n_elts+i],
                                    &var_z[n_elts+i]);
        }

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

          for (i = start_ext; i < end_ext; i++) {

            x_in = var_x[n_elts + i];
            y_in = var_y[n_elts + i];
            z_in = var_z[n_elts + i];

            _apply_vector_rotation_ni(matrix,
                                      x_in,
                                      y_in,
                                      z_in,
                                      &var_x[n_elts+i],
                                      &var_y[n_elts+i],
                                      &var_z[n_elts+i]);

          }

        }

      } /* End of loop on ranks */

    } /* End of the treatment of rotation */

  } /* End of loop on transformations */
}

/*----------------------------------------------------------------------------
 * Synchronize values for a real tensor between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   var11     <-> component of the tensor to update
 *   var12     <-> component of the tensor to update
 *   var13     <-> component of the tensor to update
 *   var21     <-> component of the tensor to update
 *   var22     <-> component of the tensor to update
 *   var23     <-> component of the tensor to update
 *   var31     <-> component of the tensor to update
 *   var32     <-> component of the tensor to update
 *   var33     <-> component of the tensor to update
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_tens_ni(const cs_halo_t  *halo,
                               cs_halo_type_t    sync_mode,
                               cs_real_t         var11[],
                               cs_real_t         var12[],
                               cs_real_t         var13[],
                               cs_real_t         var21[],
                               cs_real_t         var22[],
                               cs_real_t         var23[],
                               cs_real_t         var31[],
                               cs_real_t         var32[],
                               cs_real_t         var33[])
{
  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;
  cs_real_t  v11, v12, v13, v21, v22, v23, v31, v32, v33;

  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = cs_glob_mesh->periodicity;
  const int  have_rotation = cs_glob_mesh->have_rotation_perio;

  if (sync_mode == CS_HALO_N_TYPES || have_rotation == 0)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

        for (i = start_std; i < end_std; i++) {

          v11 = var11[n_elts + i];
          v12 = var12[n_elts + i];
          v13 = var13[n_elts + i];
          v21 = var21[n_elts + i];
          v22 = var22[n_elts + i];
          v23 = var23[n_elts + i];
          v31 = var31[n_elts + i];
          v32 = var32[n_elts + i];
          v33 = var33[n_elts + i];

          _apply_tensor_rotation_ni(matrix,
                                    v11, v12, v13, v21, v22, v23,
                                    v31, v32, v33,
                                    &var11[n_elts + i], &var12[n_elts + i],
                                    &var13[n_elts + i], &var21[n_elts + i],
                                    &var22[n_elts + i], &var23[n_elts + i],
                                    &var31[n_elts + i], &var32[n_elts + i],
                                    &var33[n_elts + i]);

        }

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

          for (i = start_ext; i < end_ext; i++) {

            v11 = var11[n_elts + i];
            v12 = var12[n_elts + i];
            v13 = var13[n_elts + i];
            v21 = var21[n_elts + i];
            v22 = var22[n_elts + i];
            v23 = var23[n_elts + i];
            v31 = var31[n_elts + i];
            v32 = var32[n_elts + i];
            v33 = var33[n_elts + i];

            _apply_tensor_rotation_ni(matrix,
                                      v11, v12, v13, v21, v22, v23,
                                      v31, v32, v33,
                                      &var11[n_elts + i], &var12[n_elts + i],
                                      &var13[n_elts + i], &var21[n_elts + i],
                                      &var22[n_elts + i], &var23[n_elts + i],
                                      &var31[n_elts + i], &var32[n_elts + i],
                                      &var33[n_elts + i]);

          }

        } /* End of the treatment of rotation */

      } /* End if halo is extended */

    } /* End of loop on ranks */

  } /* End of loop on transformations for the local rank */
}

/*----------------------------------------------------------------------------
 * Synchronize values for a real tensor (interleaved) between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   var       <-> tensor to update
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_tens(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            cs_real_t         var[])
{
  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = cs_glob_mesh->periodicity;
  const int  have_rotation = cs_glob_mesh->have_rotation_perio;

  if (sync_mode == CS_HALO_N_TYPES || have_rotation == 0)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

        for (i = start_std; i < end_std; i++)
          _apply_tensor_rotation(matrix, var + 9*(n_elts+i));

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

          for (i = start_ext; i < end_ext; i++)
            _apply_tensor_rotation(matrix, var + 9*(n_elts+i));

        } /* End of the treatment of rotation */

      } /* End if halo is extended */

    } /* End of loop on ranks */

  } /* End of loop on transformations for the local rank */
}

/*----------------------------------------------------------------------------
 * Synchronize values for a real tensor (symmetric interleaved) between
 * periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   var       <-> symmetric tensor to update (6 values)
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_sym_tens(const cs_halo_t  *halo,
                                cs_halo_type_t    sync_mode,
                                cs_real_t         var[])
{
  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = cs_glob_mesh->periodicity;
  const int  have_rotation = cs_glob_mesh->have_rotation_perio;

  if (sync_mode == CS_HALO_N_TYPES || have_rotation == 0)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

        for (i = start_std; i < end_std; i++)
          _apply_sym_tensor_rotation(matrix, var + 6*(n_elts+i));

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

          for (i = start_ext; i < end_ext; i++)
            _apply_sym_tensor_rotation(matrix, var + 6*(n_elts+i));

        } /* End of the treatment of rotation */

      } /* End if halo is extended */

    } /* End of loop on ranks */

  } /* End of loop on transformations for the local rank */
}

/*----------------------------------------------------------------------------
 * Synchronize values for a real diagonal tensor between periodic cells.
 *
 * We only know the diagonal of the tensor.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   var11     <-> component of the tensor to update
 *   var22     <-> component of the tensor to update
 *   var33     <-> component of the tensor to update
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_diag_ni(const cs_halo_t  *halo,
                               cs_halo_type_t    sync_mode,
                               cs_real_t         var11[],
                               cs_real_t         var22[],
                               cs_real_t         var33[])
{
  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;
  cs_real_t  v11, v22, v33;
  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = cs_glob_mesh->periodicity;
  const int  have_rotation = cs_glob_mesh->have_rotation_perio;

  if (sync_mode == CS_HALO_N_TYPES || have_rotation == 0)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

        for (i = start_std; i < end_std; i++) {

          v11 = var11[n_elts + i];
          v22 = var22[n_elts + i];
          v33 = var33[n_elts + i];

          _apply_tensor_rotation_ni(matrix, v11, 0, 0,
                                    0, v22, 0, 0, 0, v33,
                                    &var11[n_elts + i], NULL, NULL,
                                    NULL, &var22[n_elts + i], NULL,
                                    NULL, NULL, &var33[n_elts + i]);

        }

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

          for (i = start_ext; i < end_ext; i++) {

            v11 = var11[n_elts + i];
            v22 = var22[n_elts + i];
            v33 = var33[n_elts + i];

            _apply_tensor_rotation_ni(matrix, v11, 0, 0,
                                      0, v22, 0, 0, 0, v33,
                                      &var11[n_elts + i], NULL, NULL,
                                      NULL, &var22[n_elts + i], NULL,
                                      NULL, NULL, &var33[n_elts + i]);

          }

        } /* End if halo is extended */

      } /* End of loop on ranks */

    } /* End of the treatment of rotation */

  } /* End of loop on transformations */
}

/*----------------------------------------------------------------------------
 * Apply rotation on the gradient of Reynolds stress tensor
 *
 * parameters:
 *   drdxyz     <-> gradient on the variable (size: 3*6*n_ghost_cells)
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_rotate_rij(cs_real_t  *drdxyz)
{
  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;
  fvm_periodicity_type_t  perio_type;

  cs_real_t  matrix[3][4];

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_halo_t  *halo = mesh->halo;

  const int  n_transforms = mesh->n_transforms;
  const cs_halo_type_t  halo_type = mesh->halo_type;
  const fvm_periodicity_t  *periodicity = mesh->periodicity;

  if (mesh->halo_type == CS_HALO_N_TYPES || halo == NULL)
    return;

  assert(halo != NULL);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

        for (i = start_std; i < end_std; i++)
          _apply_tensor3sym_rotation_old(matrix, drdxyz + 18*i);

        if (halo_type == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

          for (i = start_ext; i < end_ext; i++)
            _apply_tensor3sym_rotation_old(matrix, drdxyz + 18*i);

        } /* End if an extended halo exists */

      } /* End of loop on ranks */

    } /* If the transformation is a rotation */

  } /* End of loop on transformations */
}

/*----------------------------------------------------------------------------
 * Synchronize values for a real gradient of a tensor (symmetric interleaved)
 * between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   var       <-> symmetric tensor to update (6 values)
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_sym_tens_grad(const cs_halo_t  *halo,
                                     cs_halo_type_t    sync_mode,
                                     cs_real_t         var[])
{
  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t  perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = cs_glob_mesh->periodicity;
  const int  have_rotation = cs_glob_mesh->have_rotation_perio;

  if (sync_mode == CS_HALO_N_TYPES || have_rotation == 0)
    return;

  assert(halo != NULL);

  _test_halo_compatibility(halo);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std =n_elts +  halo->perio_lst[shift + 4*rank_id];
        end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

        for (i = start_std; i < end_std; i++)
          _apply_tensor3sym_rotation(matrix, var + 18*i);

        if (sync_mode == CS_HALO_EXTENDED) {

          start_ext = n_elts + halo->perio_lst[shift + 4*rank_id + 2];
          end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

          for (i = start_ext; i < end_ext; i++)
            _apply_tensor3sym_rotation(matrix, var + 18*i);

        } /* End of the treatment of rotation */

      } /* End if halo is extended */

    } /* End of loop on ranks */

  } /* End of loop on transformations for the local rank */
}


/*----------------------------------------------------------------------------*/

END_C_DECLS
