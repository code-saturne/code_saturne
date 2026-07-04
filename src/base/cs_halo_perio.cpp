/*============================================================================
 * Functions handling periodicity.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm/fvm_periodicity.h"

#include "base/cs_base.h"
#include "base/cs_halo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_halo_perio.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Table giving the Reynolds stress component for [i][j] */

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
 *   cell_id      --> cell_id of the parent cell.
 *   xyz          <-> array of coordinates
 *----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
_apply_vector_transfo(const cs_real_t    matrix[3][4],
                      cs_lnum_t          cell_id,
                      cs_real_t         *xyz)
{
  cs_real_t  xyz_a[4];
  cs_real_t  xyz_b[3];

  /* Copy to homogeneous coordinates before transformation */

  for (cs_lnum_t j = 0; j < 3; j++)
    xyz_a[j] = xyz[cell_id*3 + j];
  xyz_a[3] = 1;

  /* Initialize output */

  for (cs_lnum_t i = 0; i < 3; i++)
    xyz_b[i] = 0.;

  for (cs_lnum_t i = 0; i < 3; i++) {
    for (cs_lnum_t j = 0; j < 4; j++)
      xyz_b[i] += matrix[i][j]*xyz_a[j];
  }

  /* Store updated cell center */

  for (cs_lnum_t j = 0; j < 3; j++)
    xyz[cell_id*3 + j] = xyz_b[j];
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

template <typename T>
static void
_apply_vector_rotation(cs_real_t    matrix[3][4],
                       T           *xyz)
{
  cs_real_t  t[3];
  for (cs_lnum_t i = 0; i < 3; i++)
    t[i] = xyz[i];

  /* Initialize output */

  for (cs_lnum_t i = 0; i < 3; i++)
    xyz[i] = matrix[i][0]*t[0] + matrix[i][1]*t[1] + matrix[i][2]*t[2];
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

template <typename T>
static void
_apply_tensor_rotation(cs_real_t   matrix[3][4],
                       T          *tensor)
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

template <typename T>
static void
_apply_sym_tensor_rotation(cs_real_t   matrix[3][4],
                           T          *tensor)
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
 * Warning: Rij stored as (11, 22, 33, 12, 23, 13)
 *
 * parameters:
 *   matrix[3][4]        --> transformation matrix in homogeneous coords.
 *                           last line = [0; 0; 0; 1] (Not used here)
 *   tensor              <-> incoming 3x3x3 tensor
 *                           (in fact 3x6 due to symmetry)
 *----------------------------------------------------------------------------*/

template <typename T>
static void
_apply_tensor3sym_rotation(cs_real_t   matrix[3][4],
                           T          *tensor)
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
  if (   halo == nullptr
      || sync_mode == CS_HALO_N_TYPES)
    return;

  const fvm_periodicity_t  *periodicity = halo->periodicity;
  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts = halo->n_local_elts;
  const int n_h = (sync_mode == CS_HALO_EXTENDED) ? 2 : 1;

  assert(halo != nullptr);

  /* Compute the new cell centers through periodicity */

  for (int t_id = 0; t_id < n_transforms; t_id++) {

    cs_lnum_t shift = 4 * halo->n_c_domains * t_id;

    cs_real_t  matrix[3][4];
    fvm_periodicity_get_matrix(periodicity, t_id, matrix);

    for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
      for (int h_idx = 0; h_idx < n_h; h_idx++) {

        int h_shift = h_idx*2;
        cs_lnum_t start_std = halo->perio_lst[shift + 4*rank_id + h_shift];
        cs_lnum_t end_std = start_std
          + halo->perio_lst[shift + 4*rank_id + h_shift + 1];

        for (cs_lnum_t i = start_std; i < end_std; i++)
          _apply_vector_transfo(matrix, n_elts+i, coords);

      }  /* End of loop on halo type */
    } /* End of loop on ranks */

  } /* End of loop on transformation */
}

/*----------------------------------------------------------------------------
 * Apply transformation on coordinates.
 *
 * parameters:
 *   halo      <-> halo associated with coordinates to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   coords    <-- coordinates on which transformation have to be done.
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_coords(const cs_halo_t      *halo,
                          cs_dispatch_context  &ctx,
                          cs_halo_type_t        sync_mode,
                          cs_real_t            *coords)
{
  if (   halo == nullptr
      || sync_mode == CS_HALO_N_TYPES)
    return;

  const fvm_periodicity_t  *periodicity = halo->periodicity;
  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts = halo->n_local_elts;

  assert(halo != nullptr);

  /* Compute the new cell centers through periodicity */

  for (int t_id = 0; t_id < n_transforms; t_id++) {

    cs_lnum_t shift = 4 * halo->n_c_domains * t_id;

    // Allow for grouping calls for n_bounds_max non-empty exchanges.
    // We want to use lambda-capture to pass these bounds to the device,
    // Without requiring an additional CPU/GPU memory exchange, and avoiding
    // the overhead of an additional exchange, while still minimizing the
    // number of kernel launches, as they have their own overhead.

    constexpr int n_bounds_max = 8;
    const int n_h = (sync_mode == CS_HALO_EXTENDED) ? 2 : 1;
    cs_lnum_t  bounds[2][n_bounds_max];

    cs_real_t  matrix[3][4];
    fvm_periodicity_get_matrix(periodicity, t_id, matrix);

    int idx = 0;
    int idx_e = halo->n_c_domains*n_h;

    while (idx < idx_e) {

      int n_bounds = 0;
      cs_lnum_t n = 0;

      while (idx < idx_e && n_bounds < n_bounds_max) {
        int rank_id = idx/n_h;
        int h_shift = (idx%n_h)*2;

        cs_lnum_t start_std = halo->perio_lst[shift + 4*rank_id + h_shift];
        cs_lnum_t end_std = start_std
          + halo->perio_lst[shift + 4*rank_id + h_shift + 1];
        cs_lnum_t n0 = end_std - start_std;
        if (n0 > 0) {
          bounds[0][n_bounds] = start_std;
          bounds[1][n_bounds] = n0;
          n += n0;
          n_bounds++;
        }
        idx += 1;
      }

      /* apply transformation for standard halo */
      ctx.parallel_for(n, [=] CS_F_HOST_DEVICE (cs_lnum_t idxb) {
        // Determine matching ghost element id from index and shifs.
        cs_lnum_t i = idxb;
        int k = 0;
        while (i >= bounds[1][k]) {
          i -= bounds[1][k];
          k++;
        }
        assert(k < n_bounds);
        i += bounds[0][k];

        // Now apply transformation
        _apply_vector_transfo(matrix, n_elts+i, coords);
      });

      n_bounds = 0;  // reset for next block

    } /* End of loop on ranks and halo types */

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

template <typename T>
void
cs_halo_perio_sync_var_vect(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            T                 var[],
                            int               incvar)
{
  if (   halo == nullptr
      || sync_mode == CS_HALO_N_TYPES)
    return;

  if (halo->n_rotations == 0)
    return;

  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

  cs_real_t matrix[3][4];

  fvm_periodicity_type_t perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t  *periodicity = halo->periodicity;

  assert(halo != nullptr);
  assert(incvar == 3);

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

// Force instanciation

template void
cs_halo_perio_sync_var_vect(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            double            var[],
                            int               incvar);

template void
cs_halo_perio_sync_var_vect(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            float             var[],
                            int               incvar);

/*----------------------------------------------------------------------------
 * Synchronize values for a real tensor (interleaved) between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   var       <-> tensor to update
 *----------------------------------------------------------------------------*/

template <typename T>
void
cs_halo_perio_sync_var_tens(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            T                 var[])
{
  if (   halo == nullptr
      || sync_mode == CS_HALO_N_TYPES)
    return;

  if (halo->n_rotations == 0)
    return;

  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = halo->periodicity;

  assert(halo != nullptr);

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

// Force instanciation

template void
cs_halo_perio_sync_var_tens(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            double            var[]);

template void
cs_halo_perio_sync_var_tens(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            float             var[]);

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
  if (   halo == nullptr
      || sync_mode == CS_HALO_N_TYPES)
    return;

  if (halo->n_rotations == 0)
    return;

  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = halo->periodicity;

  assert(halo != nullptr);

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
 * Synchronize values for a real gradient of a tensor (symmetric interleaved)
 * between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode <-- kind of halo treatment (standard or extended)
 *   var       <-> symmetric tensor to update (6 values)
 *----------------------------------------------------------------------------*/

template <typename T>
void
cs_halo_perio_sync_var_sym_tens_grad(const cs_halo_t  *halo,
                                     cs_halo_type_t    sync_mode,
                                     T                 var[])
{
  if (   halo == nullptr
      || sync_mode == CS_HALO_N_TYPES)
    return;

  if (halo->n_rotations == 0)
    return;

  int  rank_id, t_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

  cs_real_t  matrix[3][4];

  fvm_periodicity_type_t perio_type = FVM_PERIODICITY_NULL;

  const int  n_transforms = halo->n_transforms;
  const cs_lnum_t  n_elts   = halo->n_local_elts;
  const fvm_periodicity_t *periodicity = halo->periodicity;

  assert(halo != nullptr);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    shift = 4 * halo->n_c_domains * t_id;

    perio_type = fvm_periodicity_get_type(periodicity, t_id);

    if (perio_type >= FVM_PERIODICITY_ROTATION) {

      fvm_periodicity_get_matrix(periodicity, t_id, matrix);

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = n_elts +  halo->perio_lst[shift + 4*rank_id];
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

// Force instanciation

template void
cs_halo_perio_sync_var_sym_tens_grad(const cs_halo_t  *halo,
                                     cs_halo_type_t    sync_mode,
                                     double            var[]);

template void
cs_halo_perio_sync_var_sym_tens_grad(const cs_halo_t  *halo,
                                     cs_halo_type_t    sync_mode,
                                     float             var[]);

/*----------------------------------------------------------------------------*/
