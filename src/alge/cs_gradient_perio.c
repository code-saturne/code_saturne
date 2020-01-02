/*============================================================================
 * Gradient reconstruction.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_blas.h"
#include "cs_field.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_ext_neighborhood.h"
#include "cs_mesh_quantities.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gradient.h"
#include "cs_gradient_perio.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_gradient_perio.c
        Specific functions for gradient reconstruction of the Reynolds stress
        tensor in the case of periodicity of rotation.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Specific handling for Reynolds stresses with rotational periodicity */

static cs_real_t  *_drdxyz = NULL, *_wdrdxy = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update drdxyz and wdrdxy for periodic ghost cells.
 *
 * Called by PERMAS.
 *
 * parameters:
 *   h_cell_id --> cell id in halo
 *   cell_id   --> cell id
 *   rom       --> density array
 *   call_id   --> first or second call
 *   drdxyz    <-> Gradient on components of Rij (Reynolds stress tensor)
 *   wdrdxy    <-> associated working array.
 *----------------------------------------------------------------------------*/

static void
_update_drdxyz(cs_lnum_t         h_cell_id,
               cs_lnum_t         cell_id,
               const cs_real_t   rom[],
               int               call_id,
               cs_real_t         drdxyz[],
               cs_real_t         wdrdxy[])
{
  cs_lnum_t  i, j, id;

  if (call_id == 1) { /* First call */

    for (i = 0; i < 2*3; i++) {
      for (j = 0; j < 3; j++) {
        id = j + 3*i + 3*6*h_cell_id;
        wdrdxy[id] = drdxyz[id];
        drdxyz[id] *= rom[cell_id];
      }
    }

  }
  else if (call_id == 2) { /* Second call */

    for (i = 0; i < 2*3; i++) {
      for (j = 0; j < 3; j++) {
        id = j + 3*i + 3*6*h_cell_id;
        drdxyz[id] = wdrdxy[id];
      }
    }

  } /* End if second call */
}

/*----------------------------------------------------------------------------
 * Save ghost cell values of an initial Rij component gradient.
 *
 * Note that ghost cell values of comp_grad are synchronized by this function.
 *
 * parameters:
 *   comp_id    <-- Rij component id
 *   dxyz       <-> gradient on the variable (drdxy)
 *   comp_grad  <-> component gradient
 *----------------------------------------------------------------------------*/

static void
_save_halo_perio_rij(int           comp_id,
                     cs_real_t    *dxyz,
                     cs_real_3_t  *comp_grad)
{
  int  t_id, rank_id;
  cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_halo_t  *halo = mesh->halo;

  const int  strid_v = 3*comp_id;
  const int  strid_e = 3*6;

  const cs_lnum_t  n_cells = mesh->n_cells;
  const int  n_transforms = mesh->n_transforms;
  const fvm_periodicity_t  *periodicity = mesh->periodicity;

  cs_halo_sync_var_strided(mesh->halo,
                           mesh->halo_type,
                           (cs_real_t *)comp_grad,
                           3);

  for (t_id = 0; t_id < n_transforms; t_id++) {

    if (   fvm_periodicity_get_type(periodicity, t_id)
        >= FVM_PERIODICITY_ROTATION) {

      shift = 4 * halo->n_c_domains * t_id;

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

        for (i = start_std; i < end_std; i++) {
          dxyz[0 + strid_v + strid_e*i] = comp_grad[n_cells + i][0];
          dxyz[1 + strid_v + strid_e*i] = comp_grad[n_cells + i][1];
          dxyz[2 + strid_v + strid_e*i] = comp_grad[n_cells + i][2];
        }

        if (mesh->halo_type == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

          for (i = start_ext; i < end_ext; i++) {
            dxyz[0 + strid_v + strid_e*i] = comp_grad[n_cells + i][0];
            dxyz[1 + strid_v + strid_e*i] = comp_grad[n_cells + i][1];
            dxyz[2 + strid_v + strid_e*i] = comp_grad[n_cells + i][2];
          }

        } /* End if extended halo */

      } /* End of loop on ranks */

    } /* End of test on rotation */

  } /* End of loop on transformations */

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Periodicity management for INIMAS
 *
 * If INIMAS is called by NAVSTO :
 *    We assume that gradient on ghost cells given by a rotation is known
 *    and is equal to the velocity one for the previous time step.
 * If INIMAS is called by DIVRIJ
 *    We assume that (more justifiable than in previous case) gradient on
 *    ghost cells given by rotation is equal to Rij gradient for the previous
 *    time step.
 *
 * Fortran Interface:
 *
 * subroutine permas
 * *****************
 *
 * integer          iappel      :  -> : indicateur d'appel dans inimas
 *                                          = 1 si appel au debut
 *                                          = 2 si appel a la fin
 * double precision rom(ncelet) :  -> : masse volumique aux cellules
 *
 * Size of DRDXYZ and WDRDXY = n_ghost_cells*6*3
 *----------------------------------------------------------------------------*/

void
CS_PROCF (permas, PERMAS)(const cs_int_t  *iappel,
                          cs_real_t        rom[])
{
  int  rank_id, t_id;
  cs_lnum_t  i, cell_id, shift;
  cs_lnum_t  start_std, end_std, start_ext, end_ext;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_halo_t  *halo = mesh->halo;
  cs_halo_type_t  halo_type = mesh->halo_type;
  const fvm_periodicity_t  *periodicity = mesh->periodicity;

  if (halo_type == CS_HALO_N_TYPES)
    return;

  if (*iappel == 1)
    cs_halo_sync_var(mesh->halo, mesh->halo_type, rom);

  for (t_id = 0; t_id < mesh->n_transforms; t_id++) {

    if (   fvm_periodicity_get_type(periodicity, t_id)
        >= FVM_PERIODICITY_ROTATION) {

      shift = 4 * halo->n_c_domains * t_id;

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start_std = halo->perio_lst[shift + 4*rank_id];
        end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

        for (i = start_std; i < end_std; i++) {

          cell_id = mesh->n_cells + i;

          _update_drdxyz(i, cell_id, rom, *iappel, _drdxyz, _wdrdxy);

        } /* End of loop on halo elements */

        if (halo_type == CS_HALO_EXTENDED) {

          start_ext = halo->perio_lst[shift + 4*rank_id + 2];
          end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

          for (i = start_ext; i < end_ext; i++) {

            cell_id = mesh->n_cells + i;

            _update_drdxyz(i, cell_id, rom, *iappel, _drdxyz, _wdrdxy);

          } /* End of loop on halo elements */

        } /* End if extended halo */

      } /* End of loop on ranks */

    } /* End of test on rotation */

  } /* End of loop on transformation */
}

/*----------------------------------------------------------------------------
 * Preparation rotation periodicity for Reynolds stresses.
 *
 * Compute an estimation of the velocity gradient.
 * The gradient is not reconstructed (as we would need a gradient,
 * thus periodicity). I seems possible to use a least-squares gradient.
 *
 * The gradient is then saved in an array representing the ghost cells, then
 * rotated where necessary to be ready for use (cf. cs_gradient_perio_init_rij).
 *
 * Compute cell gradient of vector field.
 *----------------------------------------------------------------------------*/

void CS_PROCF (perinr, PERINR)
(
 const cs_int_t   *const imrgra,  /* <-- gradient computation mode            */
 const cs_int_t   *const iwarnp,  /* <-- verbosity level                      */
 const cs_real_t  *const epsrgp,  /* <-- precision for iterative gradient
                                         calculation                          */
 const cs_real_t  *const extrap   /* <-- extrapolate gradient at boundary     */
)
{
  int i;
  cs_real_3_t *grad;
  const char *r_name[] = {"r11", "r22", "r33", "r12", "r13", "r23"};
  cs_field_t  *f = NULL;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  cs_gradient_type_by_imrgra(*imrgra,
                             &gradient_type,
                             &halo_type);

  BFT_MALLOC(grad, mesh->n_cells_with_ghosts, cs_real_3_t);

  /* Loop on R components */

  for (i = 0; i < 6; i++) {

    int tr_dim = 0;

    f = cs_field_by_name_try(r_name[i]);

    if (f == NULL) {
      BFT_FREE(grad);
      return;
    }

    /* We do not reconstruct and do not limit as we do not know the gradient of
       neighbors (we are trying to compute it here). Starting from the second
       time step, we have previous values. */

    /* Caution: after the gradient calculation here, the rotation halo contains:
       - nothing at the first relative time step
       - otherwise, the gradient computed at the previous time step
       We will be careful not to use the halo in the general case. */

    cs_gradient_perio_init_rij(f, &tr_dim, grad);

    assert(tr_dim == 2);

    cs_gradient_scalar(f->name,
                       gradient_type,
                       halo_type,
                       0,               /* inc */
                       true,            /* recompute_cocg */
                       1,               /* n_r_sweeps */
                       tr_dim,
                       false,           /* hyd_p_flag, */
                       1,               /* w_stride */
                       *iwarnp,
                       -1,              /* clip_mode */
                       *epsrgp,
                       *extrap,
                       1.5,             /* clip_coeff (ignored here) */
                       NULL,            /* f_ext */
                       f->bc_coeffs->a,
                       f->bc_coeffs->b,
                       f->val,
                       NULL,            /* weight_var */
                       NULL,            /* internal coupling */
                       grad);

    if (_drdxyz == NULL) {
      BFT_MALLOC(_drdxyz, mesh->n_ghost_cells * 3 * 6, cs_real_t);
      BFT_MALLOC(_wdrdxy, mesh->n_ghost_cells * 3 * 6, cs_real_t);
    }

    _save_halo_perio_rij(i, _drdxyz, grad);

  }

  cs_halo_perio_rotate_rij(_drdxyz);

  BFT_FREE(grad);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize gradient rotational periodicity computation API.
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_perio_initialize(void)
{
  BFT_FREE(_drdxyz);
  BFT_FREE(_wdrdxy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize gradient rotational periodicity computation API.
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_perio_finalize(void)
{
  BFT_FREE(_drdxyz);
  BFT_FREE(_wdrdxy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update gradient rotational periodicity computation API in case of
 *         mesh modification.
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_perio_update_mesh(void)
{
  BFT_FREE(_drdxyz);
  BFT_FREE(_wdrdxy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize ghost cell values for Reynolds stress tensor gradient.
 *
 * We retrieve the gradient given by perinr (phyvar) for the Reynolds
 * stress tensor in a buffer on ghost cells. then we define
 * dpdx, dpdy and dpdz gradient (1 -> n_cells_with_ghosts).
 *
 * We can't implicitly take into account rotation of a gradient of a tensor
 * variable because we have to know all components.
 *
 * We set idimtr to 2 for the Reynolds stress tensor.
 *
 * We assume that is correct to treat periodicities implicitly for the other
 * variables when reconstructing gradients.
 *
 * \param[in]       f        pointer to field
 * \param[out]      tr_dim   2 for tensor (Rij) in case of rotation, 0 otherwise
 * \param[in, out]  grad     gradient of field
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_perio_init_rij(const cs_field_t  *f,
                           int               *tr_dim,
                           cs_real_3_t        grad[])
{
  cs_lnum_t  d_var = -1;
  cs_mesh_t  *mesh = cs_glob_mesh;

  if (f->name[0] == 'r') {
    if (strlen(f->name) == 3) {
      if (f->name[1] == '1') {
        if (f->name[2] == '1')
          d_var = 0;
        else if (f->name[2] == '2')
          d_var = 3;
        else if (f->name[2] == '3')
          d_var = 4;
      }
      else if (f->name[1] == '2') {
        if (f->name[2] == '2')
          d_var = 1;
        else if (f->name[2] == '3')
          d_var = 5;
      }
      else if (f->name[1] == '3' && f->name[2] == '3')
        d_var = 2;
    }
  }

  if (mesh->halo == NULL || d_var < 0) {
    *tr_dim = 0;
    return;
  }
  else
    *tr_dim = 2;

  /*
    When there is periodicity of rotation :
      - Retrieve gradient values for ghost cells and for the previous
        time step without reconstruction
      - Ghost cells without rotation keep their value.
  */

  if (_drdxyz != NULL) {

    int  rank_id, t_id;
    cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

    cs_halo_t  *halo = mesh->halo;

    const cs_lnum_t  n_cells   = mesh->n_cells;
    const cs_lnum_t  n_transforms = mesh->n_transforms;
    const fvm_periodicity_t  *periodicity = mesh->periodicity;

    for (t_id = 0; t_id < n_transforms; t_id++) {

      if (   fvm_periodicity_get_type(periodicity, t_id)
          >= FVM_PERIODICITY_ROTATION) {

        shift = 4 * halo->n_c_domains * t_id;

        for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

          start_std = halo->perio_lst[shift + 4*rank_id];
          end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

          for (i = start_std; i < end_std; i++) {
            grad[n_cells + i][0] = _drdxyz[0 + d_var*3 + 18*i];
            grad[n_cells + i][1] = _drdxyz[1 + d_var*3 + 18*i];
            grad[n_cells + i][2] = _drdxyz[2 + d_var*3 + 18*i];
          }

          if (mesh->halo_type == CS_HALO_EXTENDED) {

            start_ext = halo->perio_lst[shift + 4*rank_id + 2];
            end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

            for (i = start_ext; i < end_ext; i++) {
              grad[n_cells + i][0] = _drdxyz[0 + d_var*3 + 18*i];
              grad[n_cells + i][1] = _drdxyz[1 + d_var*3 + 18*i];
              grad[n_cells + i][2] = _drdxyz[2 + d_var*3 + 18*i];
            }

          } /* End if extended halo */

        } /* End of loop on ranks */

      } /* End of test on rotation */

    } /* End of loop on transformations */

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize ghost cell values for Reynolds stress tensor gradient.
 *
 * We retrieve the gradient given by perinr (phyvar) for the Reynolds
 * stress tensor in a buffer on ghost cells. then we define
 * dpdx, dpdy and dpdz gradient (1 -> n_cells_with_ghosts).
 *
 * We can't implicitly take into account rotation of a gradient of a tensor
 * variable because we have to know all components.
 *
 * We set idimtr to 2 for the Reynolds stress tensor.
 *
 * We assume that is correct to treat periodicities implicitly for the other
 * variables when reconstructing gradients.
 *
 * \param[out]      tr_dim   2 for tensor (Rij) in case of rotation, 0 otherwise
 * \param[in, out]  grad     gradient of field
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_perio_init_rij_tensor(int           *tr_dim,
                                  cs_real_63_t   grad[])
{
  cs_mesh_t  *mesh = cs_glob_mesh;
  int isou;

  if (mesh->halo == NULL ) {
    *tr_dim = 0;
    return;
  }
  else
    *tr_dim = 2;

  /*
    When there is periodicity of rotation :
      - Retrieve gradient values for ghost cells and for the previous
        time step without reconstruction
      - Ghost cells without rotation keep their value.
  */

  if (_drdxyz != NULL) {

    int  rank_id, t_id;
    cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

    cs_halo_t  *halo = mesh->halo;

    const cs_lnum_t  n_cells   = mesh->n_cells;
    const cs_lnum_t  n_transforms = mesh->n_transforms;
    const fvm_periodicity_t  *periodicity = mesh->periodicity;

    for (t_id = 0; t_id < n_transforms; t_id++) {

      if (   fvm_periodicity_get_type(periodicity, t_id)
          >= FVM_PERIODICITY_ROTATION) {

        shift = 4 * halo->n_c_domains * t_id;

        for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

          start_std = halo->perio_lst[shift + 4*rank_id];
          end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

          for (i = start_std; i < end_std; i++) {
            for (isou = 0; isou < 6 ; isou++){
              grad[n_cells + i][isou][0] = _drdxyz[0 + isou*3 + 18*i];
              grad[n_cells + i][isou][1] = _drdxyz[1 + isou*3 + 18*i];
              grad[n_cells + i][isou][2] = _drdxyz[2 + isou*3 + 18*i];
            }
          }

          if (mesh->halo_type == CS_HALO_EXTENDED) {

            start_ext = halo->perio_lst[shift + 4*rank_id + 2];
            end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

            for (i = start_ext; i < end_ext; i++) {
              for (isou = 0; isou < 6 ; isou++){
                grad[n_cells + i][isou][0] = _drdxyz[0 + isou*3 + 18*i];
                grad[n_cells + i][isou][1] = _drdxyz[1 + isou*3 + 18*i];
                grad[n_cells + i][isou][2] = _drdxyz[2 + isou*3 + 18*i];
              }
            }

          } /* End if extended halo */

        } /* End of loop on ranks */

      } /* End of test on rotation */

    } /* End of loop on transformations */

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Process grad buffers in case of rotation on Reynolds stress tensor.
 *
 * We retrieve the gradient given by cs_gradient_perio_init_rij (phyvar)
 * for the Reynolds stress tensor in a buffer on ghost cells. Then we define
 * grad gradient (1 -> n_cells_with_ghosts).
 *
 * We can't implicitly take into account rotation of a gradient of a tensor
 * variable because we have to know all components.
 *
 * We assume that is correct to treat periodicities implicitly for the other
 * variables when reconstructing gradients.
 *
 * \param[in]       f_id     field index
 * \param[in, out]  grad     gradient of field
 *
 * size of _drdxyz and _wdrdxy = n_ghost_cells*6*3
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_perio_process_rij(const cs_int_t    *f_id,
                              cs_real_3_t        grad[])
{
  cs_lnum_t  d_var = -1;
  cs_mesh_t  *mesh = cs_glob_mesh;

  cs_field_t  *f = cs_field_by_id(*f_id);

  if (f->name[0] == 'r') {
    if (strlen(f->name) == 3) {
      if (f->name[1] == '1') {
        if (f->name[2] == '1')
          d_var = 0;
        else if (f->name[2] == '2')
          d_var = 3;
        else if (f->name[2] == '3')
          d_var = 4;
      }
      else if (f->name[1] == '2') {
        if (f->name[2] == '2')
          d_var = 1;
        else if (f->name[2] == '3')
          d_var = 5;
      }
      else if (f->name[1] == '3' && f->name[2] == '3')
        d_var = 2;
    }
  }

  if (mesh->halo == NULL || d_var < 0) {
    return;
  }

  /*
    When there is periodicity of rotation :
      - Retrieve gradient values for ghost cells and for the previous
        time step without reconstruction
      - Ghost cells without rotation keep their value.
  */

  if (_drdxyz != NULL) {

    int  rank_id, t_id;
    cs_lnum_t  i, shift, start_std, end_std, start_ext, end_ext;

    cs_halo_t  *halo = mesh->halo;
    const cs_lnum_t  n_cells   = mesh->n_cells;
    const cs_lnum_t  n_transforms = mesh->n_transforms;
    const fvm_periodicity_t  *periodicity = mesh->periodicity;

    for (t_id = 0; t_id < n_transforms; t_id++) {

      if (   fvm_periodicity_get_type(periodicity, t_id)
          >= FVM_PERIODICITY_ROTATION) {

        shift = 4 * halo->n_c_domains * t_id;

        for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

          start_std = halo->perio_lst[shift + 4*rank_id];
          end_std = start_std + halo->perio_lst[shift + 4*rank_id + 1];

          for (i = start_std; i < end_std; i++) {
            grad[n_cells + i][0] = _drdxyz[0 + d_var*3 + 18*i];
            grad[n_cells + i][1] = _drdxyz[1 + d_var*3 + 18*i];
            grad[n_cells + i][2] = _drdxyz[2 + d_var*3 + 18*i];
          }

          if (mesh->halo_type == CS_HALO_EXTENDED) {

            start_ext = halo->perio_lst[shift + 4*rank_id + 2];
            end_ext = start_ext + halo->perio_lst[shift + 4*rank_id + 3];

            for (i = start_ext; i < end_ext; i++) {
              grad[n_cells + i][0] = _drdxyz[0 + d_var*3 + 18*i];
              grad[n_cells + i][1] = _drdxyz[1 + d_var*3 + 18*i];
              grad[n_cells + i][2] = _drdxyz[2 + d_var*3 + 18*i];
            }

          } /* End if extended halo */

        } /* End of loop on ranks */

      } /* End of test on rotation */

    } /* End of loop on transformations */

  }
}
/*----------------------------------------------------------------------------*/

END_C_DECLS
