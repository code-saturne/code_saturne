#ifndef __CS_HALO_PERIO_H__
#define __CS_HALO_PERIO_H__

/*============================================================================
 * Structure and function headers associated to periodicity
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 *  Public function header for Fortran API
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
                           cs_real_t  var33[]);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Apply transformation on coordinates.
 *
 * parameters:
 *   halo      <-> halo associated with coordinates to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   coords    --> coordinates on which transformation have to be done.
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_coords(const cs_halo_t  *halo,
                          cs_halo_type_t    sync_mode,
                          cs_real_t        *coords);

/*----------------------------------------------------------------------------
 * Synchronize values for a real vector (interleaved) between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> type of halo treatment (standard or extended)
 *   var       <-> vector to update
 *   incvar    <-- specifies the increment for the elements of var
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_vect(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            cs_real_t         var[],
                            int               incvar);

/*----------------------------------------------------------------------------
 * Synchronize values for a real vector between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   var_x     <-> component of the vector to update
 *   var_y     <-> component of the vector to update
 *   var_z     <-> component of the vector to update
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_vect_ni(const cs_halo_t     *halo,
                               cs_halo_type_t       sync_mode,
                               cs_real_t            var_x[],
                               cs_real_t            var_y[],
                               cs_real_t            var_z[]);

/*----------------------------------------------------------------------------
 * Synchronize values for a real tensor between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
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
cs_halo_perio_sync_var_tens_ni(const cs_halo_t *halo,
                               cs_halo_type_t   sync_mode,
                               cs_real_t        var11[],
                               cs_real_t        var12[],
                               cs_real_t        var13[],
                               cs_real_t        var21[],
                               cs_real_t        var22[],
                               cs_real_t        var23[],
                               cs_real_t        var31[],
                               cs_real_t        var32[],
                               cs_real_t        var33[]);

/*----------------------------------------------------------------------------
 * Synchronize values for a real tensor (interleaved) between periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   var       <-> tensor to update
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_tens(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            cs_real_t         var[]);

/*----------------------------------------------------------------------------
 * Synchronize values for a real tensor (symmetric interleaved) between
 * periodic cells.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   var       <-> symmetric tensor to update (6 values)
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_sym_tens(const cs_halo_t  *halo,
                                cs_halo_type_t    sync_mode,
                                cs_real_t         var[]);

/*----------------------------------------------------------------------------
 * Synchronize values for a real diagonal tensor between periodic cells.
 *
 * We only know the diagonal of the tensor.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   var11     <-> component of the tensor to update
 *   var22     <-> component of the tensor to update
 *   var33     <-> component of the tensor to update
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_diag_ni(const cs_halo_t *halo,
                               cs_halo_type_t   sync_mode,
                               cs_real_t        var11[],
                               cs_real_t        var22[],
                               cs_real_t        var33[]);

/*----------------------------------------------------------------------------
 * Apply rotation on the gradient of Reynolds stress tensor
 *
 * parameters:
 *   drdxyz     <-> gradient on the variable (size: 3*6*n_ghost_cells)
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_rotate_rij(cs_real_t  *drdxyz);

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
                                     cs_real_t         var[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HALO_PERIO_H__ */

