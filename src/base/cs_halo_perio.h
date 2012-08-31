#ifndef __CS_HALO_PERIO_H__
#define __CS_HALO_PERIO_H__

/*============================================================================
 * Structure and function headers associated to periodicity
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
 * Process dpdx, dpdy, dpdz buffers in case of rotation on velocity vector and
 * Reynolds stress tensor.
 *
 * We retrieve the gradient given by perinu and perinr (phyvar) for the
 * velocity and the reynolds stress tensor in a buffer on ghost cells. then
 * we define dpdx, dpdy and dpdz gradient (1 -> n_cells_with_ghosts).
 *
 * We can't implicitly take into account rotation of a gradient of a non-scalar
 * variable because we have to know the all three components in grad*c.
 *
 * Otherwise, we can implicitly treat values given by translation. There will
 * be replace further in grad*c.
 *
 * We assume that is correct to treat periodicities implicitly for the other
 * variables.
 *
 * Fortran Interface:
 *
 * subroutine persyn
 * *****************
 *
 * integer          ivar         :  -> : numero de la variable
 * integer          iu           : <-- : position de la Vitesse(x,y,z)
 * integer          iv           : <-- : dans RTP, RTPA
 * integer          iw           : <-- :     "                   "
 * integer          itytur       : <-- : turbulence (Rij-epsilon ITYTUR = 3)
 * integer          ir11         : <-- : position des Tensions de Reynolds
 * integer          ir22         : <-- : en Rij dans RTP, RTPA
 * integer          ir33         : <-- :     "                   "
 * integer          ir12         : <-- :     "                   "
 * integer          ir13         : <-- :     "                   "
 * integer          ir23         : <-- :     "                   "
 * double precision dpdx(ncelet) : <-> : gradient de IVAR
 * double precision dpdy(ncelet) : <-> :    "        "
 * double precision dpdz(ncelet) : <-> :    "        "
 *----------------------------------------------------------------------------*/

void
CS_PROCF (persyn, PERSYN)(const cs_int_t    *ivar,
                          const cs_int_t    *iu,
                          const cs_int_t    *iv,
                          const cs_int_t    *iw,
                          const cs_int_t    *itytur,
                          const cs_int_t    *ir11,
                          const cs_int_t    *ir22,
                          const cs_int_t    *ir33,
                          const cs_int_t    *ir12,
                          const cs_int_t    *ir13,
                          const cs_int_t    *ir23,
                          cs_real_t          dpdx[],
                          cs_real_t          dpdy[],
                          cs_real_t          dpdz[]);

/*----------------------------------------------------------------------------
 * Rotate vector values for periodic cells on extended halos.
 *
 * Fortran API:
 *
 * subroutine perrve
 * *****************
 *
 * double precision var1(ncelet) : <-> : component 1 of vector
 * double precision var2(ncelet) : <-> : component 2 of vector
 * double precision var3(ncelet) : <-> : component 3 of vector
 *----------------------------------------------------------------------------*/

void
CS_PROCF (perrve, PERRVE) (cs_real_t  var1[],
                           cs_real_t  var2[],
                           cs_real_t  var3[]);

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
 * SUBROUTINE PERMAS
 * *****************
 *
 * INTEGER          IMASPE      :  -> : suivant l'appel de INIMAS
 *                                          = 1 si appel de RESOLP ou NAVSTO
 *                                          = 2 si appel de DIVRIJ
 * INTEGER          IMASPE      :  -> : indicateur d'appel dans INIMAS
 *                                          = 1 si appel au debut
 *                                          = 2 si appel a la fin
 * DOUBLE PRECISION ROM(NCELET) :  -> : masse volumique aux cellules
 * DOUBLE PRECISION DUDXYZ      :  -> : gradient de U aux cellules halo pour
 *                                      l'approche explicite en periodicite
 * DOUBLE PRECISION DRDXYZ      :  -> : gradient de R aux cellules halo pour
 *                                      l'approche explicite en periodicite
 * DOUBLE PRECISION WDUDXY      :  -  : tableau de travail pour DUDXYZ
 * DOUBLE PRECISION WDRDXY      :  -  : tableau de travail pour DRDXYZ
 *
 * Size of DUDXYZ and WDUDXY = n_ghost_cells*3*3
 * Size of DRDXYZ and WDRDXY = n_ghost_cells*6*3
 *----------------------------------------------------------------------------*/

void
CS_PROCF (permas, PERMAS)(const cs_int_t    *imaspe,
                          const cs_int_t    *iappel,
                          cs_real_t          rom[],
                          cs_real_t         *dudxyz,
                          cs_real_t         *drdxyz,
                          cs_real_t         *wdudxy,
                          cs_real_t         *wdrdxy);

/*----------------------------------------------------------------------------
 * Process dpdx, dpdy, dpdz buffers in case of rotation on velocity vector and
 * Reynolds stress tensor.
 *
 * We retrieve the gradient given by perinu and perinr (phyvar) for the
 * velocity and the Reynolds stress tensor in a buffer on ghost cells. then
 * we define dpdx, dpdy and dpdz gradient (1 -> n_cells_with_ghosts).
 *
 * We can't implicitly take into account rotation of a gradient of alnon-scalar
 * variable because we have to know the all three components in gradrc.
 *
 * Otherwise, we can implicitly treat values given by translation. There will
 * be replace further in GRADRC.
 *
 * We set idimtr to 1 and 2 respectively for the velocity vector and the
 * Reynolds stress tensor.
 *
 * We assume that is correct to treat periodicities implicitly for the other
 * variables in gradrc, for which we set idimtr to 0.
 *
 * Fortran Interface:
 *
 * subroutine pering
 * *****************
 *
 * integer          ivar         : <-- : variable number
 * integer          idimtr       : <-- : 0 if ivar does not match a vector
 *                               :     :   or tensor or there is no periodicity
 *                               :     :   of rotation
 *                               :     : 1 for velocity, 2 for Reynolds stress
 *                               :     :   in case of periodicity of rotation
 * integer          irpvar       :     : -1 if ivar does not match a vector or
 *                               :     :   or tensor or there is no periodicity
 *                               :     :   of rotation; otherwise:
 *                               :     : 0 for iu, 1 for iv, 2 for iw
 *                               :     : 0 for ir11, 1 for ir22, 2 for ir33,
 *                               :     : 3 for ir12, 4 for ir13, 5 for ir23
 * integer          iguper       : <-- : 0/1 indicates we have not computed
 *                               :     :   gradients in dudxyz
 * integer          igrper       : <-- : 0/1 indicates we have not computed
 *                               :     :   gradients in drdxyz
 * double precision dpdx(ncelet) : <-> : gradient of ivar
 * double precision dpdy(ncelet) : <-> :    "        "
 * double precision dpdz(ncelet) : <-> :    "        "
 * double precision dudxyz       :  -> : gradient of u at ghost cells for
 *                                       explicit periodicity of rotation
 * double precision drdxyz       :  -> : gradient of r at ghost cells for
 *                                       explicit periodicity of rotation
 *
 * size of dudxyz and wdudxy = n_ghost_cells*3*3
 * size of drdxyz and wdrdxy = n_ghost_cells*6*3
 *----------------------------------------------------------------------------*/

void
CS_PROCF (pering, PERING)(const cs_int_t    *idimtr,
                          const cs_int_t    *irpvar,
                          const cs_int_t    *iguper,
                          const cs_int_t    *igrper,
                          cs_real_t          dpdx[],
                          cs_real_t          dpdy[],
                          cs_real_t          dpdz[],
                          const cs_real_t   *dudxyz,
                          const cs_real_t   *drdxyz);

/*----------------------------------------------------------------------------
 * Exchange buffers for PERINU
 *
 * Fortran Interface:
 *
 * SUBROUTINE PEINU1
 * *****************
 *
 * INTEGER          ISOU          :  -> : component of the velocity vector
 * DOUBLE PRECISION DUDXYZ        :  -> : gradient of the velocity vector
 *                                        for ghost cells and for an explicit
 *                                        treatment of the periodicity.
 * DOUBLE PRECISION W1..3(NCELET) :  -  : working buffers
 *
 * Size of DUDXYZ and WDUDXY = n_ghost_cells*3*3
 *----------------------------------------------------------------------------*/

void
CS_PROCF (peinu1, PEINU1)(const cs_int_t    *isou,
                          cs_real_t         *dudxyz,
                          cs_real_t          w1[],
                          cs_real_t          w2[],
                          cs_real_t          w3[]);

/*----------------------------------------------------------------------------
 * Apply rotation on DUDXYZ tensor.
 *
 * Fortran Interface:
 *
 * SUBROUTINE PEINU2 (VAR)
 * *****************
 *
 * DOUBLE PRECISION DUDXYZ        : <-> : gradient of the velocity vector
 *                                        for ghost cells and for an explicit
 *                                        treatment of the periodicity.
 *
 * Size of DUDXYZ and WDUDXY = n_ghost_cells*3*3
 *----------------------------------------------------------------------------*/

void
CS_PROCF (peinu2, PEINU2)(cs_real_t         *dudxyz);

/*----------------------------------------------------------------------------
 * Exchange buffers for PERINR
 *
 * Fortran Interface
 *
 * SUBROUTINE PEINR1 (VAR)
 * *****************
 *
 * INTEGER          ISOU          : -> : component of the Reynolds stress tensor
 * DOUBLE PRECISION DRDXYZ        : -> : gradient of the Reynolds stress tensor
 *                                       for ghost cells and for an explicit
 *                                       treatment of the periodicity.
 * DOUBLE PRECISION GRADX(NCELET)
 *                  GRADY(NCELET)
 *                  GRADZ(NCELET) : -  : x, y, z components of the gradient of
 *                                       the current component of the Reynolds
 *                                       stress tensor.
 *
 * Size of DRDXYZ and WDRDXY = n_ghost_cells*6*3
 *----------------------------------------------------------------------------*/

void
CS_PROCF (peinr1, PEINR1)(const cs_int_t    *isou,
                          cs_real_t         *drdxyz,
                          cs_real_t          gradx[],
                          cs_real_t          grady[],
                          cs_real_t          gradz[]);

/*----------------------------------------------------------------------------
 * Apply rotation on the gradient of Reynolds stress tensor
 *
 * Fortran Interface:
 *
 * SUBROUTINE PEINR2 (VAR)
 * *****************
 *
 * DOUBLE PRECISION DRDXYZ        :  -> : gradient of the Reynolds stress tensor
 *                                       for ghost cells and for an explicit
 *                                       treatment of the periodicity.
 *
 * Size of DRDXYZ and WDRDXY = n_ghost_cells*6*3
 *----------------------------------------------------------------------------*/

void
CS_PROCF (peinr2, PEINR2)(cs_real_t         *drdxyz);

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
cs_halo_perio_sync_var_tens(const cs_halo_t *halo,
                            cs_halo_type_t   sync_mode,
                            cs_real_t        var[]);

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
 * Synchronize values for a real diagonal tensor (interleaved)
 * between periodic cells.
 *
 * We only know the diagonal of the tensor.
 *
 * parameters:
 *   halo      <-> halo associated with variable to synchronize
 *   sync_mode --> kind of halo treatment (standard or extended)
 *   var       <-> diagonal tensor to update
 *----------------------------------------------------------------------------*/

void
cs_halo_perio_sync_var_diag(const cs_halo_t *halo,
                            cs_halo_type_t   sync_mode,
                            cs_real_t        var[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HALO_PERIO_H__ */

