#ifndef __CS_PARALL_H__
#define __CS_PARALL_H__

/*============================================================================
 * Functions dealing with parallelism
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the maximum value of a counter (int) for the entire domain in
 * case of parallelism.
 *
 * Fortran Interface
 *
 * SUBROUTINE PARCMX (IND)
 * *****************
 *
 * INTEGER          COUNTER       <-> : input = local counter
 *                                      output = global max counter
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcmx, PARCMX)(cs_int_t  *counter);

/*----------------------------------------------------------------------------
 * Compute the minimum value of a counter (int) for the entire domain in
 * case of parallelism.
 *
 * Fortran Interface
 *
 * SUBROUTINE PARCMN (IND)
 * *****************
 *
 * INTEGER          COUNTER       <-> : input = local counter
 *                                      output = global min counter
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcmn, PARCMN)(cs_int_t  *counter);

/*----------------------------------------------------------------------------
 * Compute the global sum of a counter (int) for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARCPT (CPT)
 * *****************
 *
 * INTEGER          COUNTER     : <-> : input = counter to sum
 *                                      output = global sum
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcpt, PARCPT)(cs_int_t  *counter);

/*----------------------------------------------------------------------------
 * Compute the global sum of a real for the entire domain in case of parellism
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARSOM (VAR)
 * *****************
 *
 * DOUBLE PRECISION VAR         : <-> : input = value to sum
 *                                      output = global sum
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parsom, PARSOM)(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Compute the maximum value of a real variable for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARMAX (VAR)
 * *****************
 *
 * DOUBLE PRECISION VAR         : <-> : input = local maximum value
 *                                      output = global maximum value
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmax, PARMAX)(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Compute the minimum value of a real variable for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARMIN (VAR)
 * *****************
 *
 * DOUBLE PRECISION VAR         : <-> : input = local minimum value
 *                                      output = global minimum value
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmin, PARMIN)(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Maximum value of a real and its related rank for the entire domain in
 * case of parallelism.
 *
 * Fortran Interface
 *
 * SUBROUTINE PARMXL (NBR, VAR, XYZVAR)
 * *****************
 *
 * INTEGER          NBR         :  -> : nombre de valeurs associees
 * DOUBLE PRECISION VAR         : <-> : input: local max. value
 *                                      output: global max. value
 * DOUBLE PRECISION XYZVAR(NBR) : <-> : input: value related to local max.
 *                                      output: value related to global max.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmxl, PARMXL)(cs_int_t   *nbr,
                          cs_real_t  *var,
                          cs_real_t   xyzvar[]);

/*----------------------------------------------------------------------------
 * Minimum value of a real and its related rank for the entire domain in
 * case of parallelism.
 *
 * Fortran Interface
 *
 * Interface Fortran :
 *
 * SUBROUTINE PARMNL (NBR, VAR, XYZVAR)
 * *****************
 *
 * INTEGER          NBR         : --> : size of the related variable
 * DOUBLE PRECISION VAR         : <-> : input: local max. value
 *                                      output: global max. value
 * DOUBLE PRECISION XYZVAR(NBR) : <-> : input: value related to local max.
 *                                      output: value related to global max.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmnl, PARMNL)(cs_int_t   *nbr,
                          cs_real_t  *var,
                          cs_real_t   xyzvar[]);

/*----------------------------------------------------------------------------
 * Compute the global sum for each element of an array of int in case of
 * parallelism.
 *
 * Fortran Interface
 *
 * SUBROUTINE PARISM (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS       : --> : size of the array.
 * INTEGER          ARRAY(*)     : <-> : input = local array
 *                                       output = array of global sum values.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parism, PARISM)(cs_int_t  *n_elts,
                          cs_int_t   array[]);

/*----------------------------------------------------------------------------
 * Compute the global maximum value for each element of an array of int in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARIMX (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS       : --> : size of the array.
 * INTEGER          ARRAY(*)     : <-> : input = local array
 *                                       output = array of global max. values.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parimx, PARIMX)(cs_int_t  *n_elts,
                          cs_int_t   array[]);

/*----------------------------------------------------------------------------
 * Compute the global minimum value for each element of an array of int in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARIMN (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS       : --> : size of the array.
 * INTEGER          ARRAY(*)     : <-> : input = local array
 *                                       output = array of global min. values.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parimn, PARIMN)(cs_int_t  *n_elts,
                          cs_int_t   array[]);

/*----------------------------------------------------------------------------
 * Compute the global sum for each element of an array of real in case of
 * parallelism.
 *
 * Fortran Interface
 *
 * SUBROUTINE PARRSM (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS       : --> : size of the array.
 * DOUBLE PRECISION ARRAY(*)     : <-> : input = local array
 *                                       output = array of global sum values.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parrsm, PARRSM)(cs_int_t   *n_elts,
                          cs_real_t   array[]);

/*----------------------------------------------------------------------------
 * Compute the global maximum value for each element of an array of real in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARRMX (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS        : --> : size of the array
 * DOUBLE PRECISION ARRAY(*)      : <-> : input = local array
 *                                        output = array of global max. values.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parrmx, PARRMX)(cs_int_t   *n_elts,
                          cs_real_t   array[]);

/*----------------------------------------------------------------------------
 * Compute the global minimum value for each element of an array of real in
 * case of parallelism.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARRMN (N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          N_ELTS        : --> : size of the array
 * DOUBLE PRECISION ARRAY(*)      : <-> : input = local array
 *                                        output = array of global min. values.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parrmn, PARRMN)(cs_int_t   *n_elts,
                          cs_real_t   array[]);

/*----------------------------------------------------------------------------
 * Broadcast to all the ranks the value of each element of an array of int.
 * (encapsulation of MPI_Bcast())
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARBCI (IRANK, N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER          IRANK       : --> : rank related to the sending process
 * INTEGER          N_ELTS      : --> : size of the array
 * INTEGER          ARRAY(*)    : <-> : array of int
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parbci, PARBCI)(cs_int_t   *irank,
                          cs_int_t   *n_elts,
                          cs_int_t    array[]);

/*----------------------------------------------------------------------------
 * Broadcast to all the ranks the value of each element of an array of real.
 * (encapsulation of MPI_Bcast())
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARBCR (IRANK, N_ELTS, ARRAY)
 * *****************
 *
 * INTEGER            IRANK     : --> : rank related to the sending process
 * INTEGER            N_ELTS    : --> : size of the array
 * DOUBLE PRECISION   ARRAY(*)  : <-> : array of real
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parbcr, PARBCR)(cs_int_t   *irank,
                          cs_int_t   *n_elts,
                          cs_real_t   array[]);

/*----------------------------------------------------------------------------
 * Build a global array from each local array in each domain. The size of
 * each local array can be different.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARAGV (NVAR, NVARGB, VAR, VARGB)
 * *****************
 *
 * INTEGER           N_ELTS      : --> : size of the local array
 * INTEGER           N_G_ELTS    : --> : size of the global array
 * DOUBLE PRECISION  ARRAY(*)    : --> : local array
 * DOUBLE PRECISION  G_ARRAY(*)  : <-- : global array
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (paragv, PARAGV)(cs_int_t   *n_elts,
                          cs_int_t   *n_g_elts,
                          cs_real_t   array[],
                          cs_real_t  *g_array);

/*----------------------------------------------------------------------------
 * Find a node which minimizes a given distance and its related rank.
 * May be used to locate a node among several domains.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARFPT (NODE, NDRANG, DIS2MN)
 * *****************
 *
 * INTEGER          NODE        : <-> : local number of the closest node
 * INTEGER          NDRANG      : <-- : rank id for which the distance is the
 *                                      smallest
 * DOUBLE PRECISION DIS2MN      : --> : square distance between the closest node
 *                                      and the wanted node.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parfpt, PARFPT)(cs_int_t   *node,
                          cs_int_t   *ndrang,
                          cs_real_t  *dis2mn);

/*----------------------------------------------------------------------------
 * Return the value associated to a probe.
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARHIS (NODE, NDRANG, VAR, VARCAP)
 * *****************
 *
 * INTEGER          NODE        : --> : local number of the element related to
 *                                      a measure node
 * INTEGER          NDRANG      : --> : rank of the process owning the closest
 *                                      node from the measure node
 * DOUBLE PRECISION VAR(*)      : --> : values of the variable on local elements
 * DOUBLE PRECISION VARCAP      : <-- : value of the variable for the element
 *                                      related to the measure node
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parhis, PARHIS)(cs_int_t   *node,
                          cs_int_t   *ndrang,
                          cs_real_t   var[],
                          cs_real_t  *varcap);

/*----------------------------------------------------------------------------
 * Return the global cell number of a local cell.
 * (send to all the processes)
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARCEL (LNUM, RANKID, GNUM)
 * *****************
 *
 * INTEGER          LNUM        : --> : local cell number
 * INTEGER          RANKID      : --> : rank of the domain (0 to N-1)
 * INTEGER          GNUM        : <-- : global cell number
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcel, PARCEL)(cs_int_t   *lnum,
                          cs_int_t   *rankid,
                          cs_int_t   *gnum);

/*----------------------------------------------------------------------------
 * Return the global cell number knowing its related local cell number. No
 * communication is useful.
 * Return the local cell number in serial mode.
 * Return 0 if the local cell number > mesh->n_cells
 * Return 0 if the current rank domain != RANKID
 *
 * Fortran interface :
 *
 * SUBROUTINE PARCLG (LNUM, RANKID, GNUM)
 * *****************
 *
 * INTEGER          LNUM      : --> : local cell number
 * INTEGER          RANKID    : --> : rank of the current domain (0 to N-1)
 * INTEGER          GNUM      : <-- : global cell number
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parclg, PARCLG)(cs_int_t   *lnum,
                          cs_int_t   *rankid,
                          cs_int_t   *gnum);

/*----------------------------------------------------------------------------
 * Return the global internal face number knowing its related local internal
 * face number. No communication is useful.
 * Return the local internal face number in serial mode.
 * Return 0 if the local internal face number > mesh->n_i_faces
 * Return 0 if the current rank domain != RANKID
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARFIG (LNUM, RANKID, GNUM)
 * *****************
 *
 * INTEGER          LNUM      : --> : local internal face number
 * INTEGER          RANKID    : --> : rank of the current domain (0 to N-1)
 * INTEGER          GNUM      : <-- : global internal face number
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parfig, PARFIG)(cs_int_t   *lnum,
                          cs_int_t   *rankid,
                          cs_int_t   *gnum);

/*----------------------------------------------------------------------------
 * Return the global border face number knowing its related local border face
 * number. No communication is useful.
 * Return the local border face number in serial mode.
 * Return 0 if the local border face number > mesh->n_b_faces
 * Return 0 if the current rank domain != RANKID
 *
 * Fortran Interface :
 *
 * SUBROUTINE PARFBG (LNUM, RANKID, GNUM)
 * *****************
 *
 * INTEGER          LNUM      : --> : local border face number
 * INTEGER          RANKID    : --> : rank of the current domain (0 to N-1)
 * INTEGER          GNUM      : <-- : global border face number
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parfbg, PARFBG)(cs_int_t   *lnum,
                          cs_int_t   *rankid,
                          cs_int_t   *gnum);

/*----------------------------------------------------------------------------
 * Call a barrier in case of parallelism
 *
 * This function should not be necessary in production code,
 * but it may be useful for debugging purposes.
 *
 * Fortran interface :
 *
 * SUBROUTINE PARBAR
 * *****************
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parbar, PARBAR)(void);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARALL_H__ */
