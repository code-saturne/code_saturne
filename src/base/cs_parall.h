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
 * subroutine parcmx (counter)
 * *****************
 *
 * integer          counter       <-> : input = local counter
 *                                      output = global max counter
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcmx, PARCMX)(cs_int_t  *counter);

/*----------------------------------------------------------------------------
 * Compute the minimum value of a counter (int) for the entire domain in
 * case of parallelism.
 *
 * Fortran Interface
 *
 * subroutine parcmn (counter)
 * *****************
 *
 * integer          counter       <-> : input = local counter
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
 * subroutine parcpt (counter)
 * *****************
 *
 * integer          counter     : <-> : input = counter to sum
 *                                      output = global sum
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parcpt, PARCPT)(cs_int_t  *counter);

/*----------------------------------------------------------------------------
 * Compute the global sum of a real for the entire domain in case of parellism
 *
 * Fortran Interface :
 *
 * subroutine parsom (var)
 * *****************
 *
 * double precision var         : <-> : input = value to sum
 *                                      output = global sum
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parsom, PARSOM)(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Compute the maximum value of a real variable for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * subroutine parmax (var)
 * *****************
 *
 * double precision var         : <-> : input = local maximum value
 *                                      output = global maximum value
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmax, PARMAX)(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Compute the minimum value of a real variable for the entire domain in case
 * of parallelism.
 *
 * Fortran Interface :
 *
 * subroutine parmin (var)
 * *****************
 *
 * double precision var         : <-> : input = local minimum value
 *                                      output = global minimum value
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parmin, PARMIN)(cs_real_t  *var);

/*----------------------------------------------------------------------------
 * Maximum value of a real and its related rank for the entire domain in
 * case of parallelism.
 *
 * Fortran Interface
 *
 * subroutine parmxl (nbr, var, xyzvar)
 * *****************
 *
 * integer          nbr         : <-- : size of the related variable
 * double precision var         : <-> : input: local max. value
 *                                      output: global max. value
 * double precision xyzvar(nbr) : <-> : input: value related to local max.
 *                                      output: value related to global max.
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
 * subroutine parmnl (nbr, var, xyzvar)
 * *****************
 *
 * integer          nbr         : <-- : size of the related variable
 * double precision var         : <-> : input: local max. value
 *                                      output: global max. value
 * double precision xyzvar(nbr) : <-> : input: value related to local max.
 *                                      output: value related to global max.
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
 * subroutine parism (n_elts, array)
 * *****************
 *
 * integer          n_elts       : <-- : size of the array.
 * integer          array(*)     : <-> : input = local array
 *                                       output = array of global sum values.
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
 * subroutine parimx (n_elts, array)
 * *****************
 *
 * integer          n_elts       : <-- : size of the array.
 * integer          array(*)     : <-> : input = local array
 *                                       output = array of global max. values.
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
 * subroutine parimn (n_elts, array)
 * *****************
 *
 * integer          n_elts       : <-- : size of the array.
 * integer          array(*)     : <-> : input = local array
 *                                       output = array of global min. values.
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
 * subroutine parrsm (n_elts, array)
 * *****************
 *
 * integer          n_elts       : <-- : size of the array.
 * double precision array(*)     : <-> : input = local array
 *                                       output = array of global sum values.
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
 * subroutine parrmx (n_elts, array)
 * *****************
 *
 * integer          n_elts        : <-- : size of the array
 * double precision array(*)      : <-> : input = local array
 *                                        output = array of global max. values.
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
 * subroutine parrmn (n_elts, array)
 * *****************
 *
 * integer          n_elts        : <-- : size of the array
 * double precision array(*)      : <-> : input = local array
 *                                        output = array of global min. values.
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
 * subroutine parbci (irank, n_elts, array)
 * *****************
 *
 * integer          irank       : <-- : rank related to the sending process
 * integer          n_elts      : <-- : size of the array
 * integer          array(*)    : <-> : array to broadcast
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
 * subroutine parbcr (irank, n_elts, array)
 * *****************
 *
 * integer            irank     : <-- : rank related to the sending process
 * integer            n_elts    : <-- : size of the array
 * double precision   array(*)  : <-> : array to broadcast
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
 * subroutine paragv (nvar, nvargb, var, vargb)
 * *****************
 *
 * integer           n_elts      : <-- : size of the local array
 * integer           n_g_elts    : <-- : size of the global array
 * double precision  array(*)    : <-- : local array
 * double precision  g_array(*)  : --> : global array
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
 * subroutine parfpt (node, ndrang, dis2mn)
 * *****************
 *
 * integer          node        : <-> : local number of the closest node
 * integer          ndrang      : --> : rank id for which the distance is the
 *                                      smallest
 * double precision dis2mn      : <-- : square distance between the closest node
 *                                      and the wanted node.
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
 * subroutine parhis (node, ndrang, var, varcap)
 * *****************
 *
 * integer          node        : <-- : local number of the element related to
 *                                      a measure node
 * integer          ndrang      : <-- : rank of the process owning the closest
 *                                      node from the measure node
 * double precision var(*)      : <-- : values of the variable on local elements
 * double precision varcap      : --> : value of the variable for the element
 *                                      related to the measure node
 *----------------------------------------------------------------------------*/

void
CS_PROCF (parhis, PARHIS)(cs_int_t   *node,
                          cs_int_t   *ndrang,
                          cs_real_t   var[],
                          cs_real_t  *varcap);

/*----------------------------------------------------------------------------
 * Call a barrier in case of parallelism
 *
 * This function should not be necessary in production code,
 * but it may be useful for debugging purposes.
 *
 * Fortran interface :
 *
 * subroutine parbar
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
