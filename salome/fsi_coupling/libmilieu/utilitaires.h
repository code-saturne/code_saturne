#ifndef __UTILITAIRES_H__
#define __UTILITAIRES_H__

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#ifdef __cplusplus
extern "C" {
#endif

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Predict displacement or forces based on values of the current and
 * previous time step(s)
 *
 * valpre = c1 * val1 + c2 * val2 + c3 * val3
 *----------------------------------------------------------------------------*/

void
pred(double *valpre,
     double *val1,
     double *val2,
     double *val3,
     double c1,
     double c2,
     double c3,
     int nbpts);

/*----------------------------------------------------------------------------
 * Compute the L2 norm of the difference between vectors vect1 and vect2
 *
 * dinorm = sqrt(sum on nbpts i
 *                 (sum on component j
 *                    ((vect1[i,j]-vect2[i,j])^2)))
 *----------------------------------------------------------------------------*/

double
dinorm(double *vect1,
       double *vect2,
       double nbpts);

/*----------------------------------------------------------------------------
 * Convergence test for implicit calculation case
 *
 * returns:
 *   0 if not converged
 *   1 if     converged
 *----------------------------------------------------------------------------*/

int
conv(int *icv);

/*----------------------------------------------------------------------------
 * Overwrites data from sub-iteration k-1 with data from sub-iteration k
 * dynamic data: velocities
 * efforts:      forces
 *----------------------------------------------------------------------------*/

void
val_ant(void);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif /* __UTILITAIRES_H__ */
