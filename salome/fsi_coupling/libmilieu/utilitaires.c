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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "runmilieu.h"
#include "donnees.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "utilitaires.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Public function definitions
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
     int nbpts)
{
  /* Local variables */
  int i;

  /* Update prediction array */
  for (i = 0; i < nbpts; i++) {
    valpre[3*i]     = c1*val1[3*i]     + c2*val2[3*i]     + c3*val3[3*i];
    valpre[(3*i)+1] = c1*val1[(3*i)+1] + c2*val2[(3*i)+1] + c3*val3[(3*i)+1];
    valpre[(3*i)+2] = c1*val1[(3*i)+2] + c2*val2[(3*i)+2] + c3*val3[(3*i)+2];
  }
}

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
       double nbpts)
{
  /* Local variables */
  int i;
  double norme;

  /* Compute the norm of the difference */
  norme = 0.;
  for (i = 0; i < nbpts; i++) {
    norme = norme + (vect1[3*i]-vect2[3*i])*(vect1[3*i]-vect2[3*i]);
    norme = norme + (vect1[3*i+1]-vect2[3*i+1])*(vect1[3*i+1]-vect2[3*i+1]);
    norme = norme + (vect1[3*i+2]-vect2[3*i+2])*(vect1[3*i+2]-vect2[3*i+2]);
  }
  norme = sqrt(norme/nbpts);
  return norme;
}

/*----------------------------------------------------------------------------
 * Allocate and initialize dynamic vectors (double) based on the 'nb_dyn'
 * number of points.
 *----------------------------------------------------------------------------*/

void
alldyn(void)
{
  /* Local variables */
  int k;

  xast =(double *)calloc(3*nb_dyn,sizeof(double));
  xvast=(double *)calloc(3*nb_dyn,sizeof(double));
  xvasa=(double *)calloc(3*nb_dyn,sizeof(double));
  xastp=(double *)calloc(3*nb_dyn,sizeof(double));

  for (k = 0; k < nb_dyn; k++) {

    xast[3*k]   = 0.;
    xast[3*k+1] = 0.;
    xast[3*k+2] = 0.;

    xvast[3*k]   = 0.;
    xvast[3*k+1] = 0.;
    xvast[3*k+2] = 0.;

    xvasa[3*k]   = 0.;
    xvasa[3*k+1] = 0.;
    xvasa[3*k+2] = 0.;

    xastp[3*k]   = 0.;
    xastp[3*k+1] = 0.;
    xastp[3*k+2] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * Allocate and initialize dynamic vectors (double) based on the 'nb_for'
 * number of points.
 *----------------------------------------------------------------------------*/

void
allfor(void)
{
  /* Local variables */
  int k;

  foras =(double *)calloc(3*nb_for,sizeof(double));
  foaas =(double *)calloc(3*nb_for,sizeof(double));
  fopas =(double *)calloc(3*nb_for,sizeof(double));

  for (k = 0; k < nb_for; k++) {

    foras[3*k]   = 0.;
    foras[3*k+1] = 0.;
    foras[3*k+2] = 0.;

    foaas[3*k]   = 0.;
    foaas[3*k+1] = 0.;
    foaas[3*k+2] = 0.;

    fopas[3*k]   = 0.;
    fopas[3*k+1] = 0.;
    fopas[3*k+2] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * Convergence test for implicit calculation case
 *
 * returns:
 *   0 if not converged
 *   1 if     converged
 *----------------------------------------------------------------------------*/

int
conv(int *icv)
{
  /* Local variables */
  int iret;
  double delast = 0.;

  if (lref > 0.) {

    delast = (dinorm(xast, xastp, nb_dyn))/lref;

    printf("--------------------------------\n");
    printf("convergence test:\n");
    printf("delast = %4.2le\n", delast);

    if (delast <= epsilo) {
      *(icv) = 1;
      printf("icv = %i\n", *(icv));
      printf("convergence of sub iteration\n");
      printf("--------------------------------\n");
    }
    else {
      printf("icv = %i\n", *(icv));
      printf("non convergence of sub iteration\n");
      printf("--------------------------------\n");
    }

    iret = 0;
  }
  else {
    printf("Value of lref is negative or zero\n");
    printf("calculation is aborted\n");
    printf("---------------------------------\n");
    iret = -1;
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Overwrites data from sub-iteration k-1 with data from sub-iteration k
 * dynamic data: velocities
 * efforts:      forces
 *----------------------------------------------------------------------------*/

void
val_ant(void)
{
  /* Local variables */
  int i;

  /* record efforts */
  for (i = 0; i< nb_for; i++) {
    foaas[3*i]   = foras[3*i];
    foaas[3*i+1] = foras[3*i+1];
    foaas[3*i+2] = foras[3*i+2];
  }

  /* record dynamic data */
  for (i = 0; i< nb_dyn; i++) {
    xvasa[3*i]   = xvast[3*i];
    xvasa[3*i+1] = xvast[3*i+1];
    xvasa[3*i+2] = xvast[3*i+2];
  }
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif
