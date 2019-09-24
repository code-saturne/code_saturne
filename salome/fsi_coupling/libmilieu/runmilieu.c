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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "utilitaires.h"
#include "donnees.h"
#include "communication.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "runmilieu.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 *  Global variables
 *============================================================================*/

int     nb_dyn = 0;
int     nb_for = 0;
int     ntcast = 0;
double  lref = 0.0;

double  *xast = NULL;
double  *xvast = NULL;
double  *xvasa = NULL;
double  *xastp = NULL;

double  *foras = NULL;
double  *foaas = NULL;
double  *fopas = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocate and initialize dynamic vectors (double) based on the 'nb_dyn'
 * number of points.
 *----------------------------------------------------------------------------*/

void
_allocate_arrays(void)
{
  xast =(double *)calloc(3*nb_dyn,sizeof(double));
  xvast=(double *)calloc(3*nb_dyn,sizeof(double));
  xvasa=(double *)calloc(3*nb_dyn,sizeof(double));
  xastp=(double *)calloc(3*nb_dyn,sizeof(double));

  for (int k = 0; k < nb_dyn; k++) {

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

  foras =(double *)calloc(3*nb_for,sizeof(double));
  foaas =(double *)calloc(3*nb_for,sizeof(double));
  fopas =(double *)calloc(3*nb_for,sizeof(double));

  for (int k = 0; k < nb_for; k++) {

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

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * "runmilieu" function
 *----------------------------------------------------------------------------*/

void
runmilieu(void *icompo)
{
  /* local variables */

  int i,j;
  int ierr = 0;
  int icv = 0;
  double c1, c2, c3;
  double alpha, beta;
  double dt_ast, dt_sat;
  double dtold = 0.;
  double dt = dtref;

  /* Input data for the coupled calculation are defined in SALOME by
     the study case's XML file */

  /* Initialize communication */
  if ((ierr = inicom(icompo)) >= 0) {
    printf(" Initializing communication\n");
  }

  /* Send parameters to codes */
  if ((ierr = send_param(icompo)) >= 0) {
    printf(" Send calculation parameters to codes\n");
  }

  /* Compute array sizes and initialize */

  /* Receive geometric data (nb_for and nb_dyn) */
  if (ierr >= 0) {
    ierr = recv_geom(icompo);

    printf("----------------------------------\n");
    printf(" Geometric parameters\n");
    printf("   number of coupled faces: %i\n", nb_for);
    printf("   number of coupled nodes: %i\n", nb_dyn);
    printf("   reference length (m): %4.2le\n", lref  );
    printf("----------------------------------\n");

    _allocate_arrays();

    /* Prediction coefficients */
    c1    = 0.;
    c2    = 0.;
    c3    = 0.;
    beta  = 0.;
    alpha = 0.;
  }

  /* Initialize time step */
  dt = 0.;
  dt_ast = 0.;
  dt_sat = 0.;

  /* Initialize coupling iteration */
  ntcast = 0;

 /* Main loop */
  i = 1;
  while (ierr >= 0) {
    printf("\n");
    printf("----------------------------------------------------\n");
    printf("\n");
    printf("*********************************\n");
    printf("*         iteration %i          *\n", i);
    printf("*********************************\n");
    printf("\n");
    printf("\n");

    /* Info on time scheme */
    if (nbssit <= 1) {
      printf("Explicit time-stepping scheme\n");
    }
    else {
      printf("Implicit time-stepping scheme\n");
      printf("  number of sub-iterations: %i\n", nbssit);
    }

    /* Manage time steps */

    /* Receive time steps from Code_Aster and Code_Saturne */
    ierr = recv_pdt(icompo,&(dt_ast), &(dt_sat), i);

    printf("----------------------------------\n");
    printf("reference time step:     %4.21e \n", dtref );
    printf("Code_Saturne time step:  %4.2le \n", dt_sat);
    printf("Code_Aster time step:    %4.2le \n", dt_ast);

    /* choose the smallest time step: dt = dt_ast; */
    dt = dtref;
    if (dt > dt_ast) {
      dt = dt_ast;
    }
    if (dt > dt_sat) {
      dt = dt_sat;
    }

    /* Send the selected time step */
    if (ierr >= 0) ierr = send_pdt(icompo, dt, i);

    printf("selected time step:      %4.2le \n", dt);
    printf("----------------------------------\n");
    printf("\n\n");

    icv = 0;

    j = 1;
    while (ierr >= 0) {

      printf("*********************************\n");
      printf("*     sub - iteration %i        *\n", j);
      printf("*********************************\n");
      printf("\n\n");

      /* increment coupling iteration */
      printf("midde\n");
      ntcast = ntcast + 1;
      printf("ntcast = %i\n", ntcast);

      /* printf("***************************************\n"); */
      /* printf("*        predict displacements        *\n"); */
      /* printf("***************************************\n"); */

      /* Predict displacements */

      c1 = 0.;
      c2 = 0.;
      c3 = 0.;

      /* separate prediction for explicit/implicit cases */
      if (j == 1) {
        alpha = 0.5;
        beta  = 0.;
        c1    = 1.;
        c2    = (alpha + beta) * dt ;
        c3    = -beta * dtold ;
        pred(xastp, xast, xvast, xvasa, c1, c2, c3, nb_dyn);
      }

      if (j > 1) {
        alpha = 0.5;
        c1    = alpha;
        c2    = 1. - alpha ;
        c3    = 0.;
        pred(xastp, xast, xastp, xast, c1, c2, c3, nb_dyn);
      }

      printf("--------------------------------------------\n");
      printf("Displacement prediction coefficients\n");
      printf(" C1: %4.2le\n", c1);
      printf(" C2: %4.2le\n", c2);
      printf(" C3: %4.2le\n", c3);
      printf("--------------------------------------------\n");
      printf("\n\n");

      /* send predicted displacements */
      if (ierr >= 0) ierr = send_dyn(icompo);

      /* explicit case: no need for a convergence test */

      /* implicit case: needs a convergence test */

      /* printf("***************************************\n"); */
      /* printf("*  end of displacements prediction    *\n"); */
      /* printf("***************************************\n"); */

      /* printf("*********************************\n"); */
      /* printf("*       forces prediction       *\n"); */
      /* printf("*********************************\n"); */

      /* Receive forces */
      ierr = recv_for(icompo);

      /* No difference between explicit and implicit cases for forces */
      alpha = 2.0;
      c1    = alpha;
      c2    = 1-alpha;
      c3    = 0.;
      pred(fopas, foras, foaas, foaas, c1, c2, c3, nb_for);
      printf("--------------------------------------\n");
      printf("Forces prediction coefficients\n");
      printf(" C1: %4.2le\n",c1);
      printf(" C2: %4.2le\n",c2);
      printf(" C3: %4.2le\n",c3);
      printf("--------------------------------------\n");
      printf("\n\n");

      /* send des forces */
      if (ierr >= 0) ierr = send_for(icompo);

      /* printf("*********************************\n"); */
      /* printf("*   end of forces prediction    *\n"); */
      /* printf("*********************************\n"); */

      printf("\n");

      /* explicit case: no need fo a convergence test */
      if (nbssit <= 1) {
        /* handle convergence even when no test is done */
        icv =1;
        if (ierr >= 0) ierr = send_icv1(icompo,icv);
        if (ierr >= 0) ierr = recv_icv(icompo,&(icv));
        icv = 1;
        if (ierr >= 0) ierr = send_icv2(icompo,icv);

        /* receive displacements effectively calculated by Code_Aster */
        if (ierr >= 0) ierr = recv_dyn(icompo);

        /* record previous values */
        val_ant();

        break;
      }

      /* implicit case: requires a convergence test */
      else {
        /* compute icv */
        if (ierr >= 0) ierr = conv(&(icv));
        if (ierr >= 0) ierr = send_icv1(icompo,icv);
        if (ierr >= 0) ierr = recv_icv(icompo,&(icv));
        if (ierr >= 0) ierr = send_icv2(icompo,icv);

        if((j>=nbssit) || (icv == 1)) {
          /* receive displacements effectivemey computed by Code_Aster */
          /* Receive displacements */
          if (ierr >= 0) ierr = recv_dyn(icompo);

          /* then send to Code_Saturne ? the question remains open... */
          /* if necessary, function to send these displs. should be
             created in middle and matching receive in Code_Saturne */
          /* if (ierr >= 0) ierr = send2_dyn(); */

          /* receive displacements effectiveley calculated by Code_Aster */
          if (ierr >= 0) ierr = recv_dyn(icompo);
          break;
        }
        else {
          j = j+1;
        }
      }
    } /* end of sub-iterations loop */

    /* iterations test */
    if (i >= nbpdtm) {
      ierr = -1;
    }
    /* end of iterations test */

    i = i+1;

    /* save time step */
    dtold = dt;

  } /* en of iterations loop */

  ierr = calfin(icompo);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif
