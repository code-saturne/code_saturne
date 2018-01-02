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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * SALOME headers
 *----------------------------------------------------------------------------*/

#include <calcium.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "runmilieu.h"
#include "donnees.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "communication.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Receive geometric data and sets variables nb_for, nb_dyn and lref
 *----------------------------------------------------------------------------*/

int
recv_geom(void *component)
{
  /* local variables */
  int    i = 0;
  int    ii = 0;
  int    iret = 0;
  char   nomvar[144];
  int    geom[2] = {0, 0};
  float  tps = 0.;
  double tps2=0.;
  double almloc = 0.;

  printf("in recv_geom \n");

  /* Initializations */
  nb_for = 0;
  nb_dyn = 0;
  lref   = 0.0;

  strcpy(nomvar, "DONGEO");

  /* Receive geometric variables:
   * 1 array of size 1 * 2 is received
   * geom[0] = nb_for; geom[1] = nb_dyn */

  i = 0;
  iret = cp_len(component,
                CP_ITERATION,
                &(tps),
                &(tps),
                &(i),
                nomvar,
                2,
                &(ii),
                &(geom[0]));

  if (iret < 1) {
    strcpy(nomvar, "ALMAXI");
    /* Receive lref variable */
    i = 0;
    iret = cp_ldb(component,
                  CP_ITERATION,
                  &(tps2),
                  &(tps2),
                  &(i),
                  nomvar,
                  1,
                  &(ii),
                  &(almloc));
    lref = almloc;
  }

  nb_for = geom[0];
  nb_dyn = geom[1];

  return iret;
}

/*----------------------------------------------------------------------------
 * Send geometric data to Code_Aster
 * To be removed in initialization stage
 *----------------------------------------------------------------------------*/

int
send_geom(void* component)
{
  /* Local variables */
  int iret = 0;
  char nomvar[] = "NB_DYN";

  /* Trace of call to this function */
  printf("in recv_geom\n");

  /* Send number of Code_Saturne points to Code_Aster */
  iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(nb_dyn));
  if (iret < 1) {
    strcpy(nomvar,"NB_FOR");
    /* send nbssit variable */
    iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(nb_for));
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Sends time step computed by middle component to codes
 *----------------------------------------------------------------------------*/

int
send_pdt(void *component,
         double dt,
         int numpdt)
{
  /* Local variables */
  int iret = 0;
  char nomvar[] = "DTCALC";

  /* Trace of call to this function */
  printf("in send_pdt\n");

  /* send dtref time step to codes */
  iret = cp_edb(component, CP_ITERATION, 0.0, numpdt, nomvar, 1, &(dt));

  return iret;
}

/*----------------------------------------------------------------------------
 * Receives time steps from Code_Aster and Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_pdt(void *component,
         double *dt_ast,
         double *dt_sat,
         int numpdt)
{
  /* Local variables */
  int    iret = 0;
  char   nomvar[] = "DTSAT";
  int    i, ii;
  double tps = 0.;
  double dtloc = 0.;

  /* Trace of call to this function */
  printf("in recv_pdt\n");

  i = numpdt;

  if (iret < 1) {
    strcpy(nomvar, "DTSAT");
    /* Receive Code_Saturne time step */
    i = numpdt;
    iret = cp_ldb(component,
                  CP_ITERATION,
                  &(tps),
                  &(tps),
                  &(i),
                  nomvar,
                  1,
                  &(ii),
                  &(dtloc));
    *(dt_sat) = dtloc;
  }

  if (iret < 1) {
    strcpy(nomvar, "DTSAT");
    /* Receive Code_Aster time step */
    i = numpdt;
    iret = cp_ldb(component,
                  CP_ITERATION,
                  &(tps),
                  &(tps),
                  &(i),
                  nomvar,
                  1,
                  &(ii),
                  &(dtloc));
    *(dt_ast)=dtloc;
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Sends the following parameters:
 *                                 nbpdtm
 *                                 nbssit
 *                                 epsilo
 *                                 isyncp
 *                                 ntchr
 *                                 ttpabs
 *----------------------------------------------------------------------------*/

int
send_param(void *component)
{
  /* Local variables */
  char nomvar[] = "NBPDTM";
  int iret = 0;

  printf("in send_param \n");

  /* Send nbpdtm */
  iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(nbpdtm));

  /* Send nbssit */
  if (iret < 1) {
    strcpy(nomvar, "NBSSIT");
    iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(nbssit));
  }

  /* Send epsilo */
  if (iret < 1) {
    strcpy(nomvar, "EPSILO");
    iret = cp_edb(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(epsilo));
  }

  /* Send isyncp */
  if (iret < 1) {
    strcpy(nomvar, "ISYNCP");
    iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(isyncp));
  }

  /* Send ntchr  */
  if (iret < 1) {
    strcpy(nomvar, "NTCHRO");
    iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(ntchr));
  }

  /* Send ttinit */
  if (iret < 1) {
    strcpy(nomvar,"TTINIT");
    iret = cp_edb(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(ttinit));
  }

  /* Send dtref */
  if (iret < 1) {
    strcpy(nomvar, "PDTREF");
    iret = cp_edb(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(dtref));
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Receives displacements and velocities from Code_Aster at current time step
 *----------------------------------------------------------------------------*/

int
recv_dyn(void *component)
{
  /* Local variables */
  char nomvar[] = "DEPAST";
  int iret = 0;
  int i, ii;
  double tps = 0.;

  printf("in recv_dyn \n");

  /* Receive depast */
  i = ntcast;
  iret = cp_ldb(component,
                CP_ITERATION,
                &(tps),
                &(tps),
                &(i),
                nomvar,
                3*nb_dyn,
                &(ii),
                xast);

  /* Receive vitast */
  if (iret < 1) {
    strcpy(nomvar, "VITAST");
    i = ntcast;
    iret = cp_ldb(component,
                  CP_ITERATION,
                  &(tps),
                  &(tps),
                  &(i),
                  nomvar,
                  3*nb_dyn,
                  &(ii),
                  xvast);
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Send predicted displacements to Code_Saturne
 *----------------------------------------------------------------------------*/

int
send_dyn(void *component)
{
  /* Local variables */
  char nomvar[] = "DEPSAT";
  int iret = 0;

  printf("in send_dyn \n");

  iret = cp_edb(component,
                CP_ITERATION,
                0.0,
                ntcast,
                nomvar,
                3*nb_dyn,
                xastp);

  return iret;
}

/*----------------------------------------------------------------------------
 * Receive efforts from Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_for(void *component)
{
  /* Local variables */
  char nomvar[] = "FORSAT";
  int iret = 0;
  double tps;
  int i, ii;

  printf("in recv_for \n");

  i = ntcast;
  iret = cp_ldb(component,
                CP_ITERATION,
                &(tps),
                &(tps),
                &(i),
                nomvar,
                3*nb_for,
                &(ii),
                foras);

  return iret;
}

/*----------------------------------------------------------------------------
 * Send predicted efforts to Code_Aster
 *----------------------------------------------------------------------------*/

int
send_for(void *component)
{
  /* Local variables */
  char nomvar[] = "FORAST";
  int iret = 0;
  int i;

  printf("in send_for \n");

  i = ntcast;
  iret = cp_edb(component, CP_ITERATION, 0.0, i, nomvar, 3*nb_for, fopas);

  return iret;
}

/*----------------------------------------------------------------------------
 * Send convergence indicator to Code_Saturne
 *----------------------------------------------------------------------------*/

int
send_icv1(void *component,
          int icv)
{
  /* Local variables */
  int iret;
  char nomvar[] = "ICVEXT";

  printf("in send_icv1 \n");

  /* Send icvext */
  iret =  cp_een(component, CP_ITERATION, 0.0, ntcast, nomvar, 1, &(icv));

  return iret;
}

/*----------------------------------------------------------------------------
 * Receive convergence indicator from Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_icv(void *component,
         int *icv)
{
  /* Local variables */
  int iret;
  char nomvar[] = "ICV";
  float tps=0.;
  int i, ii;

  /* Send icv */
  i = ntcast;
  iret = cp_len(component,
                CP_ITERATION,
                &(tps),
                &(tps),
                &(i),
                nomvar,
                1,
                &(ii),
                icv);
  return iret;
}

/*----------------------------------------------------------------------------
 * Send convergence indicator to Code_Aster
 *----------------------------------------------------------------------------*/

int
send_icv2(void *component,
          int icv)
{
  /* Local variables */
  int iret;
  char nomvar[] = "ICVAST";

  printf("in send_icv2 %d \n",icv);

  /* Send icvext */
  iret = cp_een(component, CP_ITERATION, 0.0, ntcast, nomvar, 1, &(icv));

  printf("iret = %d\n", iret);

  return iret;
}

/*----------------------------------------------------------------------------
 * Initialize communication with Calcium
 *----------------------------------------------------------------------------*/

int
inicom(void *component)
{
  char instance[200];

  int iret = cp_cd(component, instance);

  return iret;
}

/*----------------------------------------------------------------------------
 * End communication with Calcium and stop calculation
 *----------------------------------------------------------------------------*/

int
calfin(void *component)
{
  int iret = cp_fin(component, CP_ARRET);

  exit(0);

  return iret;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif
