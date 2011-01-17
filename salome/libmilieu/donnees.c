/*============================================================================
 *
 *     This file is part of the Code_Saturne CFD tool.
 *
 *     Copyright (C) 2008-2011 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne CFD tool is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne CFD tool is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "runmilieu.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "donnees.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 *  Global variables
 *============================================================================*/

int  nbpdtm = 0;
int  nbssit = 0;

int  isyncp = 0;
int  ntchr = 0;

double  dtref = 0.0;
double  ttinit = 0.0;

double  epsilo = 0.0;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_nbpdtm
 *
 * Récuperation du pas de temps défini dans le XML du cas de couplage
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_nbpdtm(long i)
{
  printf("\n i = %d \n", (int)i);
  nbpdtm = i;
  printf("\n nbpdtm = %d \n", nbpdtm);
}

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_nbssit
 *
 * Récuperation du nombre de sous itérations du couplage implicite défini
 * dans le XML du cas de couplage
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_nbssit(long i)
{
  nbssit = i;
}

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_isyncp
 *
 * Récuperation de isyncp défini dans le XML du cas de couplage
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_isyncp(long i)
{
  isyncp = i;
}

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_ntchr
 *
 * Récuperation de ntchr défini dans le XML du cas de couplage
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_ntchr(long i)
{
  ntchr = i;
}

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_dtref
 *
 * Récuperation de dtref défini dans le XML du cas de couplage
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_dtref(double dt)
{
  dtref = dt;
}

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_ttinit
 *
 * Récuperation de ttinit défini dans le XML du cas de couplage
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_ttinit(double tt)
{
  ttinit = tt;
}

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_epsilo
 *
 * Récuperation de epsilo défini dans le XML du cas de couplage
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_epsilo(double eps)
{
  epsilo = eps;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif
