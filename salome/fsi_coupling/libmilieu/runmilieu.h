#ifndef __RUNMILIEU_H__
#define __RUNMILIEU_H__

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

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 *  Global variables
 *============================================================================*/

extern int    nb_dyn;
extern int    nb_for;
extern int    ntcast;
extern double lref;

extern double *xast;
extern double *xvast;
extern double *xvasa;
extern double *xastp;

extern double *foras;
extern double *foaas;
extern double *fopas;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction runmilieu
 *----------------------------------------------------------------------------*/

void
runmilieu(void *icompo);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif /* __RUNMILIEU_H__ */
