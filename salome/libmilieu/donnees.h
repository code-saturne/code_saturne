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

#ifndef __DONNEES_H__
#define __DONNEES_H__

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 *  Global variables
 *============================================================================*/

extern int nbpdtm;
extern int nbssit;

extern int isyncp;
extern int ntchr;

extern double dtref;
extern double ttinit;

extern double epsilo;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_nbpdtm
 *
 * Get time step defined in coupling case XML description
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_nbpdtm(long i);

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_nbssit
 *
 * Get number of iterations of implicit coupling defined in coupling
 * case XML description
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_nbssit(long i);

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_isyncp
 *
 * Get "isyncp" defined in coupling case XML description
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_isyncp(long i);

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_ntchr
 *
 * Get "ntchr" defined in coupling case XML description
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_ntchr(long i);

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_dtref
 *
 * Get "dtref" defined in coupling case XML description
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_dtref(double dt);

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_ttinit
 *
 * Get "ttinit" defined in coupling case XML description
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_ttinit(double tt);

/*----------------------------------------------------------------------------
 * Fonction inter_cs_ast_set_epsilo
 *
 * Get "epsilo" defined in coupling case XML description
 *----------------------------------------------------------------------------*/

void
inter_cs_ast_set_epsilo(double eps);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif /* __DONNEES_H__ */
