#ifndef __CS_PARAMETERS_H__
#define __CS_PARAMETERS_H__

/*============================================================================
 * General parameters management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure of variable calculation options
 *----------------------------------------------------------------------------*/

typedef struct {
  int     iwarni;
  int     iconv ;
  int     istat ;
  int     idiff ;
  int     idifft;
  int     idften;
  int     iswdyn;
  int     ischcv;
  int     isstpc;
  int     nswrgr;
  int     nswrsm;
  int     imligr;
  int     ircflu;
  int     inc   ;
  double  thetav;
  double  blencv;
  double  epsilo;
  double  epsrsm;
  double  epsrgr;
  double  climgr;
  double  extrag;
  double  relaxv;
} cs_var_cal_opt_t;

/*----------------------------------------------------------------------------
 * Boundary condition types
 *----------------------------------------------------------------------------*/

typedef enum {
  CS_INDEF = 1,
  CS_INLET = 2,
  CS_OUTLET = 3,
  CS_iSYMMETRY = 4,
  CS_SMOOTHWALL = 5,
  CS_ROUGHWALL = 6,
  CS_ESICF = 7,
  CS_SSPCF = 8,
  CS_SOPCF = 9,
  CS_ERUCF = 10,
  CS_EPHCF = 11,
  CS_EQHCF = 12,
  CS_COUPLED = 13
};

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define general field keys.
 *
 * A recommened practice for different submodules would be to use
 * "cs_<module>_key_init() functions to define keys specific to those modules.
 *----------------------------------------------------------------------------*/

void
cs_parameters_define_field_keys(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAMETERS_H__ */
