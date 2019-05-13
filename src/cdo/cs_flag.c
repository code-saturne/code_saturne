/*============================================================================
 * Routines to handle the definition and usage of material properties
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_flag.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/* Default locations */
const cs_flag_t  cs_flag_primal_vtx  = CS_FLAG_PRIMAL | CS_FLAG_VERTEX;
const cs_flag_t  cs_flag_primal_face = CS_FLAG_PRIMAL | CS_FLAG_FACE;
const cs_flag_t  cs_flag_primal_cell = CS_FLAG_PRIMAL | CS_FLAG_CELL;
const cs_flag_t  cs_flag_dual_vtx  = CS_FLAG_DUAL | CS_FLAG_VERTEX;
const cs_flag_t  cs_flag_dual_face = CS_FLAG_DUAL | CS_FLAG_FACE;
const cs_flag_t  cs_flag_dual_cell = CS_FLAG_DUAL | CS_FLAG_CELL;
const cs_flag_t  cs_flag_dual_face_byc =
  CS_FLAG_DUAL | CS_FLAG_FACE | CS_FLAG_BY_CELL;
const cs_flag_t  cs_flag_dual_closure_byf =
  CS_FLAG_DUAL | CS_FLAG_CELL | CS_FLAG_BORDER | CS_FLAG_BY_FACE;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the label associated to a location flag
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_flag_str_location(cs_flag_t  loc)
{
  if (cs_flag_test(loc, cs_flag_primal_vtx))
    return "vertices";
  else if (cs_flag_test(loc, cs_flag_primal_face))
    return "faces";
  else if (cs_flag_test(loc, cs_flag_primal_cell))
    return "cells";
  else if (cs_flag_test(loc, cs_flag_dual_vtx))
    return "dual vertices";
  else if (cs_flag_test(loc, cs_flag_dual_face))
    return "dual faces";
  else if (cs_flag_test(loc, cs_flag_dual_cell))
    return "dual cells";
  else if (cs_flag_test(loc, cs_flag_dual_face_byc))
    return "dual face (cellwise)";
  else if (cs_flag_test(loc, cs_flag_dual_closure_byf))
    return "dual cell closure (facewise)";
  else
    return "";
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
