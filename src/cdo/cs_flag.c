/*============================================================================
 * Functions to handle the definition and usage of material properties
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
const cs_flag_t  cs_flag_primal_edge = CS_FLAG_PRIMAL | CS_FLAG_EDGE;
const cs_flag_t  cs_flag_primal_face = CS_FLAG_PRIMAL | CS_FLAG_FACE;
const cs_flag_t  cs_flag_primal_cell = CS_FLAG_PRIMAL | CS_FLAG_CELL;

const cs_flag_t  cs_flag_dual_vtx  = CS_FLAG_DUAL | CS_FLAG_VERTEX;
const cs_flag_t  cs_flag_dual_edge = CS_FLAG_DUAL | CS_FLAG_EDGE;
const cs_flag_t  cs_flag_dual_face = CS_FLAG_DUAL | CS_FLAG_FACE;
const cs_flag_t  cs_flag_dual_cell = CS_FLAG_DUAL | CS_FLAG_CELL;
const cs_flag_t  cs_flag_dual_cell_byc =
  CS_FLAG_DUAL | CS_FLAG_CELL | CS_FLAG_BY_CELL;;
const cs_flag_t  cs_flag_dual_face_byc =
  CS_FLAG_DUAL | CS_FLAG_FACE | CS_FLAG_BY_CELL;
const cs_flag_t  cs_flag_dual_closure_byf =
  CS_FLAG_DUAL | CS_FLAG_CELL | CS_FLAG_BORDER | CS_FLAG_BY_FACE;

/* Additional flags with a more consistent naming when used with FV schemes for
   which the notion of primal/dual is not used */

const cs_flag_t  cs_flag_vertex = CS_FLAG_PRIMAL | CS_FLAG_VERTEX;
const cs_flag_t  cs_flag_cell = CS_FLAG_PRIMAL | CS_FLAG_CELL;

const cs_flag_t  cs_flag_boundary_face =
  CS_FLAG_PRIMAL | CS_FLAG_FACE | CS_FLAG_BORDER;

/* According to the extended flag defined below one can identify which set of
 * quantities or connectivities have to be built on-the-fly and stored in a
 * local structure possibly owned by each thread and with a cellwise scope
 *
 * Store predefined flags to test if one some specific computations of
 * cell quantities
 */

const cs_eflag_t  cs_flag_need_v =
  CS_FLAG_COMP_PV | CS_FLAG_COMP_PVQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_FV;
const cs_eflag_t  cs_flag_need_e =
  CS_FLAG_COMP_PE | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ | CS_FLAG_COMP_EV |
  CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EF  | CS_FLAG_COMP_SEF;
const cs_eflag_t  cs_flag_need_f =
  CS_FLAG_COMP_PF  | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_FE  |
  CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EF  | CS_FLAG_COMP_SEF | CS_FLAG_COMP_HFQ |
  CS_FLAG_COMP_FV;
const cs_eflag_t  cs_flag_need_fe =
  CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EF | CS_FLAG_COMP_SEF;
const cs_eflag_t  cs_flag_need_ef =
  CS_FLAG_COMP_EF;
const cs_eflag_t  cs_flag_need_peq =
  CS_FLAG_COMP_PEQ | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_SEF | CS_FLAG_COMP_PEC;
const cs_eflag_t  cs_flag_need_dfq =
  CS_FLAG_COMP_DFQ | CS_FLAG_COMP_SEF | CS_FLAG_COMP_PEC;
const cs_eflag_t  cs_flag_need_pfq =
  CS_FLAG_COMP_PFQ | CS_FLAG_COMP_HFQ | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_SEF |
  CS_FLAG_COMP_PFC;
const cs_eflag_t  cs_flag_need_deq =
  CS_FLAG_COMP_HFQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_SEF;
const cs_eflag_t  cs_flag_need_pfc =
  CS_FLAG_COMP_PFC | CS_FLAG_COMP_HFQ;

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
  else if (cs_flag_test(loc, cs_flag_primal_edge))
    return "edges";
  else if (cs_flag_test(loc, cs_flag_primal_face))
    return "faces";
  else if (cs_flag_test(loc, cs_flag_boundary_face))
    return "boundary faces";
  else if (cs_flag_test(loc, cs_flag_primal_cell))
    return "cells";
  else if (cs_flag_test(loc, cs_flag_dual_vtx))
    return "dual vertices";
  else if (cs_flag_test(loc, cs_flag_dual_edge))
    return "dual edges";
  else if (cs_flag_test(loc, cs_flag_dual_face))
    return "dual faces";
  else if (cs_flag_test(loc, cs_flag_dual_cell))
    return "dual cells";
  else if (cs_flag_test(loc, cs_flag_dual_face_byc))
    return "dual faces (cellwise)";
  else if (cs_flag_test(loc, cs_flag_dual_closure_byf))
    return "dual cell closure (facewise)";
  else
    return "unknown";
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
