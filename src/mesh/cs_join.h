#ifndef __CS_JOIN_H__
#define __CS_JOIN_H__

/*============================================================================
 * Structure and function headers handling with joining operation
 *===========================================================================*/

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
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_join_util.h"
#include "cs_mesh.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *===========================================================================*/

/*============================================================================
 * Type definition
 *===========================================================================*/

/*=============================================================================
 * Global variables
 *===========================================================================*/

/*=============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Add a cs_join_t structure to the list of pending joinings.
 *
 * parameters:
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *
 * returns:
 *   number (1 to n) associated with new joining
 *---------------------------------------------------------------------------*/

int
cs_join_add(const char  *sel_criteria,
            float        fraction,
            float        plane,
            int          verbosity,
            int          visualization);

/*----------------------------------------------------------------------------
 * Set advanced parameters for the joining algorithm.
 *
 * parameters:
 *   join_num       <-> joining operation number
 *   mtf            <-- merge tolerance coefficient
 *   pmf            <-- pre-merge factor
 *   tcm            <-- tolerance computation mode
 *   icm            <-- intersection computation mode
 *   max_break      <-- max number of equivalences to break (merge step)
 *   max_sub_faces  <-- max. possible number of sub-faces by splitting a face
 *   tml            <-- tree max level
 *   tmb            <-- tree max boxes
 *   tmr            <-- tree max ratio
 *   tmr_distrib    <-- tree max ratio for distribution
 *---------------------------------------------------------------------------*/

void
cs_join_set_advanced_param(int      join_num,
                           double   mtf,
                           double   pmf,
                           int      tcm,
                           int      icm,
                           int      max_break,
                           int      max_sub_faces,
                           int      tml,
                           int      tmb,
                           double   tmr,
                           double   tmr_distrib);

/*----------------------------------------------------------------------------
 * Apply all the defined joining operations.
 *
 * parameters:
 *   preprocess <-- true if we are in the preprocessing stage
 *---------------------------------------------------------------------------*/

void
cs_join_all(bool  preprocess);

/*----------------------------------------------------------------------------
 * Clear remaining memory for defined joining operations.
 *---------------------------------------------------------------------------*/

void
cs_join_finalize(void);

/*----------------------------------------------------------------------------
 * Flag boundary faces that will be selected for joining.
 *
 * parameters:
 *   mesh          <-- pointer to mesh structure
 *   preprocess    <-- true if we are in the preprocessing stage
 *   b_select_flag <-> false for boundary faces not selected for joining,
 *                     true for those selected
 *---------------------------------------------------------------------------*/

void
cs_join_mark_selected_faces(const cs_mesh_t  *mesh,
                            bool              preprocess,
                            bool              b_select_flag[]);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_H__ */
