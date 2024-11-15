#ifndef __CS_SAT_COUPLING_H__
#define __CS_SAT_COUPLING_H__

/*============================================================================
 * Functions associated with code coupling.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Global variables
 *============================================================================*/

extern int  cs_glob_sat_coupling_face_interpolation_type;

/*=============================================================================
 * Type Definitions
 *============================================================================*/

typedef struct _cs_sat_coupling_t cs_sat_coupling_t;

/*----------------------------------------------------------------------------
 * Function pointer to mesh tagging function.
 *
 * Each function of this sort may be used to tag a mesh and associated
 * points for mocatin exclusion.
 *
 * Note: if the context pointer is non-null, it must point to valid data
 * when the selection function is called, so that value or structure
 * should not be temporary (i.e. local);
 *
 * parameters:
 *   context         <-> pointer to optional (untyped) value or structure.
 *   mesh            <-> nodal mesh which should be tagged
 *   n_points        <-- number of points to tag
 *   point_list_base <-- base numbering for point_list
 *   point_list      <-- optional indirection for points
 *   point_tag       --> point tag values (size: n_tags)
 *----------------------------------------------------------------------------*/

typedef void
(cs_sat_coupling_tag_t) (void            *context,
                         fvm_nodal_t     *mesh,
                         cs_lnum_t        n_points,
                         cs_lnum_t        point_list_base,
                         const cs_lnum_t  point_list[],
                         int             *point_tag);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define new code_saturne coupling.
 *
 * The arguments to \ref cs_sat_coupling_define are:
 * \param[in] saturne_name          matching code_saturne application name
 * \param[in] boundary_cpl_criteria boundary face selection criteria for coupled
 *                                  faces, or NULL
 * \param[in] volume_cpl_criteria   cell selection criteria for coupled cells, or
                                    NULL
 * \param[in] boundary_loc_criteria boundary face selection criteria for location
 *                                  (not functional)
 * \param[in] volume_loc_criteria   cell selection criteria for location
 * \param[in] reverse               reverse mode if 1
 * \param[in] verbosity             verbosity level
 *
 * In the case of only 2 code_saturne instances, the 'saturne_name' argument
 * is ignored, as there is only one matching possibility.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * code_saturne instances based on the 'saturne_name' argument.
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_define(const char  *saturne_name,
                       const char  *boundary_cpl_criteria,
                       const char  *volume_cpl_criteria,
                       const char  *boundary_loc_criteria,
                       const char  *volume_loc_criteria,
                       int          reverse,
                       int          verbosity);

/*----------------------------------------------------------------------------
 * Get number of code_saturne couplings.
 *
 * returns:
 *   number of code_saturne couplings
 *----------------------------------------------------------------------------*/

int
cs_sat_coupling_n_couplings(void);

/*----------------------------------------------------------------------------
 * Get pointer to code_saturne coupling.
 *
 * parameters:
 *   coupling_id <-- Id (0 to n-1) of code_saturne coupling
 *
 * returns:
 *   pointer to code_saturne coupling structure
 *----------------------------------------------------------------------------*/

cs_sat_coupling_t *
cs_sat_coupling_by_id(int coupling_id);

/*----------------------------------------------------------------------------
 * Create a sat_coupling_t structure.
 *
 * parameters:
 *   ref_axis           <-- reference axis
 *   face_sel_criterion <-- criterion for selection of boundary faces
 *   cell_sel_criterion <-- criterion for selection of cells
 *   sat_name           <-- code_saturne application name
 *   reverse            <-- reverse mode if 1
 *   verbosity          <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_add(const char  *face_cpl_sel_c,
                    const char  *cell_cpl_sel_c,
                    const char  *face_loc_sel_c,
                    const char  *cell_loc_sel_c,
                    const char  *sat_name,
                    int          reverse,
                    int          verbosity);

/*----------------------------------------------------------------------------
 * Create a new internal code_saturne coupling.
 *
 * arguments:
 *   tag_func          <-- pointer to tagging function
 *   tag_context       <-- pointer to tagging function context
 *   boundary_criteria <-- boundary face selection criteria, or NULL
 *   volume_criteria   <-- volume cell selection criteria, or NULL
 *   loc_tolerance     <-- location tolerance factor (0.1 recommended)
 *   reverse           <-- reverse mode if 1
 *   verbosity         <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_add_internal(cs_sat_coupling_tag_t  *tag_func,
                             void                   *tag_context,
                             const char             *boundary_cpl_criteria,
                             const char             *volume_cpl_criteria,
                             const char             *boundary_loc_criteria,
                             const char             *volume_loc_criteria,
                             float                   loc_tolerance,
                             int                     reverse,
                             int                     verbosity);

/*----------------------------------------------------------------------------
 * Initialize code_saturne couplings.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_all_init(void);

/*----------------------------------------------------------------------------
 * Array of reals exchange, associated to a given coupling.
 *
 * It is assumed that the arrays have the same size and the same values on
 * each group of processes (local and distant).
 *
 * int          cpl_id      : --> : coupling id (0-based)
 * int          nbrdis      : --> : number of values to send
 * int          nbrloc      : --> : number of values to receive
 * cs_real_t    vardis      : --> : distant values (to send)
 * cs_real_t    varloc      : <-- : local values (to receive)
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_array_exchange(int         cpl_id,
                               cs_lnum_t   nbrdis,
                               cs_lnum_t   nbrloc,
                               cs_real_t  *vardis,
                               cs_real_t  *varloc);

/*----------------------------------------------------------------------------*/
/*
 * code_saturne/code_saturne coupling using volumic source terms.
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_exchange_at_cells
(
  cs_field_t *f,    /*!<[in] pointer to cs_field_t */
  cs_real_t  *rhs,  /*!<[out] Explicit terms (RHS) */
  cs_real_t  *fimp  /*!<[out] Implicit source terms */
);

/*----------------------------------------------------------------------------*/
/*
 * code_saturne/code_saturne boundary coupling initialization call
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_bnd_initialize
(
  int *bc_type /*!<[in] boundary face types */
);

/*----------------------------------------------------------------------------*/
/*
 * Initialization of main variables for code_saturne/code_saturne coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_initialize
(
  void
);

/*----------------------------------------------------------------------------*/
/*
 * Set the list of cells and boundary faces associated to a coupling
 * and a cloud of point.
 *
 * The local "support" cells and boundary faces are used to locate
 * the values in the distant "coupled" cells and faces.
 * Depending on the role of sender and/or receiver of the current process
 * in the coupling, some of these sets can be empty or not.
 *
 * The cell values are always located and interpolated on the distant
 * "cells" support. The face values are located and interpolated on
 * the distant "face" support if present, or on the distant "cell" support
 * if not.
 *
 * If the input arrays LCESUP and LFBSUP are not ordered, they will be
 * orderd in output.
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_locate_all
(
 void
);

/*----------------------------------------------------------------------------*/
/*
 * code_saturne/code_saturne coupling using boundary conditions
 */
/*----------------------------------------------------------------------------*/

void
cs_sat_coupling_exchange_at_bnd_faces
(
  int       *bc_type, /*!<[in] boundary face types */
  cs_real_t *dt       /*!<[in] time step (per cell) */
);

/*----------------------------------------------------------------------------
 * Destroy all couplings
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_all_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COUPLAGE_H__ */
