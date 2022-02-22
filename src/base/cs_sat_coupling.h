#ifndef __CS_SAT_COUPLING_H__
#define __CS_SAT_COUPLING_H__

/*============================================================================
 * Functions associated with code coupling.
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 * Note: if the context pointer is non-NULL, it must point to valid data
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
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get number of code couplings
 *
 * Fortran interface:
 *
 * SUBROUTINE NBCCPL
 * *****************
 *
 * INTEGER          NBRCPL         : <-- : number of code couplings
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbccpl, NBCCPL)
(
 int  *nbrcpl
);

/*----------------------------------------------------------------------------
 * Set the list of cells and boundary faces associated to a coupling
 * and a cloud of point.
 *
 * The local "support" cells and boundary faces are used to localize
 * the values in the distant "coupled" cells and faces.
 * Depending on the role of sender and/or receiver of the current process
 * in the coupling, some of these sets can be empty or not.
 *
 * The cell values are always localized and interpolated on the distant
 * "cells" support. The face values are localized and interpolated on
 * the distant "face" support if present, or on the distant "cell" support
 * if not.
 *
 * If the input arrays LCESUP and LFBSUP are not ordered, they will be
 * orderd in output.
 *
 * Fortran interface:
 *
 * SUBROUTINE DEFLOC
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 *----------------------------------------------------------------------------*/

void CS_PROCF (defloc, DEFLOC)
(
 const int  *numcpl
);

/*----------------------------------------------------------------------------
 * Get the number of cells and boundary faces, "support", coupled and not
 * localized associated to a given coupling
 *
 * Fortran interface:
 *
 * SUBROUTINE NBECPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCESUP         : <-- : number of "support" cells
 * INTEGER          NFBSUP         : <-- : number of "support" boundary faces
 * INTEGER          NCECPL         : <-- : number of coupled cells
 * INTEGER          NFBCPL         : <-- : number of coupled boundary faces
 * INTEGER          NCENCP         : <-- : number of not coupled cells
 *                                 :     : (since not localized)
 * INTEGER          NFBNCP         : <-- : number of not coupled boundary faces
 *                                 :     : (since not localized)
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbecpl, NBECPL)
(
 const int        *numcpl,
       cs_lnum_t  *ncesup,
       cs_lnum_t  *nfbsup,
       cs_lnum_t  *ncecpl,
       cs_lnum_t  *nfbcpl,
       cs_lnum_t  *ncencp,
       cs_lnum_t  *nfbncp
);

/*----------------------------------------------------------------------------
 * Get the lists of coupled cells and boundary faces (i.e. receiving)
 * associated to a given coupling
 *
 * The number of cells and boundary faces, got with NBECPL(), are used
 * for arguments coherency checks.
 *
 * Fortran interface:
 *
 * SUBROUTINE LELCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCECPL         : --> : number of coupled cells
 * INTEGER          NFBCPL         : --> : number of coupled boundary faces
 * INTEGER          LCECPL(*)      : <-- : list of coupled cells
 * INTEGER          LFBCPL(*)      : <-- : list of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (lelcpl, LELCPL)
(
 const int        *numcpl,
 const cs_lnum_t  *ncecpl,
 const cs_lnum_t  *nfbcpl,
       cs_lnum_t  *lcecpl,
       cs_lnum_t  *lfbcpl
);

/*----------------------------------------------------------------------------
 * Get the lists of not coupled cells and boundary faces (i.e. receiving but
 * not localized) associated to a given coupling
 *
 * The number of cells and boundary faces, got with NBECPL(), are used
 * for arguments coherency checks.
 *
 * Fortran interface:
 *
 * SUBROUTINE LENCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCENCP         : --> : number of not coupled cells
 * INTEGER          NFBNCP         : --> : number of not coupled boundary faces
 * INTEGER          LCENCP(*)      : <-- : list of not coupled cells
 * INTEGER          LFBNCP(*)      : <-- : list of not coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (lencpl, LENCPL)
(
 const int        *numcpl,
 const cs_lnum_t  *ncencp,
 const cs_lnum_t  *nfbncp,
       cs_lnum_t  *lcencp,
       cs_lnum_t  *lfbncp
);

/*----------------------------------------------------------------------------
 * Get the number of distant point associated to a given coupling
 * and localized on the local domain
 *
 * Fortran interface:
 *
 * SUBROUTINE NPDCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCEDIS         : <-- : number of distant cells
 * INTEGER          NFBDIS         : <-- : numbre de distant boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (npdcpl, NPDCPL)
(
 const int        *numcpl,
       cs_lnum_t  *ncedis,
       cs_lnum_t  *nfbdis
);

/*----------------------------------------------------------------------------
 * Get the distant points coordinates associated to a given coupling
 * and a list of points, and the elements number and type (cell or face)
 * "containing" this points.
 *
 * The number of distant points NBRPTS must be equal to one the arguments
 * NCEDIS or NFBDIS given by NPDCPL(), and is given here for coherency checks
 * between the arguments NUMCPL and ITYSUP.
 *
 * Fortran interface:
 *
 * SUBROUTINE COOCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRPTS         : --> : number of distant points
 * INTEGER          ITYDIS         : --> : 1 : access to the points associated
 *                                 :     :     to the distant cells
 *                                 :     : 2 : access to the points associated
 *                                 :     :     to the distant boundary faces
 * INTEGER          ITYLOC         : <-- : 1 : localization on the local cells
 *                                 :     : 2 : localization on the local faces
 * INTEGER          LOCPTS(*)      : <-- : "containing" number associated to
 *                                 :     :   each point
 * DOUBLE PRECISION COOPTS(3,*)    : <-- : distant point coordinates
 * DOUBLE PRECISION DJPPTS(3,*)    : <-- : distant vectors to the coupled face
 * DOUBLE PRECISION PNDPTS(*)      : <-- : distant weighting coefficients
 *----------------------------------------------------------------------------*/

void CS_PROCF (coocpl, COOCPL)
(
 const int        *numcpl,
 const cs_lnum_t  *nbrpts,
 const int        *itydis,
       int        *ityloc,
       cs_lnum_t  *locpts,
       cs_real_t  *coopts,
       cs_real_t  *djppts,
       cs_real_t  *dofpts,
       cs_real_t  *pndpts
);

/*----------------------------------------------------------------------------
 * Get the weighting coefficient needed for a centered-like interpolation
 * in the case of a coupling on boundary faces.
 *
 * Fortran interface:
 *
 * SUBROUTINE PONDCP
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRPTS         : --> : number of distant points
 * INTEGER          ITYLOC         : <-- : 1 : localization on the local cells
 *                                 :     : 2 : localization on the local faces
 * DOUBLE PRECISION PNDCPL(*)      : <-- : weighting coefficients
 *----------------------------------------------------------------------------*/

void CS_PROCF (pondcp, PONDCP)
(
 const int        *numcpl,
 const cs_lnum_t  *nbrpts,
       int        *ityloc,
       cs_real_t  *pndcpl,
       cs_real_t  *distof
);

/*----------------------------------------------------------------------------
 * Exchange a variable associated to a set of point and a coupling.
 *
 * Fortran interface:
 *
 * SUBROUTINE VARCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRDIS         : --> : number of values to send
 * INTEGER          NBRLOC         : --> : number of values to receive
 * INTEGER          ITYVAR         : --> : 1 : variables defined at cells
 *                                 :     : 2 : variables defined at faces
 * INTEGER          STRIDE         : --> : 1 : for scalars
 *                                 :     : 3 : for vectors
 * DOUBLE PRECISION VARDIS(*)      : --> : distant variable(to send)
 * DOUBLE PRECISION VARLOC(*)      : <-- : local variable (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (varcpl, VARCPL)
(
 const int        *numcpl,
 const cs_lnum_t  *nbrdis,
 const cs_lnum_t  *nbrloc,
 const int        *ityvar,
 const cs_lnum_t  *stride,
       cs_real_t  *vardis,
       cs_real_t  *varloc
);

/*----------------------------------------------------------------------------
 * Array of integers exchange, associated to a given coupling.
 *
 * It is assumed that the arrays have the same size and the same values on
 * each group of processus (local and distant).
 *
 * Fortran interface:
 *
 * SUBROUTINE TBICPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRDIS         : --> : number of values to send
 * INTEGER          NBRLOC         : --> : number of values to receive
 * INTEGER          TABDIS(*)      : --> : distant values (to send)
 * INTEGER          TABLOC(*)      : <-- : local values (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tbicpl, TBICPL)
(
 const int        *numcpl,
 const cs_lnum_t  *nbrdis,
 const cs_lnum_t  *nbrloc,
       cs_lnum_t  *vardis,
       cs_lnum_t  *varloc
);

/*----------------------------------------------------------------------------
 * Array of reals exchange, associated to a given coupling.
 *
 * It is assumed that the arrays have the same size and the same values on
 * each group of processus (local and distant).
 *
 * Fortran interface:
 *
 * SUBROUTINE TBRCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRDIS         : --> : number of values to send
 * INTEGER          NBRLOC         : --> : number of values to receive
 * DOUBLE PRECISION TABDIS(*)      : --> : distant values (to send)
 * DOUBLE PRECISION TABLOC(*)      : <-- : local values (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tbrcpl, TBRCPL)
(
 const int        *numcpl,
 const cs_lnum_t  *nbrdis,
 const cs_lnum_t  *nbrloc,
       cs_real_t  *vardis,
       cs_real_t  *varloc
);

/*----------------------------------------------------------------------------
 * Compute the maximum value of an integer variable associated to a coupling.
 *
 * It is assumed that the integer value is the same for each group of
 * processus (local and distant).
 *
 * Fortran interface:
 *
 * SUBROUTINE MXICPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          VALDIS         : --> : distant value (to send)
 * INTEGER          VALMAX         : <-- : local maximum (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (mxicpl, MXICPL)
(
 const int        *const numcpl,
       cs_lnum_t  *const vardis,
       cs_lnum_t  *const varmax
);

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
 *   verbosity          <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_add(const char  *face_cpl_sel_c,
                    const char  *cell_cpl_sel_c,
                    const char  *face_loc_sel_c,
                    const char  *cell_loc_sel_c,
                    const char  *sat_name,
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
 * Destroy all couplings
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_all_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COUPLAGE_H__ */
